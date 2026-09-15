"""Spliced-coordinate likelihoods and prior-regularized isoform-group EM.

One UMI supplies one likelihood: its most 3-prime aligned boundary. Other
reads constrain structural compatibility but are not multiplied as independent
molecules. Structural groups use exact spliced terminal windows, a declared
capture-window assumption rather than proof of universal identifiability.
"""
from dataclasses import dataclass
import numpy as np
from scipy.special import logsumexp


def hybrid_prior(lr, sr, tau=10.0):
    """Combine LR counts with gene-normalized SR abundance using the spec formula.

    Missing SR mass defaults to uniform; when both LR depth and tau are zero,
    return the SR distribution. No pseudocount silently alters LR support.
    """
    lr, sr = np.asarray(lr,dtype=float), np.asarray(sr,dtype=float)
    if lr.ndim != 1 or lr.size == 0 or lr.shape != sr.shape:
        raise ValueError('LR and SR must be matching nonempty vectors')
    if not np.isfinite(tau) or tau < 0 or any(np.any(~np.isfinite(x)) or np.any(x < 0) for x in (lr,sr)):
        raise ValueError('counts and tau must be finite and nonnegative')
    sr = sr/sr.sum() if sr.sum() else np.full(lr.size,1/lr.size)
    return (lr+tau*sr)/(lr.sum()+tau) if lr.sum()+tau else sr


def compatible(molecule, transcript):
    """Require all aligned bases and exact intron boundaries to fit one transcript."""
    if molecule.gene_id != transcript.gene_id or not molecule.blocks:
        return False
    if not all(any(x <= a and b <= y for x,y in transcript.exons) for a,b in molecule.blocks):
        return False
    introns = {(transcript.exons[i][1],transcript.exons[i+1][0]) for i in range(len(transcript.exons)-1)}
    return set(molecule.junctions).issubset(introns)


def distance(molecule, transcript, tes=None):
    """Return spliced distance from the molecule's 3-prime boundary to a TES.

    TES alternatives may only move within/extend the terminal exon. Moving
    across an upstream splice boundary is invalid and has zero likelihood.
    """
    endpoint = max(b for a,b in molecule.blocks) if transcript.strand == '+' else min(a for a,b in molecule.blocks)
    if transcript.strand == '+':
        base = sum(max(0,b-max(a,endpoint)) for a,b in transcript.exons)
        if tes is not None and tes < transcript.exons[-1][0]: return -1
        return base + (0 if tes is None else tes-transcript.tes)
    base = sum(max(0,min(b,endpoint)-a) for a,b in transcript.exons)
    if tes is not None and tes > transcript.exons[0][1]: return -1
    return base + (0 if tes is None else transcript.tes-tes)


def structural_groups(transcripts, window=600):
    """Group identical genomic structures within the last `window` exonic bases.

    Exact endpoints are preserved; unlike a fixed TES tolerance this never
    collapses distinct ends merely because they are nearby. Full structures
    and stable transcript identifiers remain in the membership output.
    """
    if window <= 0: raise ValueError('capture window must be positive')
    groups = {}
    for i,t in enumerate(transcripts):
        remaining, terminal = window, []
        for a,b in (reversed(t.exons) if t.strand == '+' else t.exons):
            length = min(remaining,b-a)
            terminal.append((b-length,b) if t.strand == '+' else (a,a+length))
            remaining -= length
            if remaining == 0: break
        key = (t.gene_id,t.chrom,t.strand,tuple(terminal))
        groups.setdefault(key,[]).append(i)
    return list(groups.values())


@dataclass
class DistanceModel:
    """Reflected Gaussian KDE on nonnegative distances, or exponential fallback."""
    samples: np.ndarray
    bandwidth: float
    fallback_scale: float = 300.0

    def logpdf(self, value):
        """Evaluate log density without Gaussian tail underflow."""
        if value < 0 or not np.isfinite(value): return -np.inf
        if not len(self.samples): return -value/self.fallback_scale-np.log(self.fallback_scale)
        z1=(value-self.samples)/self.bandwidth
        z2=(value+self.samples)/self.bandwidth
        return float(logsumexp(np.concatenate((-z1*z1/2,-z2*z2/2)))-np.log(len(self.samples)*self.bandwidth*np.sqrt(2*np.pi)))

    def pdf(self, value):
        """Evaluate scalar density for plotting; inference uses stable logpdf."""
        return float(np.exp(self.logpdf(value)))


def fit_distance(distances, bandwidth=30.0, max_samples=2000):
    """Fit a reflected KDE using deterministic evenly spaced subsampling.

    Inputs must come from uniquely compatible SR molecules. Empty training
    data selects a reported exponential fallback, not an LR-trained density.
    """
    values = np.asarray(distances,dtype=float)
    if values.ndim != 1 or np.any(~np.isfinite(values)) or np.any(values < 0):
        raise ValueError('training distances must be finite and nonnegative')
    if not np.isfinite(bandwidth) or bandwidth <= 0 or max_samples < 1:
        raise ValueError('KDE bandwidth and sample limit must be positive')
    if len(values)>max_samples:
        values = np.sort(values)[np.linspace(0,len(values)-1,max_samples,dtype=int)]
    return DistanceModel(values,bandwidth)


def likelihoods(molecules, transcripts, model, tes, normalize=False):
    """Marginalize read-to-TES likelihood over weighted LR endpoint observations.

    `tes` maps canonical transcript IDs to (genomic boundary, nonnegative
    weight) pairs. Missing endpoint data uses the annotation TES explicitly.
    """
    result = np.full((len(molecules),len(transcripts)), -np.inf)
    cache = {}
    for j,t in enumerate(transcripts):
        observations = tes.get(t.id,[(t.tes,1.0)])
        if not observations or any(not np.isfinite(s) or s < 0 or not np.isfinite(w) or w < 0 for s,w in observations):
            raise ValueError('TES positions and weights must be finite and nonnegative')
        total = sum(w for s,w in observations)
        if total <= 0: raise ValueError('TES weights must have positive total')
        for i,m in enumerate(molecules):
            if compatible(m,t):
                # Quantized integer distances permit reuse across cells/UMIs.
                terms=[]
                for s,w in observations:
                    if w == 0: continue
                    d = distance(m,t,s)
                    if d not in cache: cache[d] = model.logpdf(d)
                    terms.append(np.log(w/total)+cache[d])
                result[i,j] = logsumexp(terms)
    if normalize and len(molecules):
        maximum=np.max(result,axis=1,keepdims=True)
        valid=np.isfinite(maximum[:,0])
        result[valid] -= maximum[valid]
    return np.exp(result)



def _em(matrix, prior, strength, tolerance, max_iter):
    """Fit a categorical mixture and return expected counts plus convergence QC."""
    # Interior initialization lets evidence recover zero-prior components.
    # Starting at a tiny floor falsely satisfies absolute convergence early.
    abundance = np.full(len(prior), 1 / len(prior))
    converged = False
    for iteration in range(1,max_iter+1):
        posterior = matrix*abundance
        posterior /= posterior.sum(axis=1,keepdims=True)
        counts = posterior.sum(axis=0)
        updated = (counts+strength*prior)/(len(matrix)+strength)
        updated = np.maximum(updated,1e-300); updated /= updated.sum()
        if np.max(np.abs(updated-abundance)) < tolerance:
            abundance = updated; converged=True; break
        abundance = updated
    posterior = matrix*abundance
    posterior /= posterior.sum(axis=1,keepdims=True)
    return posterior.sum(axis=0), converged, iteration


def quantify_gene(molecules, transcripts, sr, tes, tau=10.0, window=600,
                  model=None, em_strength=1.0, tolerance=1e-7, max_iter=500,
                  precomputed_likelihood=None, em_batch=None):
    """Estimate cell/group counts, then decompose groups with fixed hybrid weights.

    Cells are fit separately. Counts exclude prior pseudomolecules and conserve
    the number of structurally compatible UMIs. Grouping is refined whenever
    observed molecule likelihoods separate members of a terminal-window
    group; this prevents upstream evidence from being erased by grouping.
    """
    if not transcripts: raise ValueError('a gene must contain transcripts')
    if len({t.gene_id for t in transcripts}) != 1: raise ValueError('quantify_gene requires one gene')
    if not np.isfinite(em_strength) or em_strength < 0 or not np.isfinite(tolerance) or tolerance <= 0 or max_iter < 1:
        raise ValueError('invalid EM controls')
    prior = hybrid_prior([t.lr_count for t in transcripts],[sr.get(t.id,0) for t in transcripts],tau)
    if model is None:
        training=[]
        for m in molecules:
            hits=[t for t in transcripts if compatible(m,t)]
            if len(hits)==1: training.append(distance(m,hits[0]))
        model=fit_distance(training)
    matrix=(likelihoods(molecules,transcripts,model,tes,normalize=True)
            if precomputed_likelihood is None else np.asarray(precomputed_likelihood))
    if matrix.shape != (len(molecules),len(transcripts)) or np.any(~np.isfinite(matrix)) or np.any(matrix < 0):
        raise ValueError("invalid precomputed likelihood matrix")
    groups=[]
    for group in structural_groups(transcripts,window):
        signatures={}
        for j in group: signatures.setdefault(matrix[:,j].tobytes(),[]).append(j)
        groups.extend(signatures.values())
    weights=[]
    for group in groups:
        p=prior[group]
        weights.append(p/p.sum() if p.sum() else np.full(len(group),1/len(group)))
    group_likelihood=np.column_stack([matrix[:,g]@w for g,w in zip(groups,weights)])
    group_prior=np.array([prior[g].sum() for g in groups])
    cells=sorted({m.cell for m in molecules})
    group_counts=np.zeros((len(cells),len(groups)))
    isoform_counts=np.zeros((len(cells),len(transcripts)))
    valid=group_likelihood.sum(axis=1)>0
    qc={'assigned_molecules':int(valid.sum()),'incompatible_molecules':int((~valid).sum()),
        'em_nonconverged_cells':0,'max_em_iterations':0,'distance_training_samples':len(model.samples),
        'distance_fallback':not bool(len(model.samples))}
    cell_rows={c:[] for c in cells}
    for i,m in enumerate(molecules):
        if valid[i]: cell_rows[m.cell].append(i)
    def cell_solutions():
        """Keep cell fit batches bounded and preserve scalar normalization."""
        batch = []
        indices = []
        for i, c in enumerate(cells):
            rows = cell_rows[c]
            if not rows:
                continue
            values = group_likelihood[rows]
            values = values / values.max(axis=1, keepdims=True)
            if em_batch is None:
                yield i, _em(values, group_prior, em_strength, tolerance, max_iter)
            else:
                indices.append(i)
                batch.append(values)
                if len(batch) == 128:
                    yield from zip(indices, em_batch(batch, group_prior, em_strength, tolerance, max_iter))
                    batch = []; indices = []
        if batch:
            yield from zip(indices, em_batch(batch, group_prior, em_strength, tolerance, max_iter))
    for i, (counts, converged, iterations) in cell_solutions():
        group_counts[i]=counts
        qc['em_nonconverged_cells']+=int(not converged)
        qc['max_em_iterations']=max(qc['max_em_iterations'],iterations)
        for k,(g,w) in enumerate(zip(groups,weights)): isoform_counts[i,g]=counts[k]*w
    return dict(cells=cells,groups=groups,prior=prior,group_counts=group_counts,isoform_counts=isoform_counts,qc=qc)
