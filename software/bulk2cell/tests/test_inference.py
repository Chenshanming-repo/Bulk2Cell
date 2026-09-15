"""Analytical tests for conservation, compatibility and prior sensitivity."""
import numpy as np
import pytest
from bulk2cell.models import Transcript, Molecule
from bulk2cell.inference import hybrid_prior, compatible, distance, structural_groups, fit_distance, likelihoods, quantify_gene


def tx(name='a', strand='+', exons=((100,200),(300,400)), count=0):
    """Construct a small two-exon transcript."""
    return Transcript(name, 'g', 'chr1', strand, exons, count)


def test_hybrid_limits_and_validation():
    """LR counts and normalized Salmon abundance obey the specified formula."""
    assert np.allclose(hybrid_prior([9,1],[1,9],10), [.5,.5])
    assert np.allclose(hybrid_prior([0,0],[1,3],10), [.25,.75])
    assert np.allclose(hybrid_prior([0,0],[0,0],0), [.5,.5])
    with pytest.raises(ValueError): hybrid_prior([-1,2],[1,1],10)
    with pytest.raises(ValueError): hybrid_prior([1,2],[1,1],-1)


def test_spliced_coordinates_and_compatibility():
    """Introns contribute no transcript distance and incompatible junctions fail."""
    m = Molecule('c','u','g',((180,200),(300,320)),((200,300),))
    assert compatible(m,tx())
    assert distance(m,tx()) == 80
    assert distance(m,tx(strand='-')) == 80
    assert not compatible(m,tx(exons=((100,210),(300,400))))


def test_grouping_retains_upstream_distinctions_outside_window():
    """Only identical spliced terminal windows share a structural group."""
    a = tx()
    b = tx('b',exons=((50,80),(300,400)))
    assert structural_groups([a,b],80) == [[0,1]]
    assert structural_groups([a,b],150) == [[0],[1]]


def test_kde_positive_and_tes_marginalization():
    """A two-point TES distribution averages the two endpoint likelihoods."""
    model = fit_distance([10,20,30],bandwidth=10)
    m = Molecule('c','u','g',((350,370),),())
    base = likelihoods([m],[tx()],model,{})[0,0]
    mixed = likelihoods([m],[tx()],model,{'a':[(400,1),(410,1)]})[0,0]
    assert mixed == pytest.approx((base + model.pdf(40))/2)
    assert model.pdf(-1) == 0


def test_group_count_conservation_and_prior_decomposition():
    """Indistinguishable isoforms retain the hybrid ratio and molecule mass."""
    ts = [tx('a',count=6),tx('b',count=3)]
    ms = [Molecule('c',str(i),'g',((350,370),),()) for i in range(10)]
    result = quantify_gene(ms,ts,{}, {},tau=0,window=80)
    assert result['group_counts'].sum() == pytest.approx(10)
    assert result['isoform_counts'][0].tolist() == pytest.approx([20/3,10/3])
    assert result['groups'] == [[0,1]]
    assert result['qc']['assigned_molecules'] == 10


def test_incompatible_molecule_is_not_rescued_by_prior():
    """A high prior cannot assign a molecule outside the annotation."""
    result=quantify_gene([Molecule('c','u','g',((500,520),),())],[tx()],{}, {})
    assert result['isoform_counts'].sum() == 0
    assert result['qc']['incompatible_molecules'] == 1


def test_zero_prior_component_can_be_recovered_by_em():
    """A numerical initialization floor must not trigger premature convergence."""
    from bulk2cell.inference import _em
    counts, converged, iterations = _em(np.tile([1.,2.],(20,1)),np.array([1.,0.]),0,1e-7,500)
    assert counts[1] > 19.99
    assert converged
    assert iterations > 1


def test_distinct_tes_likelihoods_refine_terminal_groups():
    """LR TES information capable of separating members must survive grouping."""
    ts=[tx('a'),tx('b',exons=((50,80),(300,400)))]
    ms=[Molecule('c','u','g',((350,370),),())]
    result=quantify_gene(ms,ts,{}, {'b':[(450,1)]},window=80)
    assert result['groups']==[[0],[1]]


def test_kde_tail_does_not_drop_compatible_molecule():
    """Very small capture density is evidence, not structural incompatibility."""
    t=tx(exons=((0,2000),))
    m=Molecule('c','u','g',((580,600),),())
    result=quantify_gene([m],[t],{}, {},model=fit_distance([100]))
    assert result['qc']['assigned_molecules']==1
    assert result['isoform_counts'].sum()==pytest.approx(1)
