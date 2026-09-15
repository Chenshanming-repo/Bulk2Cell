"""Immutable genomic records; all intervals use zero-based half-open coordinates."""
from dataclasses import dataclass
import math


@dataclass(frozen=True)
class Transcript:
    """A full transcript structure with a nonnegative long-read support count."""
    id: str
    gene_id: str
    chrom: str
    strand: str
    exons: tuple[tuple[int, int], ...]
    lr_count: float = 0.0
    source: str = 'reference'

    def __post_init__(self):
        """Reject ambiguous loci, overlapping exons and invalid support at ingestion."""
        if not self.id or not self.gene_id or not self.chrom or self.strand not in ('+', '-'):
            raise ValueError('transcript requires identifiers and a + or - strand')
        if not self.exons or any(a < 0 or b <= a for a,b in self.exons):
            raise ValueError('exons must be nonempty positive-length genomic intervals')
        if any(self.exons[i][1] > self.exons[i+1][0] for i in range(len(self.exons)-1)):
            raise ValueError('exons must be sorted and nonoverlapping')
        if not math.isfinite(self.lr_count) or self.lr_count < 0:
            raise ValueError('LR support must be finite and nonnegative')

    @property
    def tes(self):
        """Return the genomic interbase TES boundary, respecting transcript strand."""
        return self.exons[-1][1] if self.strand == '+' else self.exons[0][0]


@dataclass(frozen=True)
class Molecule:
    """Union of alignment evidence for one corrected cell/UMI within one gene."""
    cell: str
    umi: str
    gene_id: str
    blocks: tuple[tuple[int, int], ...]
    junctions: tuple[tuple[int, int], ...]
