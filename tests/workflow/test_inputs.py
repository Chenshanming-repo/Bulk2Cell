import importlib.util
from pathlib import Path
import csv
import pytest

MODULE = Path(__file__).parents[2] / 'bin' / 'workflow_support.py'

def load():
    assert MODULE.exists(), 'workflow input validator has not been implemented'
    spec = importlib.util.spec_from_file_location('workflow_support', MODULE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

@pytest.fixture
def sheet(tmp_path):
    (tmp_path/'reads.bam').touch()
    (tmp_path/'primers.fa').write_text('>primer_5p\nACGT\n>primer_3p\nTGCA\n')
    fastqs = tmp_path/'fastqs'; fastqs.mkdir()
    for read in (1,2):
        (fastqs/f'library_S1_L001_R{read}_001.fastq.gz').touch()
    row = dict(sample='donor1',pacbio_bams='reads.bam',pacbio_stage='hifi',primers='primers.fa',fastq_dir='fastqs',fastq_sample='library',chemistry='SC3Pv3')
    def write(**changes):
        updated = {**row, **changes}
        path = tmp_path/'samples.csv'
        with path.open('w') as handle:
            writer=csv.DictWriter(handle, fieldnames=row); writer.writeheader(); writer.writerow(updated)
        return path
    return write

def test_resolves_paths_relative_to_samplesheet(sheet):
    path=sheet(); rows=load().validate_samples(path)
    assert rows[0]['pacbio_bams']==[str(path.parent/'reads.bam')]
    assert rows[0]['sample']=='donor1'

@pytest.mark.parametrize('changes,match', [
    ({'pacbio_stage':'subreads'}, 'stage'),
    ({'primers':''}, 'primers'),
    ({'sample':'bad;command'}, 'sample'),
    ({'pacbio_bams':'missing.bam'}, 'not found'),
    ({'chemistry':'SC5P-PE'}, 'chemistry'),
    ({'fastq_sample':'wrong'}, 'FASTQ'),
])
def test_rejects_invalid_inputs(sheet,changes,match):
    with pytest.raises(ValueError,match=match):load().validate_samples(sheet(**changes))

def test_flnc_does_not_require_primers(sheet):
    assert load().validate_samples(sheet(pacbio_stage='flnc',primers=''))[0]['primers']==''

def test_rejects_unpaired_lanes(sheet):
    path=sheet(); (path.parent/'fastqs/library_S1_L001_R2_001.fastq.gz').unlink()
    with pytest.raises(ValueError,match='paired'):load().validate_samples(path)

def test_rejects_duplicate_samples(sheet):
    path=sheet(); lines=path.read_text().splitlines(); path.write_text('\n'.join(lines+[lines[1]])+'\n')
    with pytest.raises(ValueError,match='duplicate'):load().validate_samples(path)


def test_multiple_fastq_prefixes_for_one_library(sheet):
    path = sheet(fastq_sample='HBA8_3-1;HBA8_3-2;HBA8_3-3')
    for prefix, lane in [('HBA8_3-1', '001'), ('HBA8_3-2', '003'), ('HBA8_3-3', '003')]:
        for read in (1, 2):
            (path.parent / 'fastqs' / f'{prefix}_S1_L{lane}_R{read}_001.fastq.gz').touch()
    rows = load().validate_samples(path)
    assert len(rows) == 1
    assert rows[0]['fastq_sample'] == 'HBA8_3-1,HBA8_3-2,HBA8_3-3'


@pytest.mark.parametrize('prefixes,match', [
    ('library;missing', 'no matching.*missing'),
    ('library;library', 'duplicate'),
    ('library;', 'invalid fastq_sample'),
    ('library;bad$(id)', 'invalid fastq_sample'),
])
def test_rejects_invalid_fastq_prefix_lists(sheet, prefixes, match):
    with pytest.raises(ValueError, match=match):
        load().validate_samples(sheet(fastq_sample=prefixes))


def test_requires_pairs_for_each_fastq_prefix(sheet):
    path = sheet(fastq_sample='library;second')
    (path.parent / 'fastqs/second_S1_L001_R1_001.fastq.gz').touch()
    with pytest.raises(ValueError, match='paired'):
        load().validate_samples(path)
