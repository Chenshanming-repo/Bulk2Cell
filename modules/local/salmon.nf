process SALMON_INDEX {
    tag "$sample"
    label 'salmon'
    input:
    tuple val(sample), path(gtf), path(transcripts)
    tuple path(genome), path(fai), path(annotation)
    output:
    tuple val(sample), path('salmon_index')
    script:
    """
    cut -f1 '${fai}' > decoys.txt
    cat '${transcripts}' '${genome}' > gentrome.fa
    salmon index -t gentrome.fa -d decoys.txt -i salmon_index --keepDuplicates -p ${task.cpus}
    salmon --version > salmon_index/bulk2cell_salmon_version.txt
    """
    stub:
    'mkdir salmon_index; touch salmon_index/versionInfo.json'
}

process SALMON_QUANT {
    tag "$sample"
    label 'salmon_python'
    publishDir "${params.outdir}/${sample}/salmon", mode: 'copy'
    input:
    tuple val(sample), path(index), path(bam), path(bai), path(barcodes)
    path source, stageAs: 'src'
    output:
    tuple val(sample), path('quant/quant.sf')
    path 'quant/aux_info/meta_info.json', emit: metadata
    path 'pseudobulk.json', emit: selection
    script:
    """
    export PYTHONPATH="\$PWD/src"
    extract_pseudobulk.py '${bam}' '${barcodes}' pseudobulk.fastq pseudobulk.json
    salmon quant -i '${index}' -l ${params.salmon_libtype} -r pseudobulk.fastq \
      --validateMappings --noLengthCorrection \
      --fldMean ${params.frag_mean} --fldSD ${params.frag_sd} -p ${task.cpus} -o quant
    """
    stub:
    """
    mkdir -p quant/aux_info
    printf 'Name\tLength\tEffectiveLength\tTPM\tNumReads\n' > quant/quant.sf
    echo '{}' > quant/aux_info/meta_info.json
    echo '{}' > pseudobulk.json
    """
}
