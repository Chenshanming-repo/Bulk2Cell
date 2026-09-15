process VALIDATE_SAMPLES {
    label 'python'
    cache false
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'
    input:
    path samplesheet
    output:
    path 'samples.validated.json'
    script:
    """
    workflow_support.py validate '${samplesheet}' samples.validated.json
    """
}

process PREPARE_REFERENCE {
    label 'python'
    input:
    path genome, stageAs: 'input_genome.fa'
    path annotation, stageAs: 'input_annotation.gtf'
    output:
    tuple path('reference.fa'), path('reference.fa.fai'), path('annotation.gtf')
    script:
    """
    cp --reflink=auto '${genome}' reference.fa
    cp --reflink=auto '${annotation}' annotation.gtf
    samtools faidx reference.fa
    """
    stub:
    """
    touch reference.fa reference.fa.fai annotation.gtf
    """
}

process BUILD_UNION {
    tag "$sample"
    label 'python'
    publishDir "${params.outdir}/${sample}/reference", mode: 'copy'
    input:
    tuple val(sample), path(gff), path(classification)
    tuple path(genome), path(fai), path(annotation)
    path source, stageAs: 'src'
    output:
    tuple val(sample), path('union.gtf'), path('union.fa'), emit: union
    path 'union.json', emit: provenance
    script:
    """
    export PYTHONPATH="\$PWD/src"
    workflow_support.py union --reference '${annotation}' --gff '${gff}' \
      --classification '${classification}' --genome '${genome}' --prefix union
    """
    stub:
    """
    touch union.gtf union.fa
    echo '{}' > union.json
    """
}
