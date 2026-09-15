process CELLRANGER_MKREF {
    tag "$sample"
    label 'cellranger'
    publishDir "${params.outdir}/${sample}/cellranger_reference", mode: 'copy'
    input:
    val sample
    tuple path(genome), path(fai), path(annotation)
    output:
    tuple val(sample), path('cr_reference')
    script:
    """
    ${params.cellranger} mkref --genome=cr_reference --fasta='${genome}' --genes='${annotation}' \
      --nthreads=${task.cpus} --memgb=${task.memory.toGiga()}
    ${params.cellranger} --version > cr_reference/bulk2cell_cellranger_version.txt
    """
    stub:
    'mkdir cr_reference; touch cr_reference/reference.json'
}

process CELLRANGER_COUNT {
    tag "$sample"
    label 'cellranger'
    publishDir "${params.outdir}/${sample}/cellranger", mode: 'copy'
    input:
    tuple val(sample), path(fastqs), val(fastq_sample), val(chemistry), path(reference)
    output:
    tuple val(sample), path('count/outs/possorted_genome_bam.bam'), path('count/outs/possorted_genome_bam.bam.bai'), path('count/outs/filtered_feature_bc_matrix/barcodes.tsv.gz')
    path 'count/outs/web_summary.html', emit: report
    path 'count/outs/filtered_feature_bc_matrix.h5', emit: gene_matrix
    script:
    """
    ${params.cellranger} count --id=count --transcriptome='${reference}' --fastqs='${fastqs}' \
      --sample='${fastq_sample}' --chemistry='${chemistry}' --create-bam=true \
      --include-introns=false --localcores=${task.cpus} --localmem=${task.memory.toGiga()}
    """
    stub:
    """
    mkdir -p count/outs/filtered_feature_bc_matrix
    touch count/outs/possorted_genome_bam.bam count/outs/possorted_genome_bam.bam.bai
    touch count/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
    touch count/outs/web_summary.html count/outs/filtered_feature_bc_matrix.h5
    """
}
