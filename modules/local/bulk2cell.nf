process BULK2CELL_QUANTIFY {
    tag "$sample"
    label 'python_heavy'
    publishDir "${params.outdir}/${sample}/bulk2cell", mode: 'copy'
    input:
    tuple val(sample), path(gff), path(classification), path(bam), path(bai), path(barcodes), path(salmon)
    tuple path(genome), path(fai), path(annotation)
    path source, stageAs: 'src'
    output:
    tuple val(sample), path('quantification')
    script:
    """
    export PYTHONPATH="\$PWD/src"
    export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
    python -m bulk2cell quantify-full --reference '${annotation}' --isoseq-gff '${gff}' \
      --classification '${classification}' --bam '${bam}' --barcodes '${barcodes}' \
      --salmon '${salmon}' --tau ${params.tau} --workers ${task.cpus} \
      --max-likelihood-entries ${params.max_likelihood_entries} --out quantification
    """
    stub:
    """
    mkdir quantification
    echo '{"status":"stub","biological_analysis":false}' > quantification/run.json
    """
}
