process ISOSEQ_REFINE {
    tag "$sample"
    label 'isoseq'
    publishDir "${params.outdir}/${sample}/isoseq/flnc", mode: 'copy'
    input:
    tuple val(sample), val(stage), path(bams, stageAs: 'bams/part??/*'), path(primers, stageAs: 'primers/*')
    output:
    tuple val(sample), path('flnc.bam'), path('flnc.bam.pbi')
    script:
    def inputs = bams.collect { "'${it}'" }.join(' ')
    def refinement = stage == 'hifi' ? """
        mkdir lima
        lima merged.bam '${primers[0]}' lima/fl.bam --isoseq --peek-guess -j ${task.cpus}
        parts=(lima/*.bam)
        if [ "\${#parts[@]}" -ne 1 ]; then
            echo 'Expected one demultiplexed primer pair; split multiplexed biological samples first' >&2
            exit 1
        fi
        isoseq refine "\${parts[0]}" '${primers[0]}' flnc.bam --require-polya -j ${task.cpus}
    """ : 'mv merged.bam flnc.bam'
    """
    pbmerge -o merged.bam ${inputs}
    pbindex merged.bam
    ${refinement}
    pbindex flnc.bam
    isoseq --version > isoseq.version.txt
    """
    stub:
    """
    touch flnc.bam flnc.bam.pbi
    """
}

process ISOSEQ_CLUSTER {
    tag "$sample"
    label 'isoseq_heavy'
    input:
    tuple val(sample), path(flnc), path(pbi)
    output:
    tuple val(sample), path('clustered.bam')
    script:
    """
    isoseq cluster2 '${flnc}' clustered.bam -j ${task.cpus}
    """
    stub:
    'touch clustered.bam'
}

process ISOSEQ_ALIGN {
    tag "$sample"
    label 'isoseq_heavy'
    input:
    tuple val(sample), path(clustered)
    tuple path(genome), path(fai), path(annotation)
    output:
    tuple val(sample), path('mapped.bam'), path('mapped.bam.bai')
    script:
    """
    pbmm2 align '${genome}' '${clustered}' mapped.bam --preset ISOSEQ --sort -j ${task.cpus}
    test -f mapped.bam.bai || samtools index -@ ${task.cpus} mapped.bam
    """
    stub:
    'touch mapped.bam mapped.bam.bai'
}

process ISOSEQ_COLLAPSE {
    tag "$sample"
    label 'isoseq'
    publishDir "${params.outdir}/${sample}/isoseq/collapse", mode: 'copy'
    input:
    tuple val(sample), path(mapped), path(bai), path(flnc), path(pbi)
    output:
    tuple val(sample), path('collapsed.gff'), path('collapsed.flnc_count.txt')
    path 'collapsed.read_stat.txt', optional: true
    script:
    """
    isoseq collapse --do-not-collapse-extra-5exons -j ${task.cpus} '${mapped}' '${flnc}' collapsed.gff
    """
    stub:
    'touch collapsed.gff collapsed.flnc_count.txt collapsed.read_stat.txt'
}

process PIGEON_FILTER {
    tag "$sample"
    label 'isoseq'
    stageInMode 'copy'
    publishDir "${params.outdir}/${sample}/isoseq/pigeon", mode: 'copy'
    input:
    tuple val(sample), path(gff), path(counts)
    tuple path(genome), path(fai), path(annotation)
    output:
    tuple val(sample), path('collapsed.sorted.filtered_lite.gff'), path('isoforms_classification.filtered_lite_classification.txt')
    script:
    """
    pigeon prepare '${annotation}' '${genome}' '${gff}'
    pigeon classify collapsed.sorted.gff annotation.sorted.gtf '${genome}' --fl '${counts}' -d . -o isoforms -j ${task.cpus}
    pigeon filter isoforms_classification.txt --isoforms collapsed.sorted.gff -j ${task.cpus}
    """
    stub:
    'touch collapsed.sorted.filtered_lite.gff isoforms_classification.filtered_lite_classification.txt'
}
