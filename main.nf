nextflow.enable.dsl=2

include { VALIDATE_SAMPLES; PREPARE_REFERENCE; BUILD_UNION } from './modules/local/reference'
include { ISOSEQ_REFINE; ISOSEQ_CLUSTER; ISOSEQ_ALIGN; ISOSEQ_COLLAPSE; PIGEON_FILTER } from './modules/local/isoseq'
include { CELLRANGER_MKREF; CELLRANGER_COUNT } from './modules/local/cellranger'
include { SALMON_INDEX; SALMON_QUANT } from './modules/local/salmon'
include { BULK2CELL_QUANTIFY } from './modules/local/bulk2cell'

workflow {
    if (params.help) {
        log.info 'bulk2cell: nextflow run main.nf -profile conda -params-file examples/params.yaml\nRequired: --input samples.csv --genome genome.fa --annotation genes.gtf\nSee README.md for PacBio stages, primers, 10x 3-prime chemistry and Cell Ranger setup.'
    } else {
        ['input','genome','annotation'].each { key ->
            if (!params[key]) error "Missing required parameter --${key}; see README.md"
        }
        if (params.genome.toString().endsWith('.gz') || params.annotation.toString().endsWith('.gz')) error 'Provide uncompressed genome FASTA and annotation GTF'
        if (!(params.salmon_libtype in ['A','U','SF','SR'])) error 'salmon_libtype must be A, U, SF or SR'
        if (params.tau < 0 || params.frag_mean <= 0 || params.frag_sd <= 0) error 'Invalid abundance or fragment-length parameters'
        if (params.cellranger.toString().find(/[\s'"`$]/)) error 'cellranger must be an executable name or path without shell syntax or spaces'
        def checkedFile = { value ->
            if (value.toString().find(/['"`$\n\r]/)) error 'Input paths must not contain quotes or shell syntax'
            file(value, checkIfExists: true)
        }
        def source = Channel.value(file("${projectDir}/software/bulk2cell/src"))
        VALIDATE_SAMPLES(Channel.value(checkedFile(params.input)))
        def samples = VALIDATE_SAMPLES.out.flatMap { json -> new groovy.json.JsonSlurperClassic().parseText(json.text) }
        def bulk = samples.map { row ->
            tuple(row.sample, row.pacbio_stage, row.pacbio_bams.collect { file(it, checkIfExists:true) }, row.primers ? [file(row.primers,checkIfExists:true)] : [])
        }
        def shortreads = samples.map { row -> tuple(row.sample, file(row.fastq_dir,checkIfExists:true), row.fastq_sample, row.chemistry) }
        PREPARE_REFERENCE(Channel.value(checkedFile(params.genome)), Channel.value(checkedFile(params.annotation)))
        def reference = PREPARE_REFERENCE.out
        ISOSEQ_REFINE(bulk)
        ISOSEQ_CLUSTER(ISOSEQ_REFINE.out)
        ISOSEQ_ALIGN(ISOSEQ_CLUSTER.out, reference)
        ISOSEQ_COLLAPSE(ISOSEQ_ALIGN.out.join(ISOSEQ_REFINE.out))
        PIGEON_FILTER(ISOSEQ_COLLAPSE.out[0], reference)
        BUILD_UNION(PIGEON_FILTER.out, reference, source)
        // Cell Ranger uses the supplied genome/annotation independently of Iso-Seq.
        CELLRANGER_MKREF(shortreads.map { sample, fastqs, prefix, chemistry -> sample }, reference)
        CELLRANGER_COUNT(shortreads.join(CELLRANGER_MKREF.out))
        SALMON_INDEX(BUILD_UNION.out.union, reference)
        SALMON_QUANT(SALMON_INDEX.out.join(CELLRANGER_COUNT.out[0]), source)
        BULK2CELL_QUANTIFY(PIGEON_FILTER.out.join(CELLRANGER_COUNT.out[0]).join(SALMON_QUANT.out[0]), reference, source)
    }
}
