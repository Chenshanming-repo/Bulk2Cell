#!/usr/bin/env Rscript
# Recreate original unique.reads using globally reconciled single-model genes.
suppressMessages(suppressWarnings(library(data.table)))
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3) stop("usage: batch_unique.R FILTERED_RDS GENES OUTPUT_TSV")
reads <- as.data.table(readRDS(args[[1]]))
genes <- if (file.info(args[[2]])$size == 0) character() else fread(args[[2]], header=FALSE)[[1]]
columns <- c("seqnames.rd", "start.rd", "end.rd", "strand.rd", "dist_END", "frag.id",
             "start", "end", "gene_name", "transcript_name", "bulk_TPMperc")
result <- unique(reads[gene_name %chin% genes, ..columns])
fwrite(result, args[[3]], sep="\t", col.names=FALSE)
