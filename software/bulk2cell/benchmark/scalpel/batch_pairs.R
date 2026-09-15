#!/usr/bin/env Rscript
# Extract the observed gene/transcript relation after original filtering.
suppressMessages(suppressWarnings(library(data.table)))
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) stop("usage: batch_pairs.R FILTERED_RDS OUTPUT_TSV")
reads <- as.data.table(readRDS(args[[1]]))
pairs <- unique(reads[, .(gene_name, transcript_name)])
fwrite(pairs, args[[2]], sep="\t", col.names=FALSE)
