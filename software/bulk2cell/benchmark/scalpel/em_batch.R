#!/usr/bin/env Rscript
# Execute the original upstream EM function in a shared process, avoiding a new
# interpreter per barcode. Source expressions are parsed, never rewritten.
suppressPackageStartupMessages({library(dplyr); library(data.table); library(tidyr); library(stringr)})
a <- commandArgs(trailingOnly=TRUE)
expressions <- parse(file=file.path(a[1], "src", "em_algorithm.R"))
for (expr in expressions) {
  if (is.call(expr) && is.symbol(expr[[1]]) && as.character(expr[[1]]) %in% c("=", "<-") &&
      is.symbol(expr[[2]]) && as.character(expr[[2]]) %in% c("MAX_IT", "em_algorithm")) eval(expr)
}
stopifnot(exists("em_algorithm"), MAX_IT == 30)
x <- fread(a[2], col.names=c("bc","gene_name","transcript_name","umi","frag_prob_weighted"))
x <- as.data.table(na.omit(unique(x)))
# Preserve upstream per-cell, per-gene grouping and three-decimal rounding.
groups <- split(x, interaction(x$bc, x$gene_name, drop=TRUE))
result <- lapply(groups, function(tab) {
  fit <- em_algorithm(tab)
  if (!nrow(fit)) return(NULL)
  fit$bc <- tab$bc[1]
  fit$rel_abund <- round(fit$rel_abund, 3)
  # Count only molecules with nonzero model probability, unlike Cell Ranger scaling.
  supported <- tab[, .(p=sum(frag_prob_weighted)), by=umi][p>0, .N]
  fit$supported_molecules <- supported
  fit$estimated_count <- fit$rel_abund * supported
  fit
})
fwrite(rbindlist(result), a[3], sep="\t")
# The unmodified upstream function does not expose iteration/convergence metadata.
fwrite(x[, .(candidate_molecules=uniqueN(umi), positive_molecules=uniqueN(umi[frag_prob_weighted>0])), by=.(bc,gene_name)], a[4], sep="\t")
