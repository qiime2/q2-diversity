#!/usr/bin/env Rscript

### LOAD LIBRARIES ###

cat(R.version$version.string, "\n")
args <- commandArgs(TRUE)

suppressWarnings(library(vegan))


### LOAD DATA ###

distances <- read.table(file = args[[1]], sep="\t", header=TRUE, fill=TRUE, row.names=1)
sample.md <- read.table(file = args[[2]], sep="\t", header=TRUE, fill=TRUE, row.names=1, quote="\"", na.strings="")
formula <- args[[3]]
perms <- as.integer(args[[4]])
njobs <- args[[5]]
out.path <- args[[6]]

### RUN ADONIS ###

dm <- as.dist(distances)
formula <- as.formula(paste("dm ~ ", formula))

res <- adonis(formula, data=sample.md, permutations=perms, parallel=njobs)

write.table(res$aov.tab, out.path, sep="\t", append=F, quote=FALSE)

q(status=0)
