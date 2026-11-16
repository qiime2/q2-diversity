#!/usr/bin/env Rscript

### ERROR HANDLING ###
options(error = function() {
  sink(stderr())
  on.exit(sink(NULL))
  traceback(3)
  if (!interactive()) {
    q(status = 1)
  }
})

### LOAD LIBRARIES ###

cat(R.version$version.string, "\n")
args <- commandArgs(TRUE)

suppressWarnings(library(vegan))


### LOAD DATA ###

distances <- read.table(file = args[[1]], sep="\t", header=TRUE, fill=TRUE, row.names=1)
sample.md <- read.table(file = args[[2]], sep="\t", header=TRUE, fill=TRUE, row.names=1, quote="\"", na.strings="")
formula <- args[[3]]
by <- args[[4]]
perms <- as.integer(args[[5]])
njobs <- as.integer(args[[6]])
out.path <- args[[7]]

### RUN ADONIS ###

dm <- as.dist(distances)
formula <- as.formula(paste("dm ~ ", formula))
if (by == "null") {
  by <- NULL
}

res <- adonis2(formula, by=by, data=sample.md, permutations=perms, parallel=njobs)

write.table(res, out.path, sep="\t", append=F, quote=FALSE)

q(status=0)
