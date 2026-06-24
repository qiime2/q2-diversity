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
perms <- as.integer(args[[4]])
njobs <- as.integer(args[[5]])
out.path <- args[[6]]

### RUN ADONIS ###

dm <- as.dist(distances)
formula <- as.formula(paste("dm ~ ", formula))

unit.col <- NULL
target.col <- NULL
use.unit.label.permutations <- FALSE

optional.args <- character(0)

if (length(args) > 6) {
  optional.args <- args[7:length(args)]
}

unit.flag <- "--unit-label-permutations"
unit.flag.idx <- match(unit.flag, optional.args)

if (!is.na(unit.flag.idx)) {
  use.unit.label.permutations <- TRUE

  if (length(optional.args) < unit.flag.idx + 2) {
    stop(paste(
      unit.flag,
      "requires two following arguments:",
      "permutation unit column and permutation target column."
    ))
  }

  unit.col <- optional.args[[unit.flag.idx + 1]]
  target.col <- optional.args[[unit.flag.idx + 2]]
}

get_f_value <- function(aov.tab, term) {
  if (!(term %in% rownames(aov.tab))) {
    stop(paste("Target term not found in adonis2 table:", term))
  }

  if (!("F" %in% colnames(aov.tab))) {
    stop("Expected F column not found in adonis2 result table.")
  }

  return(as.numeric(aov.tab[term, "F"]))
}

if (!use.unit.label.permutations) {

  res <- adonis2(
    formula,
    by='terms',
    data=sample.md,
    permutations=perms,
    parallel=njobs
  )

  write.table(res, out.path, sep="\t", append=F, quote=FALSE)

} else {

  if (!(unit.col %in% colnames(sample.md))) {
    stop(paste("permutation unit column not found in metadata:", unit.col))
  }

  if (!(target.col %in% colnames(sample.md))) {
    stop(paste("permutation target column not found in metadata:", target.col))
  }

  unit.ids <- unique(as.character(sample.md[[unit.col]]))

  unit.labels <- lapply(unit.ids, function(unit.id) {
    rows <- as.character(sample.md[[unit.col]]) == unit.id
    values <- unique(sample.md[[target.col]][rows])

    if (length(values) != 1) {
      stop(paste(
        "permutation target column must be constant within each unit.",
        "Problematic unit:",
        unit.id
      ))
    }

    return(values[[1]])
  })

  unit.labels <- unlist(unit.labels, use.names=FALSE)
  names(unit.labels) <- unit.ids

  observed <- suppressMessages(suppressWarnings(adonis2(
    formula,
    by='terms',
    data=sample.md,
    permutations=1,
    parallel=njobs
  )))

  observed.tab <- observed
  observed.f <- get_f_value(observed.tab, target.col)

  perm.f <- numeric(perms)

  for (i in seq_len(perms)) {
    shuffled.labels <- sample(
      unit.labels,
      size=length(unit.labels),
      replace=FALSE
    )

    names(shuffled.labels) <- names(unit.labels)

    perm.md <- sample.md
    perm.md[[target.col]] <- shuffled.labels[
      match(as.character(sample.md[[unit.col]]), names(shuffled.labels))
    ]

    perm.res <- suppressMessages(suppressWarnings(adonis2(
      formula,
      by='terms',
      data=perm.md,
      permutations=1,
      parallel=njobs
    )))

    perm.f[[i]] <- get_f_value(perm.res, target.col)
  }

  custom.p <- (sum(perm.f >= observed.f) + 1) / (perms + 1)

  if (!("Pr(>F)" %in% colnames(observed.tab))) {
    stop("Expected Pr(>F) column not found in adonis2 result table.")
  }

  observed.tab[target.col, "Pr(>F)"] <- custom.p

  write.table(observed.tab, out.path, sep="\t", append=F, quote=FALSE)
}

q(status=0)
