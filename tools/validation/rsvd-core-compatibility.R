# Focused ad8f05a-to-current comparison. This is a validation experiment, not a
# performance benchmark or a portable promise of bitwise BLAS reproducibility.
# Usage: Rscript rsvd-core-compatibility.R OLD_REPO CURRENT_REPO OUTPUT_PREFIX
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3L)
pkgload::load_all(args[2L], quiet = TRUE, export_all = FALSE)
current <- asNamespace("rsvd")
previous <- new.env(parent = parent.env(current))
for (file in list.files(file.path(args[1L], "R"), pattern = "[.]R$", full.names = TRUE)) {
  sys.source(file, envir = previous)
}

old_functions <- sort(Filter(function(name) is.function(previous[[name]]), ls(previous)))
formals_results <- data.frame(
  function_name = old_functions,
  identical_formals = vapply(old_functions, function(name) {
    exists(name, current, inherits = FALSE) &&
      identical(formals(previous[[name]]), formals(current[[name]]))
  }, logical(1)),
  stringsAsFactors = FALSE
)

exports <- function(repo) {
  directives <- as.list(parse(file.path(repo, "NAMESPACE")))
  selected <- Filter(function(x) identical(x[[1L]], as.name("export")), directives)
  sort(unlist(lapply(selected, function(x) as.character(x[-1L]))))
}
previous_exports <- exports(args[1L])
current_exports <- exports(args[2L])
helpers_identical <- identical(
  readLines(file.path(args[1L], "R", "wrapper_function.R")),
  readLines(file.path(args[2L], "R", "wrapper_function.R")))

basis <- function(m, rank, complex = FALSE) {
  if (complex) return(exp(2i * pi * outer(seq_len(m) - 1, seq_len(rank) - 1) / m) /
                         sqrt(m))
  out <- sqrt(2 / m) * cos(pi * outer(seq_len(m) - 0.5, seq_len(rank) - 1) / m)
  out[, 1L] <- out[, 1L] / sqrt(2)
  out
}
fixture <- function(complex = FALSE) {
  d <- 10^seq(0, -3, length.out = 16)
  A <- basis(36, 16, complex) %*% (d * Conj(t(basis(20, 16, complex))))
  dimnames(A) <- list(paste0("row", seq_len(nrow(A))),
                      paste0("col", seq_len(ncol(A))))
  A
}
real <- fixture()
complex <- fixture(TRUE)
sparse <- matrix(0, 40, 24)
sparse[cbind(seq_len(20), seq_len(20))] <- 10^seq(0, -3, length.out = 20)
dimnames(sparse) <- list(paste0("row", seq_len(nrow(sparse))),
                         paste0("col", seq_len(ncol(sparse))))
sparse <- Matrix::Matrix(sparse, sparse = TRUE)
inputs <- list(real_tall = real, real_wide = t(real),
               complex_tall = complex, complex_wide = Conj(t(complex)),
               sparse_tall = sparse, sparse_wide = Matrix::t(sparse))
vectors <- list(default = list(nu = NULL, nv = NULL),
                omit_left = list(nu = 0, nv = 3),
                omit_right = list(nu = 2, nv = 0),
                omit_both = list(nu = 0, nv = 0),
                truncated = list(nu = 2, nv = 3))
rows <- list()
for (input_name in names(inputs)) {
  for (sdist in c("normal", "unif", "rademacher")) {
    for (q in c(0L, 2L)) {
      for (vector_case in names(vectors)) {
        call_args <- c(list(A = inputs[[input_name]], k = 4L, p = 3L,
                            q = q, sdist = sdist), vectors[[vector_case]])
        seed <- 9100L + length(rows)
        set.seed(seed)
        old <- do.call(previous$rsvd.default, call_args)
        old_rng <- .Random.seed
        set.seed(seed)
        new <- do.call(current$rsvd.default, call_args)
        new_rng <- .Random.seed
        rows[[length(rows) + 1L]] <- data.frame(
          input = input_name, sdist = sdist, q = q, vectors = vector_case,
          exact_output = identical(old, new), exact_rng = identical(old_rng, new_rng),
          stringsAsFactors = FALSE)
      }
    }
  }
}
results <- do.call(rbind, rows)
write.csv(results, paste0(args[3L], "-outputs.csv"), row.names = FALSE)
write.csv(formals_results, paste0(args[3L], "-formals.csv"), row.names = FALSE)
summary <- c(
  R.version.string,
  paste("Matrix", as.character(packageVersion("Matrix"))),
  paste("Comparison cases:", nrow(results)),
  paste("Exactly identical output objects:", sum(results$exact_output)),
  paste("Exactly identical final RNG states:", sum(results$exact_rng)),
  paste("Existing function formals:", sum(formals_results$identical_formals),
        "of", nrow(formals_results), "identical"),
  paste("Exports:", length(previous_exports), "before;", length(current_exports),
        "after; identical:", identical(previous_exports, current_exports)),
  paste("Unchanged helper source:", helpers_identical)
)
writeLines(summary, paste0(args[3L], "-summary.txt"))
cat(paste(summary, collapse = "\n"), "\n")
stopifnot(all(results$exact_output), all(results$exact_rng),
          all(formals_results$identical_formals),
          identical(previous_exports, current_exports), helpers_identical)
