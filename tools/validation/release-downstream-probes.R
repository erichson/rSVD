# Focused execution of the interfaces that the reviewed downstream packages use.
# Some large workflows need unrelated executables, optimization, or training.
# In those cases this runs an unchanged expression extracted from the installed
# function body; the output explicitly distinguishes those from whole functions.
# Usage: Rscript release-downstream-probes.R WORK_DIR RSVD_LIBRARY LABEL
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3L)
work <- normalizePath(args[[1L]])
rsvd_library <- normalizePath(args[[2L]])
label <- args[[3L]]
.libPaths(c(rsvd_library, file.path(work, "dependency-library"), .libPaths()))
stopifnot(find.package("rsvd") == file.path(rsvd_library, "rsvd"))
output_dir <- file.path(work, paste0("probes-", label))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
options(quanteda_threads = 1L)
records <- list()
payloads <- list()
probe <- function(package, name, expression, scope = "whole function") {
  expression <- substitute(expression)
  messages <- character()
  result <- tryCatch(withCallingHandlers({
    stopifnot(requireNamespace(package, quietly = TRUE))
    set.seed(217031L)
    eval(expression, envir = parent.frame())
  }, warning = function(w) {
    messages <<- c(messages, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = identity)
  ok <- !inherits(result, "error")
  records[[name]] <<- data.frame(
    package = package, version = as.character(packageVersion(package)),
    rsvd_version = as.character(packageVersion("rsvd")),
    probe = name, scope = scope, passed = ok,
    warnings = paste(messages, collapse = " | "),
    error = if (ok) "" else conditionMessage(result), stringsAsFactors = FALSE)
  if (ok) payloads[[name]] <<- result
  cat(package, name, if (ok) "PASS" else conditionMessage(result), "\n")
  write.csv(do.call(rbind, records), file.path(output_dir, "summary.csv"), row.names = FALSE)
  saveRDS(payloads, file.path(output_dir, "outputs.rds"))
}
finite <- function(x) stopifnot(all(is.finite(x)))
relative_error <- function(actual, expected)
  norm(as.matrix(actual) - as.matrix(expected), "F") /
    max(1, norm(as.matrix(expected), "F"))
check_svd <- function(fit, A, k, tolerance = 1e-10) {
  stopifnot(length(fit$d) == k,
            identical(dim(fit$u), c(nrow(A), as.integer(k))),
            identical(dim(fit$v), c(ncol(A), as.integer(k))))
  finite(c(fit$d, fit$u, fit$v))
  stopifnot(norm(crossprod(fit$u) - diag(k), "F") < tolerance,
            norm(crossprod(fit$v) - diag(k), "F") < tolerance)
  reconstruction <- sweep(fit$u, 2L, fit$d, "*") %*% t(fit$v)
  stopifnot(relative_error(reconstruction, A) < tolerance)
  list(d = fit$d, reconstruction = unname(as.matrix(reconstruction)),
       left_projector = unname(tcrossprod(fit$u)),
       right_projector = unname(tcrossprod(fit$v)))
}
extract_statement <- function(package, function_name, marker, occurrence = 1L) {
  statements <- as.list(body(getFromNamespace(function_name, package)))[-1L]
  selected <- vapply(statements, function(x)
    grepl(marker, paste(deparse(x), collapse = "\n"), fixed = TRUE), logical(1))
  hits <- statements[selected]
  stopifnot(length(hits) >= occurrence)
  hits[[occurrence]]
}

set.seed(10217)
U <- qr.Q(qr(matrix(rnorm(60L * 4L), 60L, 4L)))
V <- qr.Q(qr(matrix(rnorm(40L * 4L), 40L, 4L)))
A <- sweep(U, 2L, c(20, 10, 5, 2), "*") %*% t(V)
dimnames(A) <- list(paste0("row", seq_len(nrow(A))), paste0("col", seq_len(ncol(A))))
PSD <- tcrossprod(sweep(U, 2L, sqrt(c(20, 10, 5, 2)), "*"))
dimnames(PSD) <- list(rownames(A), rownames(A))
counts <- matrix(0, 60L, 40L)
for (j in seq_len(4L))
  counts[(15L * (j - 1L) + 1L):(15L * j),
         (10L * (j - 1L) + 1L):(10L * j)] <- j
dimnames(counts) <- dimnames(A)

probe("FPCdpca", "distributed_svd", {
  X <- A[, seq_len(12L)] + matrix(rnorm(60L * 12L, sd = 0.02), 60L)
  fit <- FPCdpca::Drsvd(X, K = 3, nk = 20, m = 4, q = 2, k = 4)
  stopifnot(identical(names(fit), c("MSEXrsvd", "MSEvrsvd", "MSESrsvd", "kopt")),
            fit[["kopt"]] %in% seq_len(3L))
  finite(fit)
  fit
})

probe("harf", "two_pca_projection_branches", {
  env <- new.env(parent = asNamespace("harf"))
  env$meta_features <- as.data.frame(A[, seq_len(8L)])
  env$num_btwn_pcs <- 2L
  env$target <- NULL
  eval(extract_statement("harf", "h_arf", "pca_meta <- rpca"), env)
  stopifnot(identical(dim(env$meta_features), c(60L, 2L)))
  env$omx_onset_data <- as.data.frame(A[, 9:14])
  env$num_onset_pcs <- 2L
  eval(extract_statement("harf", "h_arf", "pca_onset <- rpca"), env)
  stopifnot(identical(dim(env$meta_features), c(60L, 4L)),
            identical(dim(env$onset_loadings), c(6L, 2L)))
  finite(as.matrix(env$meta_features))
  list(projected = as.matrix(env$meta_features), onset = env$onset_loadings)
}, scope = "two extracted PCA branches; no forest training")

probe("jackstraw", "randomized_pca_resampling", {
  X <- t(A[, seq_len(25L)]) + matrix(rnorm(25L * 60L, sd = 0.02), 25L)
  fit <- jackstraw::jackstraw_rpca(X, r = 3, s = 3, B = 4, verbose = FALSE)
  stopifnot(length(fit$p.value) == 25L, identical(dim(fit$null.stat), c(3L, 4L)),
            all(fit$p.value >= 0 & fit$p.value <= 1))
  finite(c(fit$p.value, fit$obs.stat, fit$null.stat))
  fit[c("p.value", "obs.stat", "null.stat")]
})

probe("MUGS", "two_site_svd_initializers", {
  env <- new.env(parent = asNamespace("MUGS"))
  env$S.1 <- PSD
  env$S.2 <- 2 * PSD
  env$p <- 2L
  for (marker in c("svd.1 <- rsvd", "svd.2 <- rsvd",
                   "emb.1 <- get_embed", "emb.2 <- get_embed"))
    eval(extract_statement("MUGS", "MUGS", marker), env)
  stopifnot(identical(dim(env$emb.1), c(60L, 2L)),
            identical(dim(env$emb.2), c(60L, 2L)),
            max(abs(rowSums(env$emb.1^2) - 1)) < 1e-12)
  finite(c(env$emb.1, env$emb.2))
  list(site1_gram = tcrossprod(env$emb.1), site2_gram = tcrossprod(env$emb.2))
}, scope = "four extracted initializer statements; no MUGS optimizer")

probe("nlpembeds", "sparse_psd_embedding", {
  fit <- nlpembeds::get_svd(Matrix::Matrix(PSD, sparse = TRUE),
                            embedding_dim = 4, svd_rank = 8)
  stopifnot(identical(dim(fit), c(60L, 4L)),
            identical(rownames(fit), rownames(PSD)))
  finite(fit)
  stopifnot(relative_error(tcrossprod(fit), PSD) < 1e-10)
  unname(tcrossprod(fit))
})

probe("sparsepca", "randomized_spca_unpenalized", {
  fit <- sparsepca::rspca(A, k = 4, alpha = 0, beta = 0, center = TRUE,
                          scale = FALSE, o = 4, max_iter = 1000, verbose = FALSE)
  ref <- prcomp(A, center = TRUE, scale. = FALSE, rank. = 4)
  stopifnot(identical(dim(fit$loadings), c(40L, 4L)),
            max(abs(fit$sdev - ref$sdev[seq_len(4L)])) < 1e-8)
  reconstruction <- fit$scores %*% t(fit$transform)
  finite(reconstruction)
  stopifnot(relative_error(reconstruction, scale(A, scale = FALSE)) < 1e-8)
  list(sdev = fit$sdev, reconstruction = unname(reconstruction))
})

probe("SPECK", "manual_reduced_rank_reconstruction", {
  fit <- SPECK::randomizedRRR(counts, rank.range.end = 8,
                              manual.rank = 4, seed.rsvd = 13)
  stopifnot(fit$rrr.rank == 4, identical(dim(fit$rrr.mat), dim(counts)),
            identical(rownames(fit$rrr.mat), rownames(counts)),
            identical(colnames(fit$rrr.mat), colnames(counts)))
  expected <- log1p(counts / rowSums(counts) * 10000)
  finite(fit$rrr.mat)
  stopifnot(relative_error(fit$rrr.mat, expected) < 1e-10)
  unname(fit$rrr.mat)
})

probe("TCA", "refactor_randomized_pca", {
  X <- t(A) + matrix(rnorm(length(A), sd = 0.02), 40L)
  fit <- TCA::refactor(X, k = 3, sparsity = 15, num_comp = 3,
                       sd_threshold = NULL, rand_svd = TRUE,
                       log_file = file.path(output_dir, "TCA.log"), verbose = FALSE)
  stopifnot(identical(dim(fit$scores), c(60L, 3L)),
            identical(dim(fit$coeffs), c(15L, 3L)),
            setequal(fit$ranked_list, rownames(X)))
  finite(c(fit$scores, fit$coeffs))
  list(scores_gram = tcrossprod(fit$scores),
       coeffs_gram = tcrossprod(fit$coeffs), ranked_list = fit$ranked_list)
})

probe("text2map", "full_width_procrustes_svd", {
  X <- matrix(rnorm(60L * 6L), 60L, 6L)
  rotation <- qr.Q(qr(matrix(rnorm(36L), 6L)))
  reference <- X %*% rotation
  fit <- text2map::find_transformation(X, reference, method = "align")
  expected <- scale(reference / sqrt(rowSums(reference^2)), scale = FALSE)
  finite(fit)
  stopifnot(relative_error(fit, expected) < 1e-10)
  unname(as.matrix(fit))
})

probe("text2map", "truncated_retrofit_embedding", {
  embedding <- U
  rownames(embedding) <- rownames(PSD)
  fit <- text2map::find_transformation(embedding,
    Matrix::Matrix(PSD, sparse = TRUE), method = "retrofit")
  stopifnot(identical(dim(fit), c(60L, 4L)),
            identical(rownames(fit), rownames(PSD)),
            max(abs(rowSums(fit^2) - 1)) < 1e-12)
  finite(fit)
  unname(as.matrix(fit))
})

for (wide in c(FALSE, TRUE)) {
  X <- if (wide) t(counts) else counts
  probe("wordvector", paste0("lsa_right_vectors_", if (wide) "wide" else "tall"), {
    x <- quanteda::as.dfm(Matrix::Matrix(X, sparse = TRUE))
    fit <- wordvector::textmodel_lsa(x, dim = 4L, min_count = 0L,
                                     engine = "rsvd", weight = "count")
    vectors <- fit$values$word
    stopifnot(identical(dim(vectors), c(ncol(X), 4L)), fit$engine == "rsvd")
    finite(vectors)
    stopifnot(relative_error(X %*% tcrossprod(vectors), X) < 1e-10)
    unname(tcrossprod(vectors))
  })
  probe("LSX", paste0("lss_uncached_right_vectors_", if (wide) "wide" else "tall"), {
    x <- quanteda::as.dfm(Matrix::Matrix(X, sparse = TRUE))
    fit <- getFromNamespace("cache_svd", "LSX")(
      x, k = 4L, weight = "count", engine = "rsvd", cache = FALSE)
    stopifnot(is.null(fit$u), identical(dim(fit$v), c(ncol(X), 4L)), length(fit$d) == 4L)
    finite(c(fit$d, fit$v))
    stopifnot(relative_error(X %*% tcrossprod(fit$v), X) < 1e-10)
    list(d = fit$d, right_projector = unname(tcrossprod(fit$v)))
  }, scope = "whole internal SVD adapter, cache disabled")
}

probe("multivarious", "svd_wrapper", {
  fit <- multivarious::svd_wrapper(A, ncomp = 4, method = "rsvd")
  stopifnot(inherits(fit, "bi_projector"), identical(dim(fit$v), c(40L, 4L)),
            identical(dim(fit$s), c(60L, 4L)), length(fit$sdev) == 4L)
  finite(c(fit$v, fit$s, fit$sdev))
  stopifnot(relative_error(fit$s %*% t(fit$v), A) < 1e-10)
  list(sdev = fit$sdev, reconstruction = unname(fit$s %*% t(fit$v)))
})

probe("flashier", "low_rank_input_conversion", {
  fit <- flashier::flash_init(rsvd::rsvd(A, k = 4))
  internal <- flashier::flash_fit(fit)
  stopifnot(inherits(internal$Y, "lowrank"), length(internal$Y) == 2L)
  reconstructed <- internal$Y[[1L]] %*% t(internal$Y[[2L]])
  finite(reconstructed)
  stopifnot(relative_error(reconstructed, A) < 1e-10)
  unname(reconstructed)
})

probe("eigencore", "validation_svd_adapter", {
  fit <- getFromNamespace("run_svd_method", "eigencore")(
    "rsvd", A, rank = 4, tol = 1e-8)
  check_svd(fit, A, 4L)
}, scope = "whole internal validation adapter; no timing")

probe("Seurat", "fitsne_pca_initializer", {
  env <- new.env(parent = asNamespace("Seurat"))
  env$X <- A
  env$dims <- 2L
  env$rand_seed <- 29L
  env$initialization <- "pca"
  eval(extract_statement("Seurat", "fftRtsne", "rsvd_out <- rsvd::rsvd"), env)
  stopifnot(identical(dim(env$initialization), c(60L, 2L)),
            abs(sd(env$initialization[, 1L]) - 0.0001) < 1e-12)
  finite(env$initialization)
  unname(tcrossprod(env$initialization))
}, scope = "extracted PCA initializer; FIt-SNE executable not run")

probe("stm", "spectral_pca_initializer", {
  env <- new.env(parent = asNamespace("stm"))
  env$Qbar <- A
  env$init.dims <- 4L
  eval(extract_statement("stm", "tsneAnchor", "Xpca <- rsvd::rpca"), env)
  stopifnot(identical(dim(env$Xpca), c(60L, 4L)))
  finite(env$Xpca)
  ref <- prcomp(A, rank. = 4L)$x
  stopifnot(relative_error(tcrossprod(env$Xpca), tcrossprod(ref)) < 1e-10)
  unname(tcrossprod(env$Xpca))
}, scope = "extracted first PCA initializer; t-SNE not run")

if (requireNamespace("BiocSingular", quietly = TRUE)) {
  representations <- list(
    dense = A,
    sparse = Matrix::Matrix(counts, sparse = TRUE),
    delayed_dense = DelayedArray::DelayedArray(A),
    delayed_sparse = DelayedArray::DelayedArray(Matrix::Matrix(counts, sparse = TRUE)))
  for (representation in names(representations)) {
    X <- representations[[representation]]
    for (deferred in c(FALSE, TRUE)) {
      probe("BiocSingular", paste("random_svd", representation, deferred, sep = "_"), {
        fit <- BiocSingular::runRandomSVD(X, k = 4, center = TRUE,
          scale = TRUE, deferred = deferred, fold = Inf)
        expected <- scale(as.matrix(X), center = TRUE, scale = TRUE)
        check_svd(fit, expected, 4L)
      })
    }
  }
  probe("BiocSingular", "omitted_vectors_and_extra_values", {
    fit <- BiocSingular::runRandomSVD(DelayedArray::DelayedArray(A),
      k = 4L, nu = 0L, nv = 2L, fold = Inf)
    stopifnot(length(fit$d) == 4L, identical(dim(fit$u), c(60L, 0L)),
              identical(dim(fit$v), c(40L, 2L)))
    finite(c(fit$d, fit$v))
    list(d = fit$d, right_projector = unname(tcrossprod(fit$v)))
  })
  probe("BiocSingular", "wide_cross_product_adapter", {
    fit <- BiocSingular::runRandomSVD(t(A), k = 4L, fold = 1)
    check_svd(fit, t(A), 4L)
  })
  probe("BiocSingular", "random_param_pca", {
    fit <- BiocSingular::runPCA(DelayedArray::DelayedArray(A), rank = 4L,
      BSPARAM = BiocSingular::RandomParam(deferred = TRUE), center = TRUE,
      scale = FALSE)
    ref <- prcomp(A, rank. = 4L)
    stopifnot(identical(dim(fit$x), c(60L, 4L)),
              identical(dim(fit$rotation), c(40L, 4L)),
              max(abs(fit$sdev - ref$sdev[seq_len(4L)])) < 1e-10,
              relative_error(fit$x %*% t(fit$rotation), scale(A, scale = FALSE)) < 1e-10)
    list(sdev = fit$sdev, scores_gram = unname(tcrossprod(fit$x)))
  })
}
writeLines(capture.output(sessionInfo()), file.path(output_dir, "session-info.txt"))
cat("Passed", sum(vapply(records, function(x) x$passed, logical(1))), "of",
    length(records), "focused downstream probes\n")
if (!all(vapply(records, function(x) x$passed, logical(1)))) quit(status = 1L)
