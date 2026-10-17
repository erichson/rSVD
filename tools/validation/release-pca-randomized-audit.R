# Independent feasible primal/dual bounds for randomized robust PCA.
# Usage: Rscript tools/validation/release-pca-randomized-audit.R LIBRARY OUTPUT_CSV
args <- commandArgs(TRUE)
library(rsvd, lib.loc = args[1L])
solver <- getFromNamespace('rrpca.default', 'rsvd')
b <- body(solver)
b[[length(b)]] <- quote(list(
  L = L, S = S, err = rrpcaObj$err, Z = Z, dualErr = dualErr,
  iterations = niter - 1L, use_exact = use_exact, refinement = refinement
))
body(solver) <- b
results <- list()
for (j in 1:6) {
  set.seed(1249 + j)
  A <- matrix(rnorm(60 * 40), 60, 40)
  lambda <- if (j <= 3) 1 / sqrt(60) else c(.03, .3, 1.3)[j - 3]
  fit <- solver(A, lambda = lambda, tol = 1e-7, maxiter = 1000)
  violation <- max(svd(fit$Z, nu = 0, nv = 0)$d[1],
                   max(abs(fit$Z)) / lambda)
  lower <- sum(A * fit$Z) / max(1, violation)
  upper <- sum(svd(fit$L, nu = 0, nv = 0)$d) +
    lambda * sum(abs(A - fit$L))
  results[[j]] <- data.frame(
    seed = 1249 + j, lambda = lambda, iterations = fit$iterations,
    primal = tail(fit$err, 1), stationarity = fit$dualErr,
    normalized_gap = (upper - lower) / upper,
    dual_infeasibility = max(0, violation - 1),
    refinement = fit$refinement,
    stagnation_exact_fallback = fit$use_exact
  )
}
print(do.call(rbind, results), digits = 8)
write.csv(do.call(rbind, results), args[2L], row.names = FALSE)
