#' Simulate Standard K-Topics Model
k_topics <- function(K, V = 500, lambda = 1) {
  rdirichlet(K, rep(lambda, V)) %>%
    flatten_beta() %>%
    widen_betas()
}

simulate_lda <- function(betas, gammas, n0=NULL) {
  n <- nrow(gammas)
  if (is.null(n0)) {
    n0 <- rpois(n, 1000)
  }

  x <- matrix(nrow = n, ncol = ncol(betas))
  for (i in seq_len(n)) {
    x[i, ] <- rmultinom(1, n0[i], t(betas) %*% gammas[i, ])
  }
  rownames(x) <- seq_len(n)
  colnames(x) <- seq_len(ncol(betas))
  x
}

.simulate_replicate <- function(betas, gammas, M = 5, mechanism = NULL) {
  if (is.null(mechanism)) {
    mechanism <- simulate_lda
  }
  x <- mechanism(betas, gammas)
  wrapper <- fit_wrapper(x, M)
  endpoint_topics(wrapper$fits$betas, wrapper$alignment)
}

simulate_replicates <- function(betas, gammas, n_reps) {
  map_dfr(seq_len(n_reps), ~ .simulate_replicate(betas, gammas), .id = "replicate")
}

fit_wrapper <- function(x, M) {
  fits <- run_lda_models(x, 1:M, c(.1), "VEM", 123, reset = TRUE)
  list(fits = fits, alignment = align_topics(fits))
}

vis_wrapper <- function(x, M) {
  fits <- run_lda_models(x, 1:M, c(.1), "VEM", 123, reset = TRUE)
  alignment <- align_topics(fits)
  list(
    fits = fits,
    alignment = alignment,
    p1 = visualize_aligned_topics(alignment),
    p2 = visualize_aligned_topics(alignment, method = "beta_alignment")
  )
}
