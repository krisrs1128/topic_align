#' Simulate Standard K-Topics Model
k_topics <- function(K, V = 500, lambda = 1) {
  rdirichlet(K, rep(lambda, V)) %>%
    flatten_beta() %>%
    widen_betas()
}

widen_betas <- function(betas) {
  betas %>%
    pivot_wider(m:k_LDA, names_from = "w", values_from = "b") %>%
    mutate(k = as.numeric(as.factor(k_LDA))) %>%
    select(-K, -k_LDA)
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

#' Simulation functions for functional equivalence experiment

perturb_topic <- function(beta_k, species_subsets, nu=2) {
  beta_k[species_subsets[[1]]] <-  nu * beta_k[species_subsets[[1]]]
  beta_k[species_subsets[[2]]] <-  (1 / nu) * beta_k[species_subsets[[2]]]
  beta_k / sum(beta_k)
}

perturb_topics <- function(B, n_per_topic = NULL, subset_size = 30, nus=c(1, 4)) {
  K <- nrow(B)
  stopifnot(subset_size < ncol(B))
  if (is.null(n_per_topic)) {
    n_per_topic <- c(rep(20, 2), rep(1, K - 2))
  }
  
  B_tilde <- vector(length = K, mode = "list")
  for (k in seq_len(K)) {
    s <- list(1:(subset_size / 2), (subset_size / 2 + 1):subset_size)
    class_k <- matrix(nrow = n_per_topic[k], ncol = ncol(B))
    nu <- runif(n_per_topic[k], nus[1], nus[2])
    for (i in seq_len(n_per_topic[k])) {
      class_k[i, ] <- perturb_topic(B[k, ], s, nu[i])
    }
    B_tilde[[k]] <- as_tibble(class_k)
  }
  
  bind_rows(B_tilde, .id = "k") %>%
    pivot_longer(-k, names_to = "w", values_to = "b")
}

equivalence_data <- function(N, V, K, lambdas, n0 = NULL, ...) {
  if (is.null(n0)) {
    n0 <- rpois(N, 1000)
  }
  
  B <- rdirichlet(K, rep(lambdas$beta, V))
  B_tilde <- perturb_topics(B, ...)
  gammas <- rdirichlet(N, rep(lambdas$gamma, K))
  x <- matrix(nrow = N, ncol = V)
  
  for (i in seq_len(N)) {
    B_i <- B_tilde %>%
      group_by(k, w) %>%
      sample_n(1) %>%
      pivot_wider(k, names_from = "w", values_from = "b") %>%
      ungroup() %>%
      select(-k) %>%
      as.matrix()
    x[i, ] <- rmultinom(1, n0[i], t(B_i) %*% gammas[i, ])
  }
  
  colnames(x) <- seq_len(ncol(x))
  rownames(x) <- seq_len(nrow(x))
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
