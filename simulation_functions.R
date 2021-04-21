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
perturb_topics <- function(B, n_per_topic = NULL, subset_size = 10, nu_max = 0.3, ...) {
  K <- nrow(B)
  if (is.null(n_per_topic)) {
    n_per_topic <- c(rep(20, 2), rep(1, K - 2))
  }

  B_tilde <- list()
  for (k in seq_len(K)) {
    B_tilde[[k]] <- rep(1, n_per_topic[k]) %*% matrix(B[k, ], nrow = 1)
    if (n_per_topic[k] == 1) next

    nu <- runif(n_per_topic[k], 0, nu_max)
    for (i in seq_len(n_per_topic[k])) {
      ix <- c(1:(subset_size / 2), (subset_size / 2 + 1):subset_size)
      B_tilde[[k]][i, ix[1]] <- (1 + nu[i]) * B_tilde[[k]][i, ix[1]]
      B_tilde[[k]][i, ix[2]] <- (1 - nu[i]) * B_tilde[[k]][i, ix[2]]
      B_tilde[[k]][i, ] <- B_tilde[[k]][i, ] / sum(B_tilde[[k]][i, ])
    }
  }

  B_tilde
}

sample_topics <- function(B_tilde) {
  B <- matrix(nrow = length(B_tilde), ncol = ncol(B_tilde[[1]]))
  for (k in seq_along(B_tilde)) {
    ix <- sample(seq_len(nrow(B_tilde[[k]])), 1)
    B[k, ] <- B_tilde[[k]][ix, ]
  }
  B
}

simulate_document <- function(B, gamma_i, n0) {
  z <- rmultinom(1, n0, gamma_i)
  x <- vector(length = ncol(B))
  for (k in seq_along(z)) {
    x <- x + rmultinom(1, z[k], B[k, ])
  }
  x
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
    B_i <- sample_topics(B_tilde)
    x[i, ] <- simulate_document(B_i, gammas[i, ], n0[i])
  }

  colnames(x) <- seq_len(ncol(x))
  rownames(x) <- seq_len(nrow(x))
  list(x = x, B = B, B_tilde = B_tilde, gammas = gammas)
}

#' Multi-environment simulation functions
hierarchical_dirichlet <- function(N_e, K, lambdas_pool = 1, lambda_e = NULL) {
  if (is.null(lambda_e)) {
    lambda_e <- rep(0.1, length(N_e))
  }

  gammas <- list()
  for (e in seq_along(N_e)) {
    Gamma_e <- rdirichlet(1, rep(lambdas_pool, K))
    gammas[[e]] <- rdirichlet(N_e[e], lambda_e[e] * Gamma_e)
  }

  gammas
}

environment_shifts <- function(N_e, K, V, lambdas, ...) {
  B <- rdirichlet(K, rep(lambdas$beta, V))
  gammas <- hierarchical_dirichlet(N_e, K, lambdas$pool, lambdas$e)
  x <- list()
  for (e in seq_along(gammas)) {
    x[[e]] <- simulate_lda(B, gammas[[e]], ...) %>%
      as_tibble()
  }

  x <- bind_rows(x, .id = "environment")
  gammas <- map_dfr(gammas, ~ as_tibble(.), .id = "environment")
  list(x = x, B = B, gammas = gammas)
}

#' General simulation functions
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

merge_betas <- function(betas, beta_hats, K_filter) {
  beta_hats <- beta_hats %>%
    filter(K == K_filter) %>%
    pivot_wider(k_LDA, names_from = w, values_from = b) %>%
    mutate(estimate = TRUE)

  bind_rows(betas, beta_hats) %>%
    select(k_LDA, estimate, everything())
}

betas_umap <- function(beta_compare, ...) {
  rec <- recipe(~ ., data = beta_compare) %>%
    update_role(k_LDA, estimate, i, new_role = "id") %>%
    step_umap(all_predictors(), ...)
  umap_res <- prep(rec)
  list(umap = umap_res, scores = juice(umap_res))
}

align_pairs <- function(betas, masses, reg=1e-3) {
  # keep track of all pairs
  N <- length(betas)
  pairs <- data.frame(
    group1 = rep(seq_len(N), each = N),
    group2 = rep(seq_len(N), N)
  ) %>%
    filter(group2 > group1) %>%
    mutate(pair = row_number())

  # perform beta-alignment
  alignments <- map2_dfr(
    pairs$group1, pairs$group2,
    ~ .beta_weights(betas[c(.x, .y)], masses[c(.x, .y)], reg),
    .id = "pair"
  ) %>%
    mutate(pair = as.integer(pair)) %>%
    left_join(pairs) %>%
    unite(source, group1, k, remove = F) %>%
    unite(target, group2, k_next, remove = F)
}

graph_data <- function(alignments, masses) {
  masses_df <- masses %>%
    map_dfr(~ tibble(mass = .) %>% mutate(k = as.factor(row_number())), .id = "group")
  nodes <- bind_rows(
    alignments %>%
      select(source, group1, k) %>%
      rename(name = source, group = group1),
    alignments %>%
      select(target, group2, k_next) %>%
      rename(name = target, group = group2, k = k_next)
  ) %>%
    unique() %>%
    mutate(group = as.factor(group), k = as.factor(k)) %>%
    left_join(masses_df)

  tbl_graph(
    edges = alignments %>% select(source, target, weight),
    nodes = nodes
  )
}
