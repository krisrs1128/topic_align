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

simulate_lda <- function(betas, gammas, n0=NULL, lambda=1e4) {
  n <- nrow(gammas)
  if (is.null(n0)) {
    n0 <- rpois(n, lambda)
  }

  x <- matrix(nrow = n, ncol = ncol(betas))
  for (i in seq_len(n)) {
    x[i, ] <- rmultinom(1, n0[i], t(betas) %*% gammas[i, ])
  }
  rownames(x) <- seq_len(n)
  colnames(x) <- seq_len(ncol(betas))
  x
}

simulate_gradient <- function(N, K, V, lambdas, alpha=1, n0=NULL, lambda=NULL) {
  B <- rdirichlet(K, rep(lambdas$beta, V))
  gammas <- rdirichlet(N, rep(lambdas$gamma, K))
  n <- nrow(gammas)
  if (is.null(n0)) {
    n0 <- rpois(n, lambda)
  }

  x <- matrix(nrow = n, ncol = ncol(B))
  for (i in seq_len(n)) {
    nu_i <- matrix(rdirichlet(1, rep(lambdas$nu, ncol(B))), ncol = 1)
    mu <- alpha * t(B) %*% gammas[i, ] + (1 - alpha) * nu_i
    x[i, ] <- rmultinom(1, n0[i], mu)
  }

  rownames(x) <- seq_len(n)
  colnames(x) <- seq_len(ncol(B))
  list(B = B, gammas = gammas, x = x)
}

cosine_similarity <- function(X, Y) {
  sim <- matrix(nrow = nrow(X), ncol = nrow(Y))
  for (i in seq_len(nrow(X))) {
    for (j in seq_len(nrow(Y))) {
      sim[i, j] <- sum(X[i, ] * Y[j, ]) / sqrt(sum(X[i, ] ^ 2) * sum(Y[j, ] ^ 2))
    }
  }
  sim
}

#' Simulation functions for functional equivalence experiment
perturb_topics <- function(B, n_per_topic = NULL, subset_size = 40, alpha = 0.1) {
  K <- nrow(B)
  if (is.null(n_per_topic)) {
    n_per_topic <- c(rep(2, 2), rep(1, K - 2))
  }

  B_tilde <- list()
  for (k in seq_len(K)) {
    B_tilde[[k]] <- rep(1, n_per_topic[k]) %*% matrix(B[k, ], nrow = 1)
    if (n_per_topic[k] == 1) next

    for (i in seq_len(n_per_topic[k])) {
      ix <- c(1:(subset_size / 2), (subset_size / 2 + 1):subset_size)
      B_tilde[[k]][i, ] <- perturb_topic(B_tilde[[k]][i, ], ix, alpha)
    }
  }

  B_tilde
}

perturb_topic <- function(Bk, species_ix, alpha = 0.1) {
  alpha <- rep(alpha, length(species_ix))
  sub_community <- rdirichlet(1, alpha)
  Bk[species_ix] <- sum(Bk[species_ix]) * sub_community
  Bk
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
  fits <- run_lda_models(x, 1:M, c(.1), "VEM", 123, reset = FALSE)
  alignment <- align_topics(fits)
  list(
    fits = fits,
    alignment = alignment,
    p1 = visualize_aligned_topics(alignment),
    p2 = visualize_aligned_topics(alignment, method = "beta_alignment")
  )
}

align_pairs <- function(betas, masses, reg=1e-3) {
  N <- length(betas)
  pairs <- data.frame(
    group1 = rep(seq_len(N), each = N),
    group2 = rep(seq_len(N), N)
  ) %>%
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
    unite(target, group2, k_next, remove = F) %>%
    filter(group1 != group2) %>%
    group_by(source) %>%
    mutate(norm_weight = weight / sum(weight, na.rm = TRUE)) %>%
    ungroup()
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
    edges = alignments %>% select(source, target, norm_weight, weight),
    nodes = nodes
  )
}

multimodal_gammas <- function(n, ks, alphas, k_shared=4, alpha_shared = 1) {
  stopifnot(min(ks) > 2)

  gamma_shared <- matrix(0, n, k_shared)
  for (i in seq_len(n)) {
    v <- rbeta(k_shared, 1, alpha_shared)
    gamma_shared[i, 1] <- v[1]
    for (k in 2:k_shared) {
      gamma_shared[i, k] <- v[k] * prod(1 - v[1:(k - 1)])
    }
  }

  used_up <- rowSums(gamma_shared)
  gammas <- map(ks, ~ matrix(0, n, . - 1))
  for (j in seq_along(gammas)) {
    for (i in seq_len(n)) {
      v <- rbeta(ks[j], 1, alphas[j])
      gammas[[j]][i, 1] <- (1 - used_up[i]) * v[1]
      for (k in 2:(ks[j] - 1)) {
        gammas[[j]][i, k] <- v[k] * (1 - used_up[i]) * prod(1 - v[1:(k - 1)])
      }
    }
  }

  map(gammas, ~ cbind(gamma_shared, .)) %>%
    map(~ cbind(., 1 - rowSums(.))) %>%
    map(~ set_colnames(., letters[1:ncol(.)]))
}

multimodal_scatter_data <- function(gamma_hat, weights) {
  scatter_data <- list()
  for (i in seq_len(nrow(weights))) {
    scatter_data[[i]] <- bind_cols(
      gamma_source = gamma_hat %>%
        filter(m == 1) %>%
        pull(str_c("V", weights$source[i])),
      gamma_target = gamma_hat %>%
        filter(m == 2) %>%
        pull(str_c("V", weights$target[i])),
    ) %>%
      mutate(
        weight = weights$weight[i],
        source = weights$source[i],
        target = weights$target[i]
      ) %>%
      unite(pair, c("source", "target"), remove = FALSE)
  }

  scatter_data <- bind_rows(scatter_data) %>%
    mutate(pair = reorder_within(pair, weight, source))
}

equivalence_similarity <- function(beta_hats_mat, B_tilde, beta_hats) {
  beta_hats <- beta_hats %>%
    pivot_wider(m:k_LDA, names_from = "w", values_from = "b") %>%
    rename(topic = k_LDA) %>%
    select(K, topic)

  cosine_similarity(beta_hats_mat, do.call(rbind, B_tilde)) %>%
    as_tibble() %>%
    bind_cols(beta_hats) %>%
    pivot_longer(starts_with("V"), names_to = "beta", values_to = "sim") %>%
    mutate(
      beta = str_replace(beta, "V", ""),
      beta = str_pad(beta, 2)
    ) %>%
    filter(K < max(K))
}

beta_hat_matrix <- function(alignment, beta_hats) {
  ordered_topics <- alignment %>%
    select(m, k_LDA, topic) %>%
    unique()

  beta_hats %>%
    left_join(ordered_topics) %>%
    select(m, K, k_LDA, topic, w, b) %>%
    pivot_wider(m:topic, names_from = w, values_from = b) %>%
    select(matches("[0-9]+")) %>%
    as.matrix()
}

widen_similarity <- function(similarity, K_ = 6) {
  similarity %>%
    filter(K == K_) %>%
    pivot_wider(topic, names_from = beta, values_from = sim) %>%
    select(-topic) %>%
    as.matrix()
}
