
#' Tree with one split per level
generate_tree <- function(n_levels = 6, edge_lengths = rep(1, 6)) {
  result <- list()
  result[[1]] <- tibble(m = 1, m_next = 2, k = 1, k_next = 1:2)

  for (k in seq(n_levels - 1)) {
    conserved <- tibble(m = k + 1, m_next = k + 2, k = seq_len(k + 1))
    branched <- tibble(m = k + 1, m_next = k + 2, k = sample(k + 1, 1))
    result[[k + 1]] <- bind_rows(conserved, branched) %>%
      arrange(k) %>%
      mutate(k_next = row_number())
  }

  bind_rows(result) %>%
    mutate(weight = 1 / m_next) %>%
    mutate(across(starts_with("m"), as.factor))
}

tree_topics <- function(tree, V=10, sigmas=NULL) {
  mu <- list()
  mu[["1"]] <- matrix(0, 1, V) # simulate root
  if (is.null(sigmas)) {
    sigmas <- rep(1, n_distinct(tree$m))
  }

  for (i in seq_len(nrow(tree))) {
    m_star <- tree$m[i]
    m_next <- as.character(tree$m_next[i])

    if (!(m_next %in% names(mu))) {
      K <- tree %>%
        filter(m == m_star) %>%
        slice_max(k_next) %>%
        pull(k_next)
      mu[[m_next]] <- matrix(0, K, V)
    }
    mu[[m_next]][tree$k_next[i], ] <- child_topic(mu[[m_star]][tree$k[i], ], sigmas[m_star])
  }

  beta <- map(mu, ~ exp(.) / rowSums(exp(.)))
  map_dfr(beta, as_tibble, .id = "m") %>%
    group_by(m) %>%
    mutate(k = row_number()) %>%
    rename_at(vars(starts_with("V")), ~ gsub("V", "", .)) %>%
    select(m, k, everything()) %>%
    ungroup()
}

child_topic <- function(mu, sigma=1) {
  V <- length(mu)
  mu + rnorm(V, 0, sigma)
}

widen_betas <- function(betas) {
  betas %>%
    pivot_wider(m:k_LDA, names_from = "w", values_from = "b") %>%
    mutate(k = as.numeric(as.factor(k_LDA))) %>%
    select(-K, -k_LDA)
}

tree_links <- function(edges, betas) {
  result <- list()
  for (i in seq_len(nrow(edges))) {
    source <- betas %>%
      filter(m == edges$m[i], k == edges$k[i]) %>%
      select(-m, -k) %>%
      rename_all(~ str_c("source_", .))
    target <- betas %>%
      filter(m == edges$m_next[i], k == edges$k_next[i]) %>%
      select(-m, -k) %>%
      rename_all(~ str_c("target_", .))
    result[[i]] <- bind_cols(source, target)
  }

  edges %>%
    bind_cols(bind_rows(result))
}

#' Simulate with No Topics
#'
#' n documents each with own topic uniform on V-dimensional simplex
uniform_simplex <- function(n, V = 500, n0 = NULL, lambda = 1) {
  if (is.null(n0)) {
    n0 <- rpois(n, 1000)
  }

  x <- matrix(nrow = n, ncol = V)
  colnames(x) <- seq_len(V)
  rownames(x) <- seq_len(n)
  for (i in seq_len(n)) {
    beta_i <- rdirichlet(1, rep(lambda, V))
    x[i, ] <- rmultinom(1, n0[i], beta_i)
  }

  x
}

#' Simulate Topics along V-Dimensional Curve
topic_curve <- function(K, k_anchor, V = 500, lambda = 1) {
  stopifnot(K > k_anchor)

  anchor_topics <- rdirichlet(k_anchor, rep(lambda, V))
  topics <- matrix(nrow = K, ncol = V)
  interpolation <- seq(1, k_anchor, length.out = K)
  for (k in seq_len(K - 1)) {
    k_prev <- floor(interpolation[k])
    lambda <- interpolation[k] - k_prev
    topics[k, ] <- (1 - lambda) * anchor_topics[k_prev, ] + lambda * anchor_topics[k_prev + 1, ]
  }
  topics[K, ] <- anchor_topics[k_anchor, ]
  flatten_beta(topics)
}

flatten_beta <- function(beta) {
  beta %>%
    as_tibble() %>%
    mutate(K = row_number(), node_id = K, depth = 1) %>%
    pivot_longer(-K:-depth, names_to = "w", values_to = "b") %>%
    mutate(w = as.integer(str_extract(w, "[0-9]+")))
}

#' Simulate Standard K-Topics Model
k_topics <- function(K, V = 500, lambda = 1) {
  rdirichlet(K, rep(lambda, V)) %>%
    flatten_beta()
}


normalize_mu <- function(mu) {
  beta <- matrix(0, nrow(mu), ncol(mu))
  for (k in seq_len(nrow(mu))) {
    beta[k, ] <- exp(mu[k, ])
    beta[k, ] <- beta[k, ] / sum(beta[k, ])
  }
  beta
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

embed_edges <- function(edges, ...) {
  umap_recipe <- recipe(~ ., edges) %>%
    update_role(replicate, m, m_next, k, k_next, weight, method, estimate, new_role = "id") %>%
    step_umap(all_predictors(), ...)

  prep(umap_recipe)
}

endpoint_topics <- function(betas, alignment) {
  beta_hats <- widen_betas(betas)
  alignment[c("beta_alignment", "gamma_alignment")] %>%
    map(~ select(., m, m_next, k, k_next, weight)) %>%
    map_dfr(~ tree_links(., beta_hats), .id = "method")
}

.simulate_replicate <- function(betas, gammas, M = 8, mechanism = NULL) {
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

fit_wrapper <- function(x, K_max) {
  fits <- run_lda_models(x, 1:K_max, c(.1), "VEM", 123, reset = TRUE)
  list(fits = fits, alignment = align_topics(fits))
}

vis_wrapper <- function(x, K_max) {
  fits <- run_lda_models(x, 1:K_max, c(.1), "VEM", 123, reset = TRUE)
  alignment <- align_topics(fits)
  list(
    fits = fits,
    alignment = alignment,
    p1 = visualize_aligned_topics(alignment),
    p2 = visualize_aligned_topics(alignment, method = "beta_alignment")
  )
}
