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

simulate_gradient <- function(N, K, V, lambdas, alpha=1, n0=NULL) {
  B <- MCMCpack::rdirichlet(K, rep(lambdas$beta, V))
  gammas <- MCMCpack::rdirichlet(N, rep(lambdas$gamma, K))
  n <- nrow(gammas)
  if (is.null(n0)) {
    n0 <- rpois(n, lambdas$count)
  }

  x <- matrix(nrow = n, ncol = ncol(B))
  for (i in seq_len(n)) {
    nu_i <- matrix(MCMCpack::rdirichlet(1, rep(lambdas$nu, ncol(B))), ncol = 1)
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
      #ix <- c(1:(subset_size / 2), (subset_size / 2 + 1):subset_size)
      ix <- seq_len(subset_size)
      B_tilde[[k]][i, ] <- perturb_topic(B_tilde[[k]][i, ], ix, alpha)
    }
  }

  B_tilde
}

perturb_topic <- function(Bk, species_ix, alpha = 0.1) {
  alpha <- rep(alpha, length(species_ix))
  sub_community <- MCMCpack::rdirichlet(1, alpha)
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

simulate_document <- function(topics, gamma_i, n0) {
  z <- rmultinom(1, n0, gamma_i)
  x <- vector(length = ncol(topics))
  for (k in seq_along(z)) {
    x <- x + rmultinom(1, z[k], topics[k, ])
  }
  x
}

equivalence_data <- function(N, V, K, lambdas, subset_size, n0 = NULL) {
  if (is.null(n0)) {
    n0 <- rpois(N, lambdas$count)
  }

  topics <- MCMCpack::rdirichlet(K, rep(lambdas$beta, V))
  perturbed_topics <- perturb_topics(
    topics,
    subset_size = subset_size,
    alpha = lambdas$alpha
  )
  gammas <- MCMCpack::rdirichlet(N, rep(lambdas$gamma, K))
  x <- matrix(nrow = N, ncol = V)

  for (i in seq_len(N)) {
    local_topic <- sample_topics(perturbed_topics)
    x[i, ] <- simulate_document(local_topic, gammas[i, ], n0[i])
  }

  colnames(x) <- seq_len(ncol(x))
  rownames(x) <- seq_len(nrow(x))
  list(x = x, topics = topics, perturbed_topics = perturbed_topics,
       gammas = gammas)
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

fit_wrapper <- function(x, n_models, method) {
  lda_params <- purrr::map(1:n_models, ~ list(k = .))
  names(lda_params) <- stringr::str_c("K", 1:n_models)
  lda_models <- alto::run_lda_models(x, lda_params, reset = TRUE)
  alto::align_topics(lda_models, method = method)
}

bind_and_write <- function(scores, out_name, id) {
  dplyr::bind_rows(scores, .id = "alpha") %>%
    dplyr::mutate(alpha = alphas[as.integer(alpha)], id = id) %>%
    readr::write_csv(file = out_name)
}

my_theme <- function() {
  th <- ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "#f7f7f7"),
      panel.border = ggplot2::element_rect(fill = NA, color = "#0c0c0c", size = 0.6),
      legend.position = "bottom"
  )
  ggplot2::theme_set(th)
}

save_str <- function(prefix, id_vars, suffix = "csv") {
  out_dir <- id_vars$out_dir
  id_vars$out_dir <- NULL
  s <- stringr::str_c(paste0(id_vars, collapse = "_"), ".", suffix)
  file.path(out_dir, s)
}
