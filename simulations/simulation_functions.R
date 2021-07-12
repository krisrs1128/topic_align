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
  topics <- MCMCpack::rdirichlet(K, rep(lambdas$beta, V))
  gammas <- MCMCpack::rdirichlet(N, rep(lambdas$gamma, K))
  n <- nrow(gammas)
  if (is.null(n0)) {
    n0 <- rpois(n, lambdas$count)
  }

  x <- matrix(nrow = n, ncol = ncol(topics))
  for (i in seq_len(n)) {
    nu_i <- matrix(
      MCMCpack::rdirichlet(1, rep(lambdas$nu, ncol(topics))),
      ncol = 1
    )
    mu <- alpha * t(topics) %*% gammas[i, ] + (1 - alpha) * nu_i
    x[i, ] <- rmultinom(1, n0[i], mu)
  }

  rownames(x) <- seq_len(n)
  colnames(x) <- seq_len(ncol(topics))
  list(topics = topics, gammas = gammas, x = x)
}

cosine_similarity <- function(x, y) {
  sim <- matrix(nrow = nrow(x), ncol = nrow(y))
  for (i in seq_len(nrow(x))) {
    for (j in seq_len(nrow(y))) {
      sim[i, j] <- sum(x[i, ] * y[j, ]) /
        sqrt(sum(x[i, ] ^ 2) * sum(y[j, ] ^ 2))
    }
  }
  sim
}

#' Simulation functions for functional equivalence experiment
perturb_topics <- function(topics, n_per_topic = NULL, subset_size = 40,
  alpha = 0.1) {
  K <- nrow(topics)
  if (is.null(n_per_topic)) {
    n_per_topic <- c(rep(2, 2), rep(1, K - 2))
  }

  perturbed_topics <- list()
  for (k in seq_len(K)) {
    perturbed_topics[[k]] <- rep(1, n_per_topic[k]) %*% matrix(topics[k, ], nrow = 1)
    if (n_per_topic[k] == 1) next

    for (i in seq_len(n_per_topic[k])) {
      ix <- seq_len(subset_size)
      perturbed_topics[[k]][i, ] <- perturb_topic(perturbed_topics[[k]][i, ], ix, alpha)
    }
  }

  perturbed_topics
}

perturb_topic <- function(topic_k, species_ix, alpha = 0.1) {
  alpha <- rep(alpha, length(species_ix))
  sub_community <- MCMCpack::rdirichlet(1, alpha)
  topic_k[species_ix] <- sum(topic_k[species_ix]) * sub_community
  topic_k
}

sample_topics <- function(perturbed_topics) {
  topics <- matrix(
    nrow = length(perturbed_topics),
    ncol = ncol(perturbed_topics[[1]])
  )
  for (k in seq_along(perturbed_topics)) {
    ix <- sample(seq_len(nrow(perturbed_topics[[k]])), 1)
    topics[k, ] <- perturbed_topics[[k]][ix, ]
  }
  topics
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

my_theme <- function() {
  th <- ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "#f7f7f7"),
      panel.border = ggplot2::element_rect(
        fill = NA, color = "#0c0c0c", size = 0.6
      ),
      legend.position = "bottom"
  )
  ggplot2::theme_set(th)
}

save_str <- function(prefix, id_vars, suffix = "csv") {
  out_dir <- id_vars$out_dir
  id_vars$out_dir <- NULL
  s <- stringr::str_c(paste0(id_vars, collapse = "_"), ".", suffix)
  file.path(out_dir, stringr::str_c(prefix, "-", s))
}
