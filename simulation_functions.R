
generate_tree <- function(max_level = 4, lambda = 1) {
  tree <- data.frame(node = 1, parent = NA, leaf = TRUE, m = 0, time = 0, parent_time = NA)

  while (TRUE) {
    # each leaf has two children
    leaf_nodes <- which(tree$leaf)
    new_ids <- seq(max(tree$node) + 1, max(tree$node) + 1 + 2 * sum(leaf_nodes)) %>%
      matrix(ncol = 2, byrow = TRUE)

    # generate new leaves
    descendents <- list()
    for (i in seq_along(leaf_nodes)) {
      descendents[[i]] <- generate_leaves(tree[leaf_nodes[i], ], new_ids[i, ], lambda)
    }
    tree$leaf <- rep(FALSE, nrow(tree))
    tree <- bind_rows(descendents) %>%
      bind_rows(tree)

    # check stopping condition
    if (max(tree$m) >= max_level) break
  }

  tree %>%
    arrange(node) %>%
    group_by(m) %>%
    mutate(
      m = as.factor(m + 1),
      K = n(),
      k_LDA = letters[row_number()],
      k_LDA_ = letters[row_number()]
    ) %>%
    ungroup()
}

generate_leaves <- function(node, ids, lambda = 1) {
  result <- data.frame(node = ids, parent = node$node, leaf = TRUE, m = node$m + 1)
  result$parent_time <- node$time
  result$time <- node$time + rexp(length(ids), lambda)
  result
}

tree_topics <- function(tree, V=10) {
 mu <- matrix(0, nrow(tree), V)
 tree <- tree %>%
   filter(!is.na(parent)) # remove root

 for (i in seq_along(tree$node)) {
   time_delta <- tree$time[i] - tree$parent_time[i]
   mu[tree$node[i], ] <- mu[tree$parent[i], ] + rnorm(V, 0, sqrt(time_delta))
 }

 normalize_mu(mu)
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
  

#' Simulate topics on tree, repelling one another
repulsion_tree <- function(max_level = 4, V = 500, etas = NULL) {
  if (is.null(etas)) {
    etas <- seq(1, 0.1, length.out = max_level + 1)
  }

  G <- graph.tree(2 ^ max_level - 1) %>%
    as_tbl_graph() %>%
    mutate(
      depth = bfs_dist(1),
      id = bfs_rank()
    )
  node_ids <- pull(G, id)
  mu <- matrix(0, nrow = length(node_ids), ncol = V)

  # simulate leaves
  for (i in seq_along(node_ids)) {
    descendants <- G %>%
      activate(edges) %>%
      filter(from == i) %>%
      pull(to)
    if (length(descendants) == 0) break

    depth  <- G %>%
      filter(id == i) %>%
      pull(depth)

    epsilon <- rnorm(V, 0, etas[depth + 1])
    mu[descendants[1], ] <- mu[i, ] + epsilon / sqrt(sum(epsilon ^ 2))
    mu[descendants[2], ] <- mu[i, ] - epsilon / sqrt(sum(epsilon ^ 2))
  }

  data.frame(
    node_id = node_ids,
    depth = G %>% pull(depth),
    normalize_mu(mu)
  ) %>%
  pivot_longer(-node_id:-depth, names_to = "w", values_to = "b") %>%
    mutate(w = as.integer(str_replace(w, "X", "")))
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

beta_list <- function(tree_data, l = 2, prefix = "beta") {
  tree_data <- tree_data %>%
    filter(m == l)
  beta <- tree_data %>%
    select(matches("[0-9]+")) %>%
    t()
  colnames(beta) <- tree_data$node
  list(mass = rep(1 / ncol(beta), ncol(beta)), pos = beta)
}

embed_edges <- function(edges, ...) {
  umap_recipe <- recipe(~ ., edges) %>%
    update_role(replicate, K, k_LDA, k_LDA_next, edge_weight, method, new_role = "id") %>%
    step_umap(all_predictors(), ...)

  prep(umap_recipe)
}

endpoint_topics_ <- function(beta, edge_weights) {
  edge_weights <- edge_weights %>%
    select(m, m_next, k_LDA, k_LDA_next, norm_weight) %>%
    rename(edge_weight = norm_weight, K = m, K_next = m_next)
  
  wide_topics <- beta %>%
    select(K, k_LDA, w, b) %>%
    pivot_wider(K:k_LDA, w, values_from = b, names_prefix = "beta_") %>%
    mutate(K = as.factor(K))

  edge_weights %>%
    left_join(wide_topics) %>%
    left_join(wide_topics, c("K_next" = "K", "k_LDA_next" = "k_LDA")) %>%
    rename_with(~ gsub(".x$", "_source", .)) %>%
    rename_with(~ gsub(".y$", "_target", .))
}

endpoint_topics <- function(beta, alignment) {
  align1 <- endpoint_topics_(beta, alignment$gamma_alignment) %>%
    mutate(method = "gamma_alignment")
  align2 <- endpoint_topics_(beta, alignment$beta_alignment) %>%
    mutate(method = "beta_alignment")
  bind_rows(align1, align2)
}

simulate_leaves <- function(betas, gamma) {
  leaves <- betas %>%
    filter(depth == max(depth)) %>%
    pivot_wider(node_id, w, values_from = b) %>%
    select(-node_id)
  simulate_lda(leaves, gamma)
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