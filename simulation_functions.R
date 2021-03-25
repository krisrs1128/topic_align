
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

agglomorate_gamma <- function(G, gamma, node_id = 1) {
  descendants <- G %>%
    filter(
      is.finite(node_distance_to(node == node_id, mode = "out")),
      leaf
    ) %>%
    pull(node) %>%
    as.character()
  tibble(
    d = seq_len(nrow(gamma)),
    g = rowSums(gamma[, descendants, drop=F])
  )
}

agglomorate_gammas <- function(tree, gamma) {
  G <- tbl_graph(
    tree %>% select(node, m, leaf), 
    tree %>% filter(!is.na(parent)) %>% select(-m, -leaf)
  )
  
  nodes <- pull(G, node)
  gamma_glom <- map(seq_along(nodes), ~ agglomorate_gamma(G, gamma, nodes[.]))
  names(gamma_glom) <- nodes
  bind_rows(gamma_glom, .id = "node") %>%
    mutate(node = as.integer(node)) %>%
    left_join(tree) %>%
    select(m, K, k_LDA, d, g)
}

reshape_tree_beta <- function(tree, topics) {
  cbind(tree, topics) %>%
    pivot_longer(matches("[0-9]+"), "w", values_to = "b") %>%
    select(m, K, k_LDA, w, b, k_LDA_)
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

simulation_error <- function(n_lev, lambda, layers = c(2, 3)) {
  tree <- generate_tree(n_lev, lambda)
  topics <- tree_topics(tree)
  sources <- beta_list(cbind(tree, topics), layers[1])
  targets <- beta_list(cbind(tree, topics), layers[2])
  alignment <- transport_align_pair(sources, targets, lambda = 0.01)
  
  tree_layer <- tree %>%
    filter(m %in% layers)
  
  alignment %>%
    mutate(parent = as.integer(k_LDA), node = as.integer(k_LDA_next)) %>%
    full_join(tree_layer) %>%
    mutate(
      time_delta = time - parent_time,
      time_delta = ifelse(is.na(time_delta), 0, time_delta),
      lambda = lambda,
      linked = time_delta != 0
    ) %>%
    select(lambda, parent, node, time_delta, linked, w)
}