


align_topics = function(lda_models, m_ref = NULL, order_constrain = NULL){

  # 1. align topics based on the gamma matrices
  aligned_topics_gamma = align_sequence(
    lda_models$gammas,
    weight_fun = gamma_similarities
  )

  # 2. re-order topics based on the topic alignment
  topics_order = .order_topics(aligned_topics_gamma, order_constrain)
  aligned_topics_gamma = .join_order(aligned_topics_gamma, topics_order)

  # 3. align topics based on the beta matrices if possible
  aligned_topics_beta = align_sequence(
    lda_models$gammas,
    lda_models$betas,
    transport_similarities
  ) %>%
    .join_order(topics_order)

  list(
    lda_models = lda_models,
    gamma_alignment = aligned_topics_gamma,
    topics_order = topics_order,
    beta_alignment = aligned_topics_beta
  )
}

.order_topics = function(aligned_topics, order_constrain = NULL){

    if(is.null(order_constrain)){
      M = levels(aligned_topics$m)
      M = M %>% factor(., levels = M)

      order_constrain =
        aligned_topics %>%
        filter(m == M[1]) %>%
        select(m, k_LDA) %>%
        distinct() %>%
        arrange(k_LDA) %>%
        mutate(k = row_number())
    }else{
      # check provided order_constrain
      # 1. should have columns m, k_LDA, k
      # 2. m should be unique and one of M
      # 3. k_LDA should be distinct and matching the k_LDA of aligned_topics for the same m
      # 4. k should be distinct integers from 1:n(k_LDA)
    }



    # this function computes 3 features for each topic
    # 1. a topic order (k) which give a "rank" to each topic ranging from 1 to K for each model. If an order constrain is provided, then the computed topic order takes that into account.
    # 2. a branch, which reflects the emergence of new topics as the number of topic increases (branch). Most aligned topics are on the same branch. There might be two or more topics on a same branch for a given model.
    # 3. a topic name (topic) which reflects the alignment across models (k) (most aligned topics have the same or related topic names across models)

    t1 = .compute_topic_order(aligned_topics, order_constrain)
    t2 = .define_branches(aligned_topics, t1)
    t3 = .identify_topic_name(t2)

    topic_order = t3 %>%  select(m, k_LDA, branch, topic, k) %>%  arrange(m, branch, topic)
    topic_order
  }


.compute_topic_order = function(aligned_topics, order_constrain){

  downward_order = .compute_topic_order_downward(order_constrain = order_constrain, aligned_topics = aligned_topics)
  upward_order = .compute_topic_order_upward(order_constrain = order_constrain, aligned_topics = aligned_topics)
  ordered_topics =
    bind_rows(downward_order,
              upward_order) %>%
    distinct() %>%
    arrange(m, k_LDA, k) %>%
    select(m, k_LDA, k)

  ordered_topics
}


.compute_topic_order_upward = function(order_constrain, aligned_topics){

  ordered_topics = order_constrain
  M = levels(aligned_topics$m) %>%  factor(., levels = levels(aligned_topics$m))
  j = which(ordered_topics$m[1] == M)
  if(j == length(M)) return(ordered_topics)
  M_to_align = M[(j+1):length(M)]

  for(this_m in M_to_align){

    this_m = this_m %>% factor(., levels = levels(M))
    prev_m = next_level(this_m, n = -1)

    these_topics =
      aligned_topics %>%
      filter(m == prev_m) %>%
      left_join(.,
                ordered_topics %>%
                  filter(m == prev_m),
                by = c("m","k_LDA")
      )

    this_topic_order =
      these_topics %>%
      arrange(k_LDA_next, -weight) %>%
      group_by(k_LDA_next) %>%
      slice_head() %>%
      arrange(k, -weight) %>%
      ungroup() %>%
      select(k_LDA_next) %>%
      dplyr::rename(k_LDA = k_LDA_next) %>%
      mutate(m = this_m,
             k = row_number())

    ordered_topics =
      bind_rows(ordered_topics, this_topic_order)

  }

  ordered_topics
}


.compute_topic_order_downward = function(order_constrain, aligned_topics){

  ordered_topics = order_constrain
  M = levels(aligned_topics$m) %>%  factor(., levels = levels(aligned_topics$m))
  j = which(ordered_topics$m[1] == M)
  if(j == 1) return(ordered_topics)
  M_to_align = M[(j-1):1]

  for(this_m in M_to_align){

    this_m = this_m %>% factor(., levels = levels(M))
    next_m = next_level(this_m, n = 1)

    these_topics =
      aligned_topics %>%
      filter(m == this_m) %>%
      left_join(.,
                ordered_topics %>%
                  filter(m == next_m) %>%
                  mutate(m = this_m) %>%
                  dplyr::rename(k_LDA_next = k_LDA,
                                k_next = k),
                by = c("m","k_LDA_next")
      )

    this_topic_order =
      these_topics %>%
      arrange(k_LDA, -weight) %>%
      group_by(k_LDA) %>%
      slice_head(n = 1) %>%
      ungroup() %>%
      arrange(k_next, -weight) %>%
      mutate(k = k_next) %>%
      group_by(k) %>%
      mutate(k = k + (row_number()-1)/n())

    this_topic_order =
      this_topic_order %>%
      select(m, k_LDA, k)

    ordered_topics =
      bind_rows(ordered_topics, this_topic_order)

  }

  ordered_topics
}

.join_order <- function(aligned_topics, topics_order) {
  aligned_topics %>%
    left_join(.,
              topics_order,
              by = c("m", "k_LDA")) %>%
    left_join(.,
              topics_order %>%
                dplyr::rename(k_next = k,
                              k_LDA_next = k_LDA,
                              m_next = m,
                              branch_next = branch,
                              topic_next = topic),
              by = c("m_next", "k_LDA_next"))
}

.define_branches = function(aligned_topics, topic_order){

  M = levels(aligned_topics$m) %>%  factor(., levels = levels(aligned_topics$m))
  M = rev(M)

  branches =
    topic_order %>%
    filter(m == M[1]) %>%
    dplyr::rename(branch = k) %>%
    select(m, k_LDA, branch)

  for(this_m in M[-1]){

    this_m = this_m %>% factor(., levels = levels(M))
    next_m = next_level(this_m, n = 1)

    these_topics =
      aligned_topics %>%
      filter(m == this_m) %>%
      left_join(.,
                branches %>%
                  filter(m == next_m) %>%
                  mutate(m = this_m) %>%
                  dplyr::rename(k_LDA_next = k_LDA,
                                branch_next = branch),
                by = c("m","k_LDA_next")
      )

    this_model_branch =
      these_topics %>%
      arrange(k_LDA, -weight) %>%
      group_by(k_LDA) %>%
      slice_head(n = 1) %>%
      ungroup() %>%
      arrange(branch_next, -weight) %>%
      mutate(branch = branch_next)

    this_model_branch =
      this_model_branch %>%
      select(m, k_LDA, branch)

    branches =
      bind_rows(branches, this_model_branch)

  }

  branches %>%
    left_join(., topic_order, by = c("m","k_LDA"))

}

.identify_topic_name = function(topic_branch_and_order){

  topic_names =
    topic_branch_and_order %>%
    arrange(m, branch, k) %>%
    group_by(m, branch) %>%
    mutate(n_topic_in_branch = n(),
           topic_nb_in_branch = row_number(),
           topic =
             ifelse(n_topic_in_branch == 1,
                    branch %>% as.character(),
                    str_c(branch, letters[topic_nb_in_branch]))
             ) %>%
    select(m, k_LDA, branch, k, topic) %>%
    ungroup()

  topic_names
}

.compute_topic_number = function(topic_branch_and_order, order_constrain){

  M = topic_branch_and_order$m %>% levels()
  M = M %>% factor(., levels = M)

  max_n_per_branch =
    topic_branch_and_order %>%
    group_by(branch, m) %>%
    summarize(n = n(), .groups = "drop") %>%
    group_by(branch) %>%
    summarize(n = max(n), .groups = "drop")

  topic_number =
    max_n_per_branch[rep(1:nrow(max_n_per_branch), max_n_per_branch$n),] %>%
    mutate(min_k = row_number()) %>%
    group_by(branch) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    select(branch, min_k) %>%
    left_join(.,
              topic_branch_and_order %>% filter(m == last(M)),
              by = "branch") %>%
    mutate(k = pmax(l, min_k)) %>%
    select(m, k_LDA, l, branch, k)

  if(unique(order_constrain$m) ==  unique(topic_number$m)){
    order_constrain_tmp = order_constrain %>% arrange(k_LDA)
    topic_number_tmp = topic_number %>% arrange(k_LDA)
    if(!all(order_constrain_tmp$k == topic_number_tmp$k))
      warning("The topic order requested with the provided `order_constrain` is respected but the topic number had to be modified to account for the alignement structure.\n")
  }

  for(this_m in rev(M)[-1]){
    this_m = this_m %>% factor(., levels = levels(M))
    next_m = next_level(this_m, n = 1)

    this_m_topic_number =
      topic_branch_and_order %>%
      filter(m == this_m) %>%
      left_join(.,
                topic_number %>%
                  filter(m == next_m) %>%
                  select(branch, k) %>%
                  arrange(branch, k) %>%
                  group_by(branch) %>% slice_head(n = 1) %>% ungroup() %>%
                  dplyr::rename(k_next = k),
                by = "branch") %>%
      arrange(l) %>%
      group_by(branch) %>%
      mutate(k = k_next + row_number()-1) %>%
      select(m, k_LDA, l, branch, k)

    topic_number =
      bind_rows(topic_number, this_m_topic_number)

  }

  topic_number
}


next_level = function(f, n = 1){
  factor_levels = levels(f)
  m = match(f, factor_levels)
  new_m = ((m + n -1) %% length(factor_levels)) + 1
  new_f = factor_levels[new_m] %>% factor(., levels = factor_levels)
  new_f
}

.lengthen_weights <- function(A) {
  test <- A %>%
    rownames_to_column("k_LDA") %>%
    pivot_longer(-k_LDA, names_to = "k_LDA_next", values_to = "weight") %>%
    mutate(k_LDA_next = str_replace(k_LDA_next, "X", ""))
}

gamma_similarities <- function(gammas, ...) {
  A <- t(gammas[[1]]) %*% gammas[[2]]
  dimnames(A) <- map(gammas, ~ colnames(.))
  data.frame(A) %>%
    .lengthen_weights()
}

transport_similarities <- function(gammas, betas, reg = 1e-3, ...) {
  B <- do.call(cbind, betas)
  C <- philentropy::JSD(t(B))
  ix <- seq_len(ncol(betas[[1]]))

  A <- matrix(0, nrow = ncol(betas[[1]]), ncol = ncol(betas[[2]]))
  for (i in seq_len(nrow(gammas[[1]]))) {
    a <- t(gammas[[1]][i,, drop=F])
    b <- t(gammas[[2]][i,, drop=F])
    A <- A + Barycenter::Sinkhorn(a, b, C[ix, -ix, drop=F], lambda = reg)$Transportplan
  }

  dimnames(A) <- map(gammas, ~ colnames(.))
  data.frame(A) %>%
    .lengthen_weights()
}

align_sequence <- function(gamma_hats, beta_hats, weight_fun) {
  weights <- list()
  for (i in seq_along(gamma_hats)) {
    if (i == length(gamma_hats)) break
    weights[[i]] <- weight_fun(gamma_hats[i:(i + 1)], beta_hats[i:(i + 1)]) %>%
      mutate(m = names(gamma_hats)[i], m_next = names(gamma_hats)[i + 1])
  }
  postprocess_weights(weights, nrow(gamma_hats[[1]]))
}

postprocess_weights <- function(weights, n_docs) {
  bind_rows(weights) %>%
    group_by(m, m_next, k_LDA, k_LDA_next) %>%
    summarize(document_mass = sum(weight), .groups = "drop") %>%
    mutate(weight = document_mass / n_docs) %>%
    group_by(m_next, k_LDA_next) %>%
    mutate(norm_weight = weight / sum(weight)) %>%
    ungroup() %>%
    mutate(across(c("m", "m_next"), as.factor))
}

split_fits <- function(x, parameter = "gammas") {
  x <- x %>%
    split(.$m) %>%
    map(~ select(., -m, -K))

  if (parameter == "gammas") {
    x <- x %>%
      map(~ rename(., doc = d)) %>%
      map(~ pivot_wider(., doc, names_from = k_LDA, values_from = g)) %>%
      map(~ column_to_rownames(., "doc") %>% as.matrix())
  } else {
    x <- x %>%
      map(~ rename(., word = w)) %>%
      map(~ pivot_wider(., word, names_from = k_LDA, values_from = b)) %>%
      map(~ column_to_rownames(., "word") %>% as.matrix())
  }

  x
}
