


align_topics = function(data, lda_models, m_ref = NULL, order_constrain = NULL){
  
  # adding the reference model
  if(is.null(m_ref)) m_ref = lda_models$gammas$m %>% levels() %>% last()
  lda_models$gammas = 
    lda_models$gammas %>% 
    mutate(m_ref = m_ref %>% 
             factor(.,levels = levels(lda_models$gammas$m))
           )
  
  # 1. align topics based on the gamma matrices
    
  aligned_topics_gamma = 
    .align_topics_gamma(
      gammas = lda_models$gammas)
  
  # 2. re-order topics based on the topic alignment
  topics_order = 
    .order_topics(
      aligned_topics = aligned_topics_gamma, order_constrain = order_constrain)
  
  # we add the new order for the topics to lda_models...
  lda_models$gammas = 
    lda_models$gammas %>% 
    left_join(., topics_order, by = c("m", "k_LDA"))
  lda_models$betas = 
    lda_models$betas %>% 
    left_join(., topics_order, by = c("m", "k_LDA"))
  
  aligned_topics_gamma = 
    aligned_topics_gamma %>% 
    left_join(., 
              topics_order, 
              by = c("m", "k_LDA")) %>% 
    left_join(., 
              topics_order %>% 
                dplyr::rename(k_next = k, k_LDA_next = k_LDA, m_next = m), 
              by = c("m_next", "k_LDA_next")) %>% 
    left_join(., 
              topics_order %>% 
                dplyr::rename(k_ref = k, k_LDA_ref = k_LDA, m_ref = m), 
              by = c("m_ref", "k_LDA_ref"))
    
  ans = list(
    lda_models = lda_models,
    gamma_alignment = aligned_topics_gamma,
    topics_order = topics_order
  )
  
  # 3. align topics based on the beta matrices if possible
  condition = FALSE
  if(condition){
    betas = lda_models$betas
    # betas = betas %>%  ... (add order)
    aligned_topics_beta = .align_topics_beta(data, betas = betas)
    ans$beta_alignment = aligned_topics_beta
  }
  
  ans
}




.align_topics_gamma = function(gammas){
  
  
  # 1. reference topic for each word
  t_d = 
    gammas %>%
    filter(m == m_ref) %>% 
    arrange(d, -g) %>% 
    group_by(d) %>% 
    slice_head(n = 1) %>% 
    select(m_ref, d, k_LDA, g) %>% 
    dplyr::rename(g_ref = g,
                  k_LDA_ref = k_LDA) %>% 
    ungroup()
    
  
  
  # 2a. Gamma_m
  
  gammas_m = 
    gammas %>% 
    select(-m_ref, -K) %>% 
    filter(m != levels(gammas$m) %>% last()) %>% 
    dplyr::rename(g_m = g) %>% 
    select(d, m, k_LDA, g_m)
  
  # 2b. Gamma_m_next
  
  gammas_m_next = 
    gammas %>% 
    select(-m_ref, -K) %>% 
    filter(m != levels(gammas$m)[1]) %>% 
    dplyr::rename(
      m_next = m,
      g_next = g,
      k_LDA_next = k_LDA) %>% 
    mutate(m = next_level(m_next, n = -1)) %>% 
    select(d, m, m_next, k_LDA_next, g_next) 
  
  # 3. weights
  
  W = 
    left_join(gammas_m,
              gammas_m_next, 
              by = c("m","d")) %>% 
    left_join(., t_d, by = "d") %>% 
    select(d, m, m_next, m_ref, k_LDA, k_LDA_next, k_LDA_ref,g_m, g_next, g_ref) %>% 
    mutate(w = g_m * g_next) %>% 
    group_by(m, m_next, m_ref, k_LDA, k_LDA_next, k_LDA_ref) %>% 
    summarize(w = sum(w), .groups = "drop") %>% 
    mutate(w = w / length(unique(gammas$d)))
  
  W
}

.order_topics = 
  function(
    aligned_topics,
    order_constrain = NULL
  ){
    
    M = levels(aligned_topics$m)
    M = M %>% factor(., levels = M)
    
    aligned_topics_summ = 
      aligned_topics %>% 
      group_by(m, k_LDA, k_LDA_next) %>% 
      summarize(weight = sum(w), .groups = "drop")
    
    if(is.null(order_constrain)){
      order_constrain = 
        aligned_topics_summ %>% 
        filter(m == M[1]) %>% 
        select(m, k_LDA) %>% 
        distinct() %>% 
        arrange(k_LDA) %>% 
        mutate(k = row_number())
    }else{
      # check provided order_constrain
      # 1. should have columns m, k_LDA, k
      # 2. m should be unique and one of M
      # 3. k_LDA should be distinct and matching the k_LDA of aligned_topic_summ for the same m
      # 4. k should be distinct integers from 1:n(k_LDA)
    }
    
    upward_order = .order_upward(order_constrain = order_constrain, aligned_topics = aligned_topics_summ)
    downward_order = .order_downward(order_constrain = order_constrain, aligned_topics = aligned_topics_summ)
    
    ordered_topics = 
      bind_rows(downward_order,
                upward_order) %>% 
      distinct() %>% 
      arrange(m, k_LDA, k)
    
    ordered_topics
    
  }



.order_upward = function(order_constrain, aligned_topics){
  
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



.order_downward = function(order_constrain, aligned_topics){
  
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
      arrange(-weight) %>% 
      mutate(is_duplicated = duplicated(k_next))
    
    if(any(this_topic_order$is_duplicated)){
      J = which(this_topic_order$is_duplicated)
      remaining_k_next_topics = 
        setdiff(
          these_topics$k_next %>% unique(),
          this_topic_order$k_next %>% unique()
          )  %>% sort()
      for(j in J){
        duplicated_topic = this_topic_order$k_next[j]
        closest_available_topic = 
          remaining_k_next_topics[which.min(abs(remaining_k_next_topics - duplicated_topic))]
        this_topic_order$k_next[j] = closest_available_topic
        remaining_k_next_topics = 
          setdiff(remaining_k_next_topics,
                  closest_available_topic)
      }
    }
     
    this_topic_order = 
      this_topic_order %>% 
      mutate(
        k = k_next
      ) %>% 
      select(m, k_LDA, k)
    
    ordered_topics =
      bind_rows(ordered_topics, this_topic_order)
    
  }
  
  ordered_topics
}



next_level = function(f, n = 1){
  factor_levels = levels(f)
  m = match(f, factor_levels)
  new_m = ((m + n -1) %% length(factor_levels)) + 1
  new_f = factor_levels[new_m] %>% factor(., levels = factor_levels)
  new_f
}





