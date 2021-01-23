

# Run LDA (topicmodels package) for various values of K

run_lda_models = 
  function(
    data,
    Ks = NULL, 
    alphas = NULL,
    method = c("VEM", "Gibbs"),
    seed = NULL, 
    dir = NULL,
    reset = FALSE
  ){
    
    # check data
    # TODO
    
    # check Ks
    if(any(Ks>26)) stop("Ks must contain integers smaller than 26") # TODO: change that
    
    # directory where each LDA model is temporally stored
    if(is.null(dir)){
      delete_dir = TRUE
      dir = ".tmp/"
    }else{
      delete_dir = FALSE
    }
    
    if(reset & dir.exists(dir)){
      unlink(dir, recursive = TRUE)
    }
    
    if(!dir.exists(dir)) dir.create(dir)
    
    # running LDA
    ok = map(
      .x = Ks[Ks>1],
      .f = function(K){
        topic_model_file = str_c(dir,"_LDA_K_", K, ".Rdata")
        if(!file.exists(topic_model_file)){
          topic_model = LDA(x = asv_for_topic, 
                            k = K, method = method, 
                            control = list(seed = seed))
          save(topic_model, file = topic_model_file)
        }
      }
    )
    
    # gammas
    gammas = 
      map_dfr(
        .x = Ks[Ks>1],
        .f = function(K){
          topic_model_file = str_c(dir,"_LDA_K_", K, ".Rdata")
          load(file = topic_model_file)
          
          documents_topic_composition = 
            topic_model@gamma %>% 
            as.data.frame() %>% 
            set_colnames(1:K) %>% 
            mutate(d = topic_model@documents) %>% 
            pivot_longer(cols = -d,
                         names_to = "k_LDA",
                         values_to = "g") %>% 
            mutate(K = K,
                   k_LDA = letters[k_LDA %>% as.integer()])
          documents_topic_composition 
        }
      )
    if(1 %in% Ks){
      gammas = 
        gammas %>% 
        bind_rows(
          tibble(
            K = 1, k_LDA = "a", 
            d = unique(gammas$d), 
            g = 1),
          .
        )
    }
    gammas = gammas %>% mutate(m = K %>% factor(., levels = Ks))
    gammas = gammas %>% select(m, K, k_LDA, d, g)
    
    
    # betas
    betas = 
      map_dfr(
        .x = Ks[Ks>1],
        .f = function(K){
          
          topic_model_file = str_c(dir,"_LDA_K_", K, ".Rdata")
          load(file = topic_model_file)
          
          topic_composition = 
            topic_model@beta %>% 
            exp() %>% 
            as.data.frame() %>% 
            set_colnames(topic_model@terms) %>% 
            mutate(k_LDA = letters[1:K]) %>% 
            pivot_longer(cols = -k_LDA,
                         names_to = "w",
                         values_to = "b") %>% 
            mutate(K = K)
          topic_composition 
        }
      )
    if(1 %in% Ks){
      data_mat = as.matrix(data)
      mean_prop = 
        (data_mat/rowSums(data_mat)) %>% 
        apply(., 2, mean) %>% 
        as.vector() %>%
        tibble(K = 1, k_LDA = "a", 
               w = colnames(data_mat),
               b = .)
      betas = bind_rows(mean_prop,betas)
    }
    betas = betas %>% mutate(m = K %>% factor(., levels = Ks))
    betas = betas %>% select(m, K, k_LDA, w, b)
    
    
    # delete dir (if dir was unspecified)
    if(delete_dir) unlink(dir)
    
    # return results
    res = list(
      betas = betas,
      gammas = gammas
    )
    res
  }


trim_models = 
  function(models,
           min_prop = NULL,
           n_words = NULL){
    
    trimmed_models = models
    
    if(is.null(min_prop) & is.null(n_words)) min_prop = 0.025
    if(!is.null(min_prop)){
      trimmed_models$betas = 
        trimmed_models$betas %>% 
        group_by(w) %>% 
        mutate(max_prop = max(b)) %>% 
        filter(max_prop >= min_prop) %>% 
        select(-max_prop) %>% 
        ungroup()
      return(trimmed_models)
    }else{
      top_n_words = 
        trimmed_models$betas %>% 
        group_by(w) %>% 
        summarize(max_prop = max(b),.groups = "drop") %>% 
        arrange(-max_prop) %>% 
        slice_head(n = n_words)
      trimmed_models$betas = 
        filter(w %in% top_n_words$w)
      return(trimmed_models)
    }
  }

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
      aligned_topics = aligned_topics_gamma) # , order_constrain = order_constrain
  
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
    aligned_topics
  ){
    
    M = levels(aligned_topics$m)
    M = M %>% factor(., levels = M)
    
    aligned_topics_summ = 
      aligned_topics %>% 
      group_by(m, k_LDA, k_LDA_next) %>% 
      summarize(weight = sum(w), .groups = "drop")

    ordered_topics = 
      aligned_topics_summ %>% 
      filter(m == M[1]) %>% 
      select(m, k_LDA) %>% 
      distinct() %>% 
      arrange(k_LDA) %>% 
      mutate(k = row_number())
    
    for(this_m in M[-1]){
        
        this_m = this_m %>% factor(., levels = levels(M))
        prev_m = next_level(this_m, n = -1)
        
        these_topics = 
          aligned_topics_summ %>% 
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




visualize_aligned_topics = function(aligned_topics){
  
  gamma_alignment = 
    aligned_topics$gamma_alignment %>% 
    select(-contains("_LDA"))
  
  m_ref = unique(gamma_alignment$m_ref)
  M = levels(gamma_alignment$m)
  
  h_m = 
    bind_rows(
      gamma_alignment %>% 
        select(m, k, w),
      gamma_alignment %>% 
        filter(m_next == M %>% last()) %>% 
        select(m_next, k_next, w) %>% 
        dplyr::rename(m = m_next, k = k_next)
      ) %>% 
    group_by(m) %>% 
    summarize(h = sum(w), .groups = "drop",
              K = max(k)) %>% 
    mutate(remaining_space = 2 - h,
           delta_k = remaining_space/(as.numeric(K)+1)) 
  
  
  layout_rect = 
    gamma_alignment %>% 
    select(m, k, k_ref, k_next, w) %>% 
    arrange(m, k, k_ref, k_next) %>% 
    # we add the last model
    bind_rows(
      .,
      gamma_alignment %>% 
        filter(m_next == M %>% last()) %>% 
        select(m_next, k_next, k_ref, w) %>% 
        group_by(m_next, k_next, k_ref) %>% 
        summarize(w = sum(w), .groups = "drop") %>% 
        mutate(m = m_next,
               k = k_next) %>% 
        select(m, k, k_ref, k_next, w)
      ) %>% 
    left_join(., h_m %>% select(m, h, delta_k), by = "m") %>% 
    group_by(m) %>% 
    mutate(cum_w = cumsum(w),
           ymin = -k*delta_k - cum_w ,
           ymax = ymin + w,
           y = (ymin + ymax)/2,
           height = w)
  
  
  layout_rect = 
    layout_rect %>% 
    mutate(k_ref = k_ref %>% factor(),
           m = m %>% factor(., levels = c(levels(gamma_alignment$m),"document")),
           m_num = m %>% as.numeric())
  
  
  OUT = 
    layout_rect %>% 
    filter(m != last(M)) %>% 
    select(m, k, k_ref, k_next, ymin, ymax, m_num) %>% 
    mutate(x = m_num + 0.2)
  
  IN = 
    layout_rect %>% 
    filter(m != last(M)) %>% 
    select(m, k, k_ref, k_next, m_num, w) %>% 
    arrange(m, k_next, k_ref, k) %>% 
    mutate(m_next = next_level(m, n = 1)) %>% 
    left_join(
      ., 
      h_m %>% 
        select(m, delta_k) %>% 
        dplyr::rename(m_next = m),
      by = "m_next") %>% 
    mutate(cum_w = cumsum(w),
           ymin = -k_next*delta_k - cum_w,
           ymax = ymin + w,
           x = m_num +1 - 0.2) %>% 
    select(m, k, k_ref, k_next, ymin, ymax, m_num, x)
  
  layout_ribbons = 
    bind_rows(OUT, IN) %>% 
    arrange(m, k, k_ref, k_next)
  
  layout_ribbons = 
    layout_ribbons %>% 
    mutate(h = ymax - ymin,
           flow_id = str_c("m_",m,"|k_",k,"|k_ref_",k_ref,"|k_next_",k_next)) 
    
    
  M_nums = unique(layout_rect$m_num)  %>%  sort()
  
  g = 
    ggplot(
      layout_rect,
      aes(fill = k_ref)) +
    geom_tile(
      aes(
        x = m_num, y = y,
        width = 0.4,
        height = height)
    ) +
    guides(fill = FALSE) +
    geom_ribbon(
      data = layout_ribbons,
      aes(x = x, ymin = ymin, ymax = ymax, group = flow_id, fill = k_ref),
      alpha = 0.5
    ) +
    scale_x_continuous(breaks = M_nums, minor_breaks = NULL, labels = M) +
    scale_y_continuous(breaks = NULL) +
    xlab("models") +
    ylab("topic composition")
  
  g
  
  }

next_level = function(f, n = 1){
  factor_levels = levels(f)
  m = match(f, factor_levels)
  new_m = ((m + n -1) %% length(factor_levels)) + 1
  new_f = factor_levels[new_m] %>% factor(., levels = factor_levels)
  new_f
}





