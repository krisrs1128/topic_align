

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
          topic_model = LDA(x = data,
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

align_topics = function(lda_models, m_ref = NULL, order_constrain = NULL){

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
    aligned_topics_beta = .align_topics_beta(betas = betas)
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
      slice_head() %>%
      arrange(k_next, -weight) %>%
      ungroup() %>%
      select(m, k_LDA) %>%
      mutate(k = row_number())

    ordered_topics =
      bind_rows(ordered_topics, this_topic_order)

  }

  ordered_topics
}



# .order_topics =
#   function(
#     aligned_topics
#   ){
#
#     M = levels(aligned_topics$m)
#     M = M %>% factor(., levels = M)
#
#     aligned_topics_summ =
#       aligned_topics %>%
#       group_by(m, k_LDA, k_LDA_next) %>%
#       summarize(weight = sum(w), .groups = "drop")
#
#     ordered_topics =
#       aligned_topics_summ %>%
#       filter(m == M[1]) %>%
#       select(m, k_LDA) %>%
#       distinct() %>%
#       arrange(k_LDA) %>%
#       mutate(k = row_number())
#
#     for(this_m in M[-1]){
#
#         this_m = this_m %>% factor(., levels = levels(M))
#         prev_m = next_level(this_m, n = -1)
#
#         these_topics =
#           aligned_topics_summ %>%
#           filter(m == prev_m) %>%
#           left_join(.,
#                     ordered_topics %>%
#                       filter(m == prev_m),
#                     by = c("m","k_LDA")
#           )
#
#         this_topic_order =
#           these_topics %>%
#           arrange(k_LDA_next, -weight) %>%
#           group_by(k_LDA_next) %>%
#           slice_head() %>%
#           arrange(k, -weight) %>%
#           ungroup() %>%
#           select(k_LDA_next) %>%
#           dplyr::rename(k_LDA = k_LDA_next) %>%
#           mutate(m = this_m,
#                  k = row_number())
#
#         ordered_topics =
#           bind_rows(ordered_topics, this_topic_order)
#
#       }
#
#     ordered_topics
#   }


visualize_aligned_topics =
  function(
    aligned_topics,
    color_by = c("alignment","reference"),
    add_leaves = TRUE,
    min_beta = 0.025,
    n_words = NULL,
    add_words_labels = TRUE,
    reverse_x_axis = FALSE
  ){

  color_by = color_by[1]

  if(add_leaves){
    if(!is.null(n_words))
      aligned_topics$lda_models = trim_models(models = aligned_topics$lda_models, n_words = n_words)
    else
      aligned_topics$lda_models = trim_models(models = aligned_topics$lda_models, min_prop = min_beta)
  }

  layouts = .compute_layout(aligned_topics = aligned_topics,
                            color_by = color_by,
                            add_leaves = add_leaves)


  layout_rect = layouts$rect
  layout_ribbons = layouts$ribbons


  if(add_leaves){
    leaves_layout = .compute_leaves_layout(aligned_topics)
    layout_rect =
      bind_rows(
        layout_rect %>% select(m, k_ref, m_num, y, height) %>% mutate(m = m %>% as.character()),
        leaves_layout$rect %>% select(m, k_ref, m_num, y, height)
        )
    layout_ribbons =
      bind_rows(
        layout_ribbons %>% ungroup() %>% select(x, k_ref, ymin, ymax, flow_id),
        leaves_layout$ribbons %>% select(x, k_ref, ymin, ymax, flow_id)
      )
  }


  M_nums = unique(layout_rect$m_num)  %>%  sort()
  M = unique(layout_rect$m)

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
    scale_y_continuous(breaks = NULL) +
    xlab("models") +
    ylab("topic composition")

  if(add_leaves & add_words_labels){
    longuest_word = leaves_layout$text$w %>% str_length() %>% max()
    g =
      g +
      geom_text(data = leaves_layout$text,
                aes(x = m_num + 0.35,
                    y = y,
                    label = w,
                    col = k_ref),
                hjust = ifelse(reverse_x_axis,1,0), size = 3) +
      guides(col = FALSE) +
      expand_limits(x = max(leaves_layout$text$m_num) + longuest_word/5)
      coord_cartesian(clip = 'off')
  }

  if(reverse_x_axis){
    g = g + scale_x_reverse(breaks = M_nums, minor_breaks = NULL, labels = M)
  }else{
    g = g + scale_x_continuous(breaks = M_nums, minor_breaks = NULL, labels = M)
  }

  g
}


.compute_layout = function(aligned_topics, color_by, add_leaves){

  gammas = aligned_topics$lda_models$gammas

  gamma_alignment =
    aligned_topics$gamma_alignment %>%
    select(-contains("_LDA"))

  m_ref = unique(gammas$m_ref)
  M = unique(gammas$m)

  h_m =
    gammas %>%
    group_by(m) %>%
    summarize(K = max(k), .groups = "drop") %>%
    mutate(delta_k = 1/(K+1))

  if(color_by == "reference")
    layouts = .compute_layout_reference(gamma_alignment = gamma_alignment,
                                        h_m = h_m, M = M, m_ref = m_ref,
                                        add_leaves = add_leaves)
  else
    layouts = .compute_layout_alignment(gammas = gammas, gamma_alignment = gamma_alignment,
                                        h_m = h_m, M = M, m_ref = m_ref,
                                        add_leaves = add_leaves)


  layouts
}


.compute_layout_alignment = function(gammas, gamma_alignment, h_m, M, m_ref, add_leaves){

  k_ref_table = .compute_k_ref_table(gamma_alignment = gamma_alignment)
  k_ref_levels = k_ref_table$k_ref %>% unique() %>%  sort()

  layout_rect =
    gammas %>%
    select(m, k, g) %>%
    group_by(m, k) %>%
    summarize(g = sum(g), .groups = "drop") %>%
    group_by(m) %>%
    mutate(w = g/sum(g)) %>%
    left_join(., k_ref_table, by = c("m", "k")) %>%
    left_join(., h_m %>% select(m, delta_k), by = "m") %>%
    group_by(m) %>%
    mutate(height = w,
           cum_w = cumsum(w),
           y = - k*delta_k - cum_w + w/2,
           m_num = m %>% as.integer(),
           k_ref = k_ref %>% factor(., levels = k_ref_levels))

  no_k_ref_gamma_alignment =
    gamma_alignment %>%
    group_by(m, m_next, k, k_next) %>%
    summarize(w = sum(w), .groups = "drop")

  OUT =
    no_k_ref_gamma_alignment %>%
    mutate(m_num = m %>% as.integer(),
           x = m_num + 0.2) %>%
    left_join(.,
              k_ref_table %>%
                dplyr::rename(m_next = m, k_next = k),
              by = c("m_next", "k_next")) %>%
    left_join(., h_m %>% select(m, delta_k), by = "m") %>%
    group_by(m) %>%
    mutate(cum_w = cumsum(w),
           ymin = - k*delta_k - cum_w,
           ymax = ymin + w)

  IN =
    no_k_ref_gamma_alignment %>%
    mutate(m_num = m %>% as.integer(),
           x = m_num + 1 - 0.2) %>%
    left_join(.,
              k_ref_table %>%
                dplyr::rename(m_next = m, k_next = k),
              by = c("m_next", "k_next")) %>%
    left_join(., h_m %>% select(m, delta_k) %>% dplyr::rename(m_next = m), by = "m_next") %>%
    arrange(m, m_next, k_next, k) %>%
    group_by(m) %>%
    mutate(cum_w = cumsum(w),
           ymin = -k_next*delta_k - cum_w,
           ymax = ymin + w)


  layout_ribbons =
    bind_rows(OUT, IN) %>%
    arrange(m, k, k_ref, k_next)

  layout_ribbons =
    layout_ribbons %>%
    mutate(h = ymax - ymin,
           flow_id = str_c("m_",m,"|k_",k,"|k_ref_",k_ref,"|k_next_",k_next),
           k_ref = k_ref %>% factor(., levels = k_ref_levels))


  list(rect = layout_rect, ribbons = layout_ribbons)
}

.compute_k_ref_table = function(gamma_alignment){

  M = levels(gamma_alignment$m)

  k_ref_table =
    gamma_alignment %>%
    select(m_next, k_next) %>%
    filter(m_next == last(M)) %>%
    distinct() %>%
    mutate(k_ref = k_next,
           k = k_next,
           m = m_next) %>%
    select(m, k, k_ref)

  for(this_m in rev(M)[-1]){
    this_m = this_m %>% factor(., levels = M)
    this_k_ref_table =
      gamma_alignment %>%
      select(m, k, k_next, w) %>%
      filter(as.character(m) == this_m) %>%
      group_by(m,k, k_next) %>%
      summarize(w = sum(w), .groups = "drop") %>%
      arrange(k, -w) %>%
      group_by(k) %>%
      slice_head(n = 1) %>%
      left_join(.,
                k_ref_table %>%
                  filter(m == next_level(this_m)) %>%
                  dplyr::rename(k_next = k) %>%
                  select(-m),
                by = "k_next") %>%
      select(m, k, k_ref)

    k_ref_table =
      bind_rows(this_k_ref_table, k_ref_table)
  }
  k_ref_table =
    k_ref_table %>% arrange(m, k, k_ref)
  k_ref_table
}


.compute_leaves_layout = function(
  aligned_topics
){

  m_ref = unique(aligned_topics$lda_models$gammas$m_ref)

  adjusted_betas_ref =
    aligned_topics$lda_models$betas %>%
    select(m, k, w, b) %>%
    filter(m == m_ref) %>%
    group_by(m, k) %>%
    mutate(b = b/sum(b))

  ref_topic_for_each_word =
    adjusted_betas_ref %>%
    arrange(w, - b) %>%
    group_by(w) %>%
    slice_head(n = 1) %>%
    dplyr::rename(k_ref = k) %>%
    ungroup() %>%
    arrange(k_ref, -b) %>%
    mutate(w_index = row_number())

  ref_topic_weights =
    aligned_topics$lda_models$gammas %>%
    filter(m == m_ref) %>%
    group_by(k) %>%
    summarize(weight = sum(g), .groups = "drop") %>%
    mutate(weight = weight/sum(weight)) %>%
    dplyr::rename(k_ref = k)

  est_average_freq =
    adjusted_betas_ref %>%
    dplyr::rename(k_ref = k) %>%
    left_join(., ref_topic_weights, by = "k_ref") %>%
    mutate(f = b * weight) %>%
    group_by(w) %>%
    summarize(f = sum(f), .groups = "drop") %>%
    arrange(-f)


  # rectangles

  total_height = est_average_freq$f %>% sum()
  delta = (2-total_height)/(nrow(ref_topic_weights) + nrow(est_average_freq))

  leaves_layout_rect =
    ref_topic_for_each_word %>%
    left_join(., est_average_freq, by = "w") %>%
    mutate(m = "words",
           m_num = length(levels(aligned_topics$gamma_alignment$m))+1 + 1,
           cum_f = cumsum(f),
           y = - (k_ref + row_number() - 1)*delta - cum_f + f/2,
           height = f,
           k_ref = k_ref %>% factor())

  # Flows

  delta_out = 1/(nrow(ref_topic_weights)+1)

  leaves_OUT =
    ref_topic_weights %>%
    left_join(., adjusted_betas_ref %>% dplyr::rename(k_ref = k), by = "k_ref") %>%
    mutate(h = weight * b) %>%
    # filter(h > 0.001) %>%
    # mutate(h = h/sum(h)) %>%
    select(m, k_ref, w, h) %>%
    mutate(x = m_ref %>% as.numeric() + 0.2,
           cum_h = cumsum(h),
           ymin = - delta_out * k_ref - cum_h,
           ymax = ymin + h)


  fraction_of_word_from_each_topic =
    adjusted_betas_ref %>%
    ungroup() %>%
    group_by(w) %>%
    mutate(c = b/sum(b)) %>%
    arrange(w, k)

  leaves_IN =
    est_average_freq %>%
    left_join(., fraction_of_word_from_each_topic, by = "w") %>%
    dplyr::rename(k_ref = k) %>%
    # left_join(.,leaves_OUT %>% select(m, k_ref, w, h), by = c("w", "m", "k_ref")) %>%
    # filter(!is.na(h)) %>% select(-h) %>%
    mutate(h = f*c,
           h = h/sum(h)) %>%
    left_join(.,
              ref_topic_for_each_word %>%
                dplyr::rename(k_ref_w = k_ref) %>%
                select(w, k_ref_w, w_index),
              by = "w") %>%
    arrange(w_index, k_ref) %>%
    mutate(x = m_ref %>% as.numeric() + 0.8 + 1,
           cum_h = cumsum(h),
           ymin = - delta * (k_ref_w + w_index - 1) - cum_h,
           ymax = ymin + h)

  leaves_layout_ribbon =
    bind_rows(
      leaves_OUT %>% select(m, k_ref, w, x, ymin, ymax),
      leaves_IN %>% select(m, k_ref, w, x, ymin, ymax)
    ) %>%
    arrange(m, k_ref, w, x) %>%
    left_join(., ref_topic_for_each_word %>% select(w, w_index), by = "w") %>%
    mutate(flow_id = str_c("leaves|k_ref_",k_ref,"|w_",w_index)) %>%
    select(-w_index) %>%
    mutate(k_ref = k_ref %>% factor())

  # words

  leaves_layout_words =
    leaves_layout_rect %>%
    select(m_num, y, k_ref, w)

  # return the layouts

  list(rect = leaves_layout_rect,
       ribbons = leaves_layout_ribbon,
       text = leaves_layout_words)
}




.compute_layout_reference = function(gamma_alignment, h_m, M, m_ref, add_leaves){

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
    left_join(., h_m %>% select(m, delta_k), by = "m") %>%
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
    filter(as.character(m) != last(M)) %>%
    select(m, k, k_ref, k_next, ymin, ymax, m_num) %>%
    mutate(x = m_num + 0.2)

  IN =
    layout_rect %>%
    filter(as.character(m) != last(M)) %>%
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

  list(rect = layout_rect, ribbons = layout_ribbons)
}



# visualize_aligned_topics = function(aligned_topics, color_by = c("reference","alignment")){
#
#   gamma_alignment =
#     aligned_topics$gamma_alignment %>%
#     select(-contains("_LDA"))
#
#   m_ref = unique(gamma_alignment$m_ref)
#   M = levels(gamma_alignment$m)
#
#   h_m =
#     bind_rows(
#       gamma_alignment %>%
#         select(m, k, w),
#       gamma_alignment %>%
#         filter(m_next == M %>% last()) %>%
#         select(m_next, k_next, w) %>%
#         dplyr::rename(m = m_next, k = k_next)
#       ) %>%
#     group_by(m) %>%
#     summarize(h = sum(w), .groups = "drop",
#               K = max(k)) %>%
#     mutate(remaining_space = 2 - h,
#            delta_k = remaining_space/(as.numeric(K)+1))
#
#
#   layout_rect =
#     gamma_alignment %>%
#     select(m, k, k_ref, k_next, w) %>%
#     arrange(m, k, k_ref, k_next) %>%
#     # we add the last model
#     bind_rows(
#       .,
#       gamma_alignment %>%
#         filter(m_next == M %>% last()) %>%
#         select(m_next, k_next, k_ref, w) %>%
#         group_by(m_next, k_next, k_ref) %>%
#         summarize(w = sum(w), .groups = "drop") %>%
#         mutate(m = m_next,
#                k = k_next) %>%
#         select(m, k, k_ref, k_next, w)
#       ) %>%
#     left_join(., h_m %>% select(m, h, delta_k), by = "m") %>%
#     group_by(m) %>%
#     mutate(cum_w = cumsum(w),
#            ymin = -k*delta_k - cum_w ,
#            ymax = ymin + w,
#            y = (ymin + ymax)/2,
#            height = w)
#
#
#   layout_rect =
#     layout_rect %>%
#     mutate(k_ref = k_ref %>% factor(),
#            m = m %>% factor(., levels = c(levels(gamma_alignment$m),"document")),
#            m_num = m %>% as.numeric())
#
#
#   if(color_by == "reference"){
#     OUT =
#       layout_rect %>%
#       filter(m != last(M)) %>%
#       select(m, k, k_ref, k_next, ymin, ymax, m_num) %>%
#       mutate(x = m_num + 0.2)
#
#     IN =
#       layout_rect %>%
#       filter(m != last(M)) %>%
#       select(m, k, k_ref, k_next, m_num, w) %>%
#       arrange(m, k_next, k_ref, k) %>%
#       mutate(m_next = next_level(m, n = 1)) %>%
#       left_join(
#         .,
#         h_m %>%
#           select(m, delta_k) %>%
#           dplyr::rename(m_next = m),
#         by = "m_next") %>%
#       mutate(cum_w = cumsum(w),
#              ymin = -k_next*delta_k - cum_w,
#              ymax = ymin + w,
#              x = m_num +1 - 0.2) %>%
#       select(m, k, k_ref, k_next, ymin, ymax, m_num, x)
#
#   }
#
#   if(color_by == "alignment"){
#     layout_rect =
#       layout_rect %>%
#       group_by(m, m_num, k, k_next) %>%
#       summarize(w = sum(w),
#                 ymin = min(ymin),
#                 ymax = max(ymax),
#                 height = sum(height),
#                 .groups = "drop") %>%
#       mutate(y = (ymin+ymax)/2)
#
#     k_ref_table =
#       layout_rect %>%
#       select(m, k) %>%
#       filter(m == last(M)) %>%
#       mutate(k_ref = k)
#
#     for(this_m in rev(M)[-1]){
#       this_m = this_m %>% factor(., levels = levels(layout_rect$m))
#       this_k_ref_table =
#       layout_rect %>%
#         select(m, k, k_next, w) %>%
#         filter(m == this_m) %>%
#         arrange(k, -w) %>%
#         group_by(k) %>%
#         slice_head(n = 1) %>%
#         left_join(.,
#                   k_ref_table %>%
#                     filter(m == next_level(this_m)) %>%
#                     dplyr::rename(k_next = k) %>%
#                     select(-m),
#                   by = "k_next") %>%
#         select(m, k, k_ref)
#
#       k_ref_table =
#         bind_rows(this_k_ref_table, k_ref_table)
#     }
#
#     layout_rect =
#       layout_rect %>%
#       left_join(., k_ref_table, by = c("m", "k")) %>%
#       mutate(k_ref = k_ref %>% factor())
#
#     OUT =
#       layout_rect %>%
#       filter(m != last(M)) %>%
#       select(m, k, k_ref, k_next, ymin, ymax, m_num) %>%
#       mutate(x = m_num + 0.2)
#
#     IN =
#       layout_rect
#
#
#   }
#
#
#   layout_ribbons =
#     bind_rows(OUT, IN) %>%
#     arrange(m, k, k_ref, k_next)
#
#   layout_ribbons =
#     layout_ribbons %>%
#     mutate(h = ymax - ymin,
#            flow_id = str_c("m_",m,"|k_",k,"|k_ref_",k_ref,"|k_next_",k_next))
#
#
#   M_nums = unique(layout_rect$m_num)  %>%  sort()
#
#   g =
#     ggplot(
#       layout_rect,
#       aes(fill = k_ref)) +
#     geom_tile(
#       aes(
#         x = m_num, y = y,
#         width = 0.4,
#         height = height)
#     ) +
#     guides(fill = FALSE) +
#     geom_ribbon(
#       data = layout_ribbons,
#       aes(x = x, ymin = ymin, ymax = ymax, group = flow_id, fill = k_ref),
#       alpha = 0.5
#     ) +
#     scale_x_continuous(breaks = M_nums, minor_breaks = NULL, labels = M) +
#     scale_y_continuous(breaks = NULL) +
#     xlab("models") +
#     ylab("topic composition")
#
#   g
#
#   }

next_level = function(f, n = 1){
  factor_levels = levels(f)
  m = match(f, factor_levels)
  new_m = ((m + n -1) %% length(factor_levels)) + 1
  new_f = factor_levels[new_m] %>% factor(., levels = factor_levels)
  new_f
}
