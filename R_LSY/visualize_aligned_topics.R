

visualize_aligned_topics =
  function(
    aligned_topics,
    add_leaves = TRUE,
    min_beta = 1e-3,
    n_words = NULL,
    add_words_labels = TRUE,
    reverse_x_axis = FALSE,
    min_k = NULL,
    max_k = NULL,
    method = "gamma_alignment"
  ){

    M = levels(aligned_topics$topics_order$m)
    M = M %>% factor(., levels = M)
    all_ks = aligned_topics$topics_order$branch %>% unique() %>% sort()
    min_k = ifelse(is.null(min_k), min(all_ks),min_k)
    max_k = ifelse(is.null(max_k), max(all_ks),max_k)
    if(min_k > min(all_ks)) stop("'min_k' must be equal or smaller than the first topic.\n")
    if(max_k < max(all_ks)) stop("'max_k' must be equal or larger than the last topic.\n")
    aligned_topics$lda_models$betas <- .lengthen_betas(aligned_topics$lda_models$betas, M) %>%
      left_join(aligned_topics$topics_order)
    aligned_topics$lda_models$gammas <- .lengthen_gammas(aligned_topics$lda_models$gammas, M) %>%
      left_join(aligned_topics$topics_order)

    if(add_leaves){
      if(!is.null(n_words))
        aligned_topics$lda_models = trim_models(models = aligned_topics$lda_models, n_words = n_words)
      else
        aligned_topics$lda_models = trim_models(models = aligned_topics$lda_models, min_prop = min_beta)
    }

    aligned_topics$alignment <- aligned_topics[[method]]

    # if there are "ghost topics", we need to add them
    aligned_topics_with_ghosts =
      .add_ghost_topics_if_needed(
        aligned_topics, min_k = min_k, max_k = max_k
      )


    # finally, we need to compute a topic index (for the y-position)
    aligned_topics = .compute_topic_index(aligned_topics_with_ghosts)

    # compute layout translates the information of topic weights and flow into x, y, height, ymin, ymax, color
    layouts = .compute_layout(aligned_topics = aligned_topics,
                              min_k = min_k,
                              max_k = max_k)

    layout_rect = layouts$rect
    layout_ribbons = layouts$ribbons

    if(add_leaves){
      leaves_layout = .compute_leaves_layout(aligned_topics, min_k = min_k, max_k = max_k)

      layout_rect =
        bind_rows(
          layout_rect %>% select(m, x, y, height, color) %>% mutate(m = m %>% as.character()),
          leaves_layout$rect %>% select(m, x, y, height, color)
        )
      layout_ribbons =
        bind_rows(
          layout_ribbons %>% ungroup() %>% select(x, ymin, ymax, color, flow_id),
          leaves_layout$ribbons %>% select(x, ymin, ymax, color, flow_id)
        )
    }


    M_nums = unique(layout_rect$x)  %>%  sort()
    M = unique(layout_rect$m)
    all_ks = min_k:max_k

    g =
      ggplot(
        layout_rect,
        aes(fill = color)) +
      geom_tile(
        aes(
          x = x, y = y,
          width = 0.4,
          height = height)
      ) +
      guides(fill = FALSE) +
      geom_ribbon(
        data = layout_ribbons,
        aes(x = x, ymin = ymin, ymax = ymax, group = flow_id, fill = color),
        alpha = 0.5
      ) +
      scale_y_continuous(breaks = NULL, limits = c(-2,0)) +
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
                      col = color),
                  hjust = ifelse(reverse_x_axis,1,0), size = 3) +
        guides(col = FALSE) +
        expand_limits(x = max(leaves_layout$text$m_num) + longuest_word/5) +
        coord_cartesian(clip = 'off')
    }

    if(reverse_x_axis){
      g = g + scale_x_reverse(breaks = M_nums, minor_breaks = NULL, labels = M)
    }else{
      g = g + scale_x_continuous(breaks = M_nums, minor_breaks = NULL, labels = M)
    }

    g = g +
      scale_fill_discrete(limits = all_ks %>% factor()) +
      scale_color_discrete(limits = all_ks %>% factor())

    g
  }




.add_ghost_topics_if_needed =
  function(aligned_topics, min_k, max_k){
    alignment = aligned_topics$alignment
    topic_order = aligned_topics$topics_order
    gammas = aligned_topics$lda_models$gammas
    betas = aligned_topics$lda_models$betas

    M = unique(topic_order$m) %>% sort()
    M = M %>% factor(., levels = M)


    # Do we need to add ghost topics?

    K_max = max_k - min_k + 1 # the total number of topics, incl the "ghost" topics
    all_ks = min_k:max_k
    actual_ks = topic_order$branch %>% unique() %>% sort()

    ghost_topics = setdiff(x = all_ks, y = actual_ks) %>% sort()
    K_g = length(ghost_topics) # the number of ghost topic to add

    # if there are no ghost topic, we don't need to do anything.
    if(K_g == 0){
      gammas$is_ghost = FALSE
      if("betas" %in% names(aligned_topics$lda_models)) betas$is_ghost = FALSE
      alignment$is_ghost = FALSE
      topic_order$is_ghost = FALSE
      lda_models = list(gammas = gammas)
      if("betas" %in% names(aligned_topics$lda_models)) lda_models$betas = betas
      return(
        list(lda_models = lda_models,
             alignment =  alignment,
             topic_order = topic_order)
      )
    }


    # if there are ghost topics, we need to add them to the topic_order, the gammas, the betas and the alignment.

    # 1. adding the ghost topics to the topic order.
    # here, we simply need to create a data.frame that, for each model, has the ghost topics.
    # the k, branch and topic for each ghost topic is equal to the ghost topic number.


    topic_order_g =
      tibble(branch = ghost_topics,
             k = ghost_topics,
             topic = ghost_topics %>% as.character(),
             is_ghost = TRUE) %>%
      expand_grid(
        m = M
      )

    topic_order =
      bind_rows(topic_order %>% mutate(is_ghost = FALSE),
                topic_order_g) %>%
      arrange(m, branch, k)

    # 2. adding the ghost topics to the gammas
    # we need to do three things:
    # a. define ghost documents for the ghost topics
    # b. define the weight of the ghost topics and
    # c. re-scale the existing weights so that the sum of weights is 1 for each model.

    w_g = 1/K_max # the weight of a ghost topic
    m_g = w_g * length(unique(gammas$d)) # the document mass of each  ghost topic

    s = (K_max - K_g) / K_max # the scaling factor for the height and flows of the non-ghost topics

    gammas = gammas %>% mutate(g = s*g, m_ref = last(M))

    gammas_g =
      topic_order %>%
      filter(is_ghost) %>%
      mutate(
        d = "ghost_document",
        g = m_g,
        m_ref = last(M)
      )

    gammas =
      bind_rows(
        gammas %>%
          select(m, k_LDA, branch, topic, k, d, g, m_ref) %>%
          mutate(is_ghost = FALSE),
        gammas_g
      ) %>%
      arrange(m, branch, k, d)


    # 3. adding the ghost topics to the betas
    # we need to do two things:
    # a. define ghost words for the ghost topics
    # b. define the weight of the ghost words and
    N_w = length(unique(betas$w)) # First, we count the current number of words
    n_w_per_topic = (N_w/length(actual_ks)) %>% round() # then we compute the average number of words per topic
    N_g = K_g * n_w_per_topic # N_g is the number of ghost words (= number of ghost topic * average number of topic)
    mean_beta_per_topic = betas %>% group_by(m) %>%
      summarize(sum_b = sum(b), K = n_distinct(k_LDA), .groups = "drop") %>%
      mutate(mean_b = sum_b / K) %>%
      ungroup() %>% select(mean_b) %>% unlist() %>% mean() # that's because we have trimmed the model

    betas_g =
      topic_order %>%
      filter(is_ghost) %>%
      expand_grid(
        w_index = 1:n_w_per_topic
      ) %>%
      mutate(
        w = str_c("ghost_word_",k,"_",w_index),
        b = mean_beta_per_topic/n_w_per_topic,
      ) %>%
      select(m, branch, k, topic, w, b, is_ghost)

    betas =
      bind_rows(
        betas %>% select(m, branch, topic, k, w, b) %>% mutate(is_ghost = FALSE),
        betas_g
      ) %>%
      arrange(m, branch, k, w)


    # 4. adding the ghost topics to the topic alignment
    # for this, we consider that there is a 1-1 perfect alignment between a given ghost topic across models

    alignment_g =
      tibble(
        m = M[-length(M)],
        m_next = M[-1]) %>%
      left_join(
        .,
        tibble(
          branch = ghost_topics,
          branch_next = ghost_topics,
          k = ghost_topics,
          k_next = ghost_topics,
          topic = ghost_topics %>% as.character(),
          topic_next = ghost_topics %>% as.character(),
          document_mass = m_g,
          weight = w_g,
          norm_weight = 1,
          is_ghost = TRUE
        ),
        by = character()
      )

    # and we need to re-scale the weights from the original alignment

    alignment =
      alignment %>% mutate(weight = s * weight)

    alignment =
      bind_rows(alignment %>% mutate(is_ghost = FALSE),
                alignment_g) %>%
      arrange(m, branch, k, branch_next, k_next)



    # DONE
    # return the updated gammas and betas
    lda_models = list(gammas = gammas)
    if("betas" %in% names(aligned_topics$lda_models)) lda_models$betas = betas
    list(lda_models = lda_models,
         alignment = alignment,
         topic_order = topic_order)
  }

.lengthen_betas_ <- function(beta) {
  as_tibble(beta) %>%
    mutate(., w = as.factor(row_number())) %>%
    pivot_longer(., -w, names_to = "k_LDA", values_to = "b") %>%
    mutate(
      k_LDA = str_replace(k_LDA, "V", ""),
      K = ncol(gamma)
    )
}

.lengthen_gammas_ <- function(gamma) {
  as_tibble(gamma) %>%
    mutate(., d = as.factor(row_number())) %>%
    pivot_longer(., -d, names_to = "k_LDA", values_to = "g") %>%
    mutate(
      k_LDA = str_replace(k_LDA, "V", ""),
      K = ncol(gamma)
    )
}

.lengthen_betas <- function(betas, m_levels) {
  map_dfr(betas, ~ .lengthen_betas_(.), .id = "m") %>%
    mutate(m = factor(m, levels = m_levels))
}

.lengthen_gammas <- function(gammas, m_levels) {
  map_dfr(gammas, ~ .lengthen_gammas_(.), .id = "m") %>%
  mutate(m = factor(m, levels = m_levels))
}

.compute_topic_index = function(aligned_topics) {
  betas = aligned_topics$lda_models$betas
  gammas = aligned_topics$lda_models$gammas
  alignment = aligned_topics$alignment
  topic_order = aligned_topics$topic_order
  all_branches = topic_order$branch %>% unique() %>% sort()

  topic_order =
    topic_order %>%
    mutate(branch_upper_part = (branch < max(all_branches)/2),
           k_order_for_index = ifelse(branch_upper_part, k, -k)) %>%
    arrange(m, branch, k_order_for_index) %>%
    group_by(m) %>%
    mutate(k_i = row_number()) %>%
    select(-branch_upper_part, -k_order_for_index)

  gammas = gammas %>% left_join(topic_order %>% select(m, topic, k_i), by = c("m","topic"))
  betas = betas %>% left_join(topic_order %>% select(m, topic, k_i), by = c("m","topic"))

  alignment =
    alignment %>%
    left_join(
      topic_order %>% select(m, topic, k_i),
      by = c("m","topic")) %>%
    left_join(
      topic_order %>% select(m, topic, k_i) %>%
        dplyr::rename(m_next = m, topic_next = topic, k_i_next = k_i),
      by = c("m_next","topic_next"))

  # DONE
  # return the updated gammas and betas
  lda_models = list(gammas = gammas, betas = betas)
  list(lda_models = lda_models,
       alignment = alignment,
       topics_order = topic_order)
}



.compute_layout = function(aligned_topics, min_k, max_k) {

  gammas = aligned_topics$lda_models$gammas
  alignment = aligned_topics$alignment %>%
    select(-contains("_LDA"))

  m_ref = unique(gammas$m_ref)
  M = unique(gammas$m)

  # we compute the delta_k for each model based on the number of topics (actual + ghosts if any)
  h_m =
    gammas %>%
    group_by(m) %>%
    summarize(K = length(unique(k)), .groups = "drop") %>%
    mutate(delta_k = 1/(K+1))
  all_ks = min_k:max_k

  # RECTANGLES
  layout_rect =
    gammas %>%
    select(m, branch, k_i, g, is_ghost) %>%
    group_by(m, branch, k_i, is_ghost) %>%
    summarize(g = sum(g), .groups = "drop") %>%
    group_by(m) %>%
    mutate(w = g/sum(g)) %>%
    left_join(., h_m %>% select(m, delta_k), by = "m") %>%
    arrange(m, k_i) %>%
    group_by(m) %>%
    mutate(
      height = w,
      cum_w = cumsum(w),
      y = - k_i*delta_k - cum_w + w/2,
      m_num = match(m, M),
      x = m_num,
      color = branch %>%  factor(., levels = all_ks)
    )

  # FLOWS
  OUT =
    alignment %>%
    mutate(m_num = match(m, M),
           x = m_num + 0.2) %>%
    left_join(., h_m %>% select(m, delta_k), by = "m") %>%
    group_by(m) %>%
    arrange(m, k_i, k_i_next) %>%
    mutate(
      cum_w = cumsum(weight),
      ymin = - k_i*delta_k - cum_w,
      ymax = ymin + weight,
      color = branch_next %>%  factor(., levels = all_ks))

  IN =
    alignment %>%
    mutate(m_num = match(m, M),
           x = m_num + 1 - 0.2)  %>%
    left_join(., h_m %>% select(m, delta_k) %>% dplyr::rename(m_next = m), by = "m_next") %>%
    arrange(m, m_next, k_i_next, k_i) %>%
    group_by(m) %>%
    mutate(cum_w = cumsum(weight),
           ymin = -k_i_next*delta_k - cum_w,
           ymax = ymin + weight,
           color = branch_next %>% factor(., levels = all_ks))


  layout_ribbons =
    bind_rows(OUT, IN) %>%
    arrange(m, k_i, k_i_next)

  layout_ribbons =
    layout_ribbons %>%
    mutate(h = ymax - ymin,
           flow_id = str_c("m_",m,"|topic_",topic,"|topic_next_",topic_next))

  # removing the ghost topics
  layout_rect = layout_rect %>% filter(!is_ghost)
  layout_ribbons = layout_ribbons %>% filter(!is_ghost)

  layouts = list(rect = layout_rect, ribbons = layout_ribbons)

  layouts
}




.compute_k_ref_table = function(alignment){

  M = levels(alignment$m)

  k_ref_table =
    alignment %>%
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
      alignment %>%
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
  aligned_topics,
  min_k,
  max_k
){

  m_ref = last(aligned_topics$lda_models$gammas$m)
  adjusted_betas_ref =
    aligned_topics$lda_models$betas %>%
    filter(m == m_ref) %>%
    group_by(m, k) %>%
    mutate(b = b/sum(b))

  ref_topic_for_each_word =
    adjusted_betas_ref %>%
    arrange(w, -b) %>%
    group_by(w) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    arrange(k_i, -b) %>%
    mutate(w_index = row_number())

  ref_topic_weights =
    aligned_topics$lda_models$gammas %>%
    filter(m == m_ref) %>%
    group_by(k_i) %>%
    summarize(weight = sum(g), .groups = "drop") %>%
    mutate(weight = weight/sum(weight))

  est_average_freq =
    adjusted_betas_ref %>%
    left_join(., ref_topic_weights, by = "k_i") %>%
    mutate(f = b * weight) %>%
    group_by(w) %>%
    summarize(f = sum(f), .groups = "drop") %>%
    arrange(-f)

  # rectangles
  total_height = est_average_freq$f %>% sum()
  delta = (2-total_height)/(max_k - min_k + 1 + nrow(est_average_freq)) # the space between words and the additional space between topics

  leaves_layout_rect =
    ref_topic_for_each_word %>%
    left_join(., est_average_freq, by = "w") %>%
    mutate(m = "words",
           m_num = length(levels(aligned_topics$alignment$m))+1 + 2,
           x = m_num,
           cum_f = cumsum(f),
           y = - (k_i + row_number() - 1)*delta - cum_f + f/2,
           height = f ,
           color = branch %>% factor(., levels = min_k:max_k))

  # Flows
  delta_out = 1/(max_k - min_k + 2)
  leaves_OUT =
    ref_topic_weights %>%
    left_join(., adjusted_betas_ref , by = "k_i") %>%
    mutate(w = w %>% factor(., levels = ref_topic_for_each_word$w)) %>%
    arrange(k_i, w) %>%
    mutate(h = weight * b) %>%
    select(m, k_i, branch, w, h) %>%
    mutate(x = m_ref %>% as.numeric() + 0.2,
           cum_h = cumsum(h),
           ymin = - delta_out * k_i - cum_h,
           ymax = ymin + h,
           color = branch %>% factor(., levels = min_k:max_k))


  fraction_of_word_from_each_topic =
    adjusted_betas_ref %>%
    ungroup() %>%
    group_by(w) %>%
    mutate(c = b/sum(b)) %>%
    arrange(w, k_i)

  leaves_IN =
    est_average_freq %>%
    left_join(., fraction_of_word_from_each_topic, by = "w") %>%
    mutate(h = f*c,
           h = h/sum(h)) %>%
    left_join(.,
              ref_topic_for_each_word %>%
                dplyr::rename(k_i_w = k_i) %>%
                select(w, k_i_w, w_index),
              by = "w") %>%
    arrange(w_index, k_i) %>%
    mutate(x = m_ref %>% as.numeric() + 0.8 + 2,
           cum_h = cumsum(h),
           ymin = - delta * (k_i_w + w_index - 1) - cum_h,
           ymax = ymin + h,
           color = branch %>% factor(., levels = min_k:max_k))

  leaves_layout_ribbon =
    bind_rows(
      leaves_OUT %>% select(m, w, k_i,  x, ymin, ymax, color),
      leaves_IN %>% select(m, w, k_i, x, ymin, ymax, color)
    ) %>%
    arrange(m, k_i, w, x) %>%
    left_join(., ref_topic_for_each_word %>% select(w, w_index), by = "w") %>%
    mutate(flow_id = str_c("leaves|k_i_",k_i,"|w_",w_index)) %>%
    select(-w_index)

  # words

  leaves_layout_words =
    leaves_layout_rect %>%
    filter(!is_ghost) %>%
    select(m_num, x, y, k_i, w, color)

  # remove the ghost words if any
  leaves_layout_rect =
    leaves_layout_rect %>%
    filter(!str_detect(w, "ghost_word_"))

  leaves_layout_ribbon =
    leaves_layout_ribbon %>%
    filter(!str_detect(w, "ghost_word_"))

  # return the layouts

  list(rect = leaves_layout_rect,
       ribbons = leaves_layout_ribbon,
       text = leaves_layout_words)
}




.compute_layout_reference = function(alignment, h_m, M, m_ref){

  layout_rect =
    alignment %>%
    select(m, k, k_ref, k_next, w) %>%
    arrange(m, k, k_ref, k_next) %>%
    # we add the last model
    bind_rows(
      .,
      alignment %>%
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
           m = m %>% factor(., levels = c(levels(alignment$m),"document")),
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
