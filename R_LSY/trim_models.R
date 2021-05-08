
trim_models =
  function(models,
           min_prop = NULL,
           n_words = NULL){

    trimmed_models = models
    if(is.null(min_prop) & is.null(n_words)) min_prop = 0.025
    if(!is.null(min_prop)){
      trimmed_models$betas = trimmed_models$betas %>%
        bind_rows() %>%
        group_by(w) %>%
        mutate(max_prop = max(b)) %>%
        filter(max_prop >= min_prop) %>%
        select(-max_prop) %>%
        ungroup()
    }else{
      top_n_words =
        trimmed_models$betas %>%
        group_by(w) %>%
        summarize(max_prop = max(b),.groups = "drop") %>%
        arrange(-max_prop) %>%
        slice_head(n = n_words)
      trimmed_models$betas =
        trimmed_models$betas %>%
        filter(w %in% top_n_words$w)
    }
    trimmed_models
  }
