
plot_gamma <-  function(alignment, models) {
  gamma <- 
    map_dfr(
      .x = models,
      .f = function(m){
        this_gamma <- 
          alignment@models[[m]]$gamma
        this_gamma <- 
          this_gamma %>%  
          set_colnames(1:ncol(this_gamma))  %>% 
          as_tibble() %>% 
          mutate(
            m = names(alignment@models)[m],
            d = rownames(this_gamma),
          ) %>% 
          pivot_longer(., cols = -c(m, d),
                       names_to = "k",
                       values_to = "topic_prop") %>% 
          mutate(k = k %>% as.integer())
        this_gamma
      }
    )
  
  gamma <- 
    gamma %>% 
    mutate(m = m %>% factor(., levels = names(alignment@models))) %>% 
    left_join(., alignment@topics, by = c("m", "k")) 
  
  d_order <- 
    gamma %>% 
    filter(m == names(alignment@models)[last(models)]) %>% 
    group_by(d) %>% 
    summarize(
      main_topic = k[which.max(topic_prop)],
      avg_path = weighted.mean(as.integer(path), topic_prop), 
      .groups = "drop"
    ) %>% 
    arrange(main_topic, avg_path) 
  
  
  gamma <- 
    gamma %>% 
    mutate(d = d %>% factor(., levels = d_order$d))
  
  ggplot(gamma, aes(x = k, y = d %>% fct_rev(), fill = path, alpha = topic_prop)) +
    geom_tile() +
    facet_grid(. ~ m, scales = "free", space = "free") +
    scale_y_discrete(breaks = NULL) +
    scale_alpha(range = c(0,1), limits = c(0,1)) +
    guides(fill = "none", alpha = "none") + 
    scale_x_continuous(breaks = unique(gamma$k) %>% sort(), minor_breaks = NULL) +
    scale_fill_discrete(breaks = gamma$path %>% levels(), limits = gamma$path %>% levels()) +
    theme_bw() + 
    theme(panel.spacing.x = unit(0, "pt"), 
          strip.text.y = element_text(angle = 0, hjust = 0, color = "black")
    ) +
    ylab("") +
    xlab("")
  
  
}
