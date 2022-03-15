
plot_scores <- function(topics_score_summary, colors, labels = "auto") {
  
  df_long <- 
    topics_score_summary %>% 
    group_by(method, id, N, V, K) %>% 
    # mutate(score_n_paths = n_paths / max(n_topics)) %>% 
    pivot_longer(
      cols = c(n_paths, coherence_summary, refinement_summary),
      names_to = "score_name",
      values_to = "score_value"
    ) %>% 
    filter(!(n_topics == max(n_topics))) %>% 
    mutate(
      score_name = 
        case_when(
          score_name == "refinement_summary" ~ "minimum refinement score",
          score_name == "coherence_summary" ~ "minimum coherence score",
          score_name == "n_paths" ~ ("# paths")
        ) %>% 
        factor(
          .,
          levels = c("minimum coherence score", "minimum refinement score", "# paths")
        )
    )
  
  
  plotlist <- 
    map(
      .x = levels(df_long$score_name),
      .f = function(sn) {
        
        if (str_detect(sn, "coherence")) {
          ybreaks <- seq(0,1, by = 0.25)
        } else {
          ybreaks <- seq(1, 10, by = 2)
        }
        
        ggplot(df_long %>% filter(score_name == sn),
               aes(x = n_topics , y = score_value, 
                   color = method, group = interaction(id,method))) +
          geom_line(alpha = 0.2) +
          facet_grid(
            N ~ method ,  
            labeller = labeller(N = label_both, method = label_value)
          ) +
          guides(color = "none") +
          ylab(sn) +
          xlab("Nb of topics") +
          scale_x_continuous(breaks = 0:10, minor_breaks = NULL) +
          scale_y_continuous(breaks = ybreaks) +
          scale_color_manual(values = colors)
        
      }
    )
  
  ggarrange(plotlist = plotlist, ncol = 3, nrow = 1, labels = labels)
  
}