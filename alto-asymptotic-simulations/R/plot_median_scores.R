

plot_median_scores <- function(topics_score_summary, score = "refinement") {
  
  if (score == "refinement") {
    topics_score_summary$y = topics_score_summary$median_refinement
    topics_score_summary <-  topics_score_summary %>% 
      filter(n_topics < max(n_topics))
  } else if (score == "coherence") {
    topics_score_summary$y = topics_score_summary$median_coherence
  } else {
    topics_score_summary$y = topics_score_summary$n_paths
  }
  
  ggplot(topics_score_summary, 
         aes(x = n_topics %>% factor(), y = y,  group = id, color = log2(N))) +
    geom_line(alpha = 0.2) +
    facet_grid(N ~ method, labeller = label_both) +
    ylab(score) +
    xlab("Nb of topics") +
    scale_color_gradient(low = "gray50", high = "steelblue2") +
    guides(color = "none")
  
}