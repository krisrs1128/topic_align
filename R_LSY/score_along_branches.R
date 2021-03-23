


compute_score_along_branches = function(aligned_topics){
  
  alignment = aligned_topics$gamma_alignment
  alignment_along_branches = 
    alignment %>% 
    filter(branch == branch_next) %>% 
    group_by(m, m_next, branch, branch_next, topic_next) %>% 
    summarize(norm_weight = sum(norm_weight), .groups = "drop") %>% 
    group_by(m, m_next, branch, branch_next) %>% 
    summarize(norm_weight = mean(norm_weight), .groups = "drop")
  
  scores = 
    aligned_topics$topics_order %>% 
    select(m, branch) %>% 
    arrange(branch, m) %>% 
    group_by(branch) %>% 
    slice_head(n = 1) %>% 
    ungroup() %>% 
    mutate(score = 1) %>% 
    arrange(m)
  
  for(this_m in unique(alignment$m) %>% sort()){
    
    this_m_scores = 
      alignment_along_branches %>% 
      filter(m == this_m) %>% 
      left_join(., scores %>%  filter(m == this_m), by = c("m", "branch")) %>% 
      mutate(score = score * norm_weight) %>% 
      select(m_next, branch, score) %>% 
      dplyr::rename(m = m_next)
    
    scores = bind_rows(scores, this_m_scores) %>% arrange(m, branch)
    
    
  }
  
  scores
}