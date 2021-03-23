




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
