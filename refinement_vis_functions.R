add_refinement_scores_to_alignment <- function(fits, alignment, daughter_threshold = .99) {
    K_max = max(alignment$lda_models$betas$k)
    clusters = fits$gammas %>% select(c(K, k_LDA)) %>% unique %>% subset(K < (K_max - 2))
    scores = list()
    for(i in 1:nrow(clusters)) {
        K = clusters$K[i]
        k_LDA = clusters$k_LDA[i]
        rs_df = refinement_scores(fits$gammas, K, k_LDA, K_max = K_max, daughter_threshold = daughter_threshold)
        if(!is.na(rs_df)) {
            rs_df$k = K; rs_df$k_next = K+1; rs_df$k_LDA = k_LDA
        }
        scores[[i]] = rs_df
    }
    score_df = Reduce(rbind, scores) %>% subset(!is.na(score)) %>% mutate(k = as.factor(k), k_next = as.factor(k_next)) %>% tibble
    augmented = merge(alignment$gamma_alignment, score_df,
                      by.x = c("m", "m_next", "k_LDA", "k_LDA_next"),
                      by.y = c("k", "k_next", "k_LDA", "k_LDA_next"),
                      all.x = TRUE)
    alignment$gamma_alignment = augmented
    alignment$beta_alignment$score = NA
    return(alignment)
}

refinement_scores <- function(gammas, K, k_LDA, K_max, daughter_threshold = .99) {
    R12 = make_refinement_matrix(gammas, K, K+1)
    daughters = which(R12[, k_LDA] >= daughter_threshold)
    if(length(daughters) == 0)
        return(NA)
    daughter_names = rownames(R12)[daughters]
    scores = matrix(0, nrow = K_max - K - 1, ncol = length(daughter_names))
    colnames(scores) = daughter_names
    for(j in (K+2):K_max) {
        R1j = make_refinement_matrix(gammas, K, j)
        R2j = make_refinement_matrix(gammas, K+1, j)
        descendants = which(R1j[, k_LDA] >= daughter_threshold)
        if(length(descendants) > 0) {
            descendant_names = rownames(R1j)[descendants]
            daughter_descendant_scores = R2j[descendant_names,daughter_names,drop = FALSE]
            scores[j - K - 1,] = min(apply(daughter_descendant_scores, 1, max))
            ## If there is a cluster that never has a descendant in
            ## the set we're interested in, it needs to be treated
            ## differently
            cluster_with_no_refinement = Reduce(intersect, (apply(daughter_descendant_scores, 1, function(x) setdiff(1:length(x), which.max(x)))))
            if(length(cluster_with_no_refinement) > 0) {
                scores[j - K - 1, cluster_with_no_refinement] = apply(daughter_descendant_scores[,cluster_with_no_refinement, drop = FALSE], 2, mean)
            }
        } else {
            scores[j - K - 1,] = rep(NA, ncol(scores))
        }
    }
    return(data.frame(score = colMeans(scores, na.rm = TRUE), k_LDA_next = colnames(scores)))
}

make_refinement_matrix <- function(gammas, m_current, m_next, version = "diagonal_scaling") {
    C1 = gammas %>%
        subset(m == m_current) %>%
        rename(document = d) %>%
        pivot_wider(id_cols = document, names_from = k_LDA, values_from = g)

    C2 = gammas %>%
        subset(m == m_next) %>%
        rename(document = d) %>%
        pivot_wider(id_cols = document, names_from = k_LDA, values_from = g)

    C1_m = C1 %>% select(-document) %>% as.matrix
    C2_m = C2 %>% select(-document) %>% as.matrix

    ## refinement matrix estimated by least squares
    if(version == "leastsquares") {
        R = solve(t(C2_m) %*% C2_m) %*% t(C2_m) %*% C1_m
    } else if (version == "laura"){
        R = t(C2_m) %*% C1_m / nrow(C1_m)
    } else if (version == "diagonal_scaling") {
        R = t(C2_m) %*% C1_m
        R = sweep(R, 1, rowSums(R), "/")
    }
    return(R)
}
