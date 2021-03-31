make_refinements_from_lda_models = function(lda_models, k1, k2, version = "diagonal_scaling") {
    C1 = lda_models$gammas %>%
        subset(K == k1) %>%
        pivot_wider(id_cols = d, names_from = k_LDA, names_prefix = "topic_", values_from = g)

    C2 = lda_models$gammas %>%
        subset(K == k2) %>%
        pivot_wider(id_cols = d, names_from = k_LDA, names_prefix = "topic_", values_from = g)

    C1_m = as.matrix(C1[,-1])
    C2_m = as.matrix(C2[,-1])

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
get_split_score_lda = function(lda_models, k, cluster_id, k_max, daughter_threshold = .99) {
    R12 = make_refinements_from_lda_models(lda_models, k, k+1)
    daughters = which(R12[,cluster_id] >= daughter_threshold)
    daughter_names = rownames(R12)[daughters]
    scores = numeric(k_max - k - 1)
    for(j in (k+2):k_max) {
        R1j = make_refinements_from_lda_models(lda_models, k, j)
        R2j = make_refinements_from_lda_models(lda_models, k+1, j)
        descendants = which(R1j[,cluster_id] >= 0)
        descendant_names = rownames(R1j)[descendants]
        daughter_descendant_scores = R2j[descendant_names,daughter_names,drop = FALSE]
        ## OR logic = min, we want everything to be a refinement
        ## abs(.5 - x) is a function that is maximal when x is 0 or 1 (indicating a refinement) and minimal when x = .5 (indicating not a refinement)
        scores[j - k - 1] = min(2 * abs(.5 - daughter_descendant_scores))
    }
    return(mean(scores))
}
