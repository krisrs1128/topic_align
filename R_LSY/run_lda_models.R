
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

    if(is.null(seed)) seed = sample(1:1000, 1)

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
        get(load(topic_model_file))
      }
    )

    names(ok) <- Ks[Ks > 1]
    gammas <- map(ok, ~ .@gamma)
    betas <- map(ok, ~ exp(t(.@beta)))
    if (1 %in% Ks) {
      gammas = c(list("1" = matrix(1, nrow = nrow(data), ncol = 1)), gammas)
      betas = c(list("1" = matrix(colSums(data) / sum(colSums(data)), ncol = 1)), betas)
    }

    # log-Likelihood of each sample
    ll =
      map_dfr(
        .x = Ks[Ks>1],
        .f = function(K){

          topic_model_file = str_c(dir,"_LDA_K_", K, ".Rdata")
          load(file = topic_model_file)

          loglik =
            tibble(
              K = K,
              document = topic_model@documents,
              log_likelihood = topic_model@loglikelihood )
          loglik
        }
      )
    ll = ll %>% mutate(m = K %>% factor(., levels = Ks))

    # delete dir (if dir was unspecified)
    if(delete_dir) unlink(dir)

    # return results
    res = list(
      betas = betas,
      gammas = gammas,
      log_likelihood = ll
    )
    res
  }
