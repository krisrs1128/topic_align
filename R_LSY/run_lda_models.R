
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
      }
    )
    
    # gammas
    gammas = 
      map_dfr(
        .x = Ks[Ks>1],
        .f = function(K){
          topic_model_file = str_c(dir,"_LDA_K_", K, ".Rdata")
          load(file = topic_model_file)
          
          documents_topic_composition = 
            topic_model@gamma %>% 
            as.data.frame() %>% 
            set_colnames(1:K) %>% 
            mutate(d = topic_model@documents) %>% 
            pivot_longer(cols = -d,
                         names_to = "k_LDA",
                         values_to = "g") %>% 
            mutate(K = K,
                   k_LDA = letters[k_LDA %>% as.integer()])
          documents_topic_composition 
        }
      )
    if(1 %in% Ks){
      gammas = 
        gammas %>% 
        bind_rows(
          tibble(
            K = 1, k_LDA = "a", 
            d = unique(gammas$d), 
            g = 1),
          .
        )
    }
    gammas = gammas %>% mutate(m = K %>% factor(., levels = Ks))
    gammas = gammas %>% select(m, K, k_LDA, d, g)
    
    
    # betas
    betas = 
      map_dfr(
        .x = Ks[Ks>1],
        .f = function(K){
          
          topic_model_file = str_c(dir,"_LDA_K_", K, ".Rdata")
          load(file = topic_model_file)
          
          topic_composition = 
            topic_model@beta %>% 
            exp() %>% 
            as.data.frame() %>% 
            set_colnames(topic_model@terms) %>% 
            mutate(k_LDA = letters[1:K]) %>% 
            pivot_longer(cols = -k_LDA,
                         names_to = "w",
                         values_to = "b") %>% 
            mutate(K = K)
          topic_composition 
        }
      )
    if(1 %in% Ks){
      data_mat = as.matrix(data)
      mean_prop = 
        (data_mat/rowSums(data_mat)) %>% 
        apply(., 2, mean) %>% 
        as.vector() %>%
        tibble(K = 1, k_LDA = "a", 
               w = colnames(data_mat),
               b = .)
      betas = bind_rows(mean_prop,betas)
    }
    betas = betas %>% mutate(m = K %>% factor(., levels = Ks))
    betas = betas %>% select(m, K, k_LDA, w, b)
    
    
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

