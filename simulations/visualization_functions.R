read_and_process <- function(pattern, grouping, directory, max_K=10) {
  paths <- list.files(directory, pattern, full = TRUE)
  map_dfr(paths, read_csv, .id = "file") %>%
    suppressMessages() %>%
    mutate(
      m = factor(m, levels = str_c("K", 1:max_K)),
      file = basename(paths[as.integer(file)]),
      method = ifelse(grepl("product", file), "product", "transport")
    ) %>%
    unite("group", grouping, remove = FALSE)
}

named_list <- function(v) {
  setNames(vector(length = length(v), mode = "list"), basename(v))
}

read_betas <- function(paths, ix = c(4, 2), beta_id=NULL) {
  betas <- named_list(paths)
  beta_hats <- named_list(paths)

  for (p in paths) {
    exper <- get(load(p))
    beta_hats[[basename(p)]] <- map(models(exper[[ix[1]]]), ~ exp(.$beta))

    if (is.null(beta_id)) {
      betas[[basename(p)]] <- exper[[ix[2]]]
    } else {
      betas[[basename(p)]] <- exper[[ix[2]]][[beta_id]]
    }
  }

  list(beta_hats = beta_hats, betas = betas)
}

topic_cossim <- function(betas, beta_hats) {
  similarity <- list()
  m_names <- names(beta_hats[[1]])

  l <- 1
  for (i in seq_along(betas)) {
    for (j in seq_along(m_names)) {
      max_sim <- cosine_similarity(betas[[i]], beta_hats[[i]][[j]]) %>%
        apply(2, max)
      similarity[[l]] <- tibble(
        file = names(betas)[i],
        m = m_names[j],
        k = seq_len(nrow(beta_hats[[i]][[j]])),
        sim = max_sim
      )
      l <- l + 1
    }
  }

  bind_rows(similarity)
}

GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin,
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
  data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
  grp <- data[1, "group"]
  newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
  newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
  newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])

  if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
    stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
      1))
    quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
    aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
    aesthetics$alpha <- rep(1, nrow(quantiles))
    both <- cbind(quantiles, aesthetics)
    quantile_grob <- GeomPath$draw_panel(both, ...)
    ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
  }
  else {
    ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
  }
})

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ...,
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE,
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin,
        position = position, show.legend = show.legend, inherit.aes = inherit.aes,
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}
