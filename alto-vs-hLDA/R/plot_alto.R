
plot_alto <- function(alignment) {
  ggarrange(
    plot(alignment),
    ggarrange(
      plot_beta(alignment, models = c(3,6,11), threshold = "0.005"), 
      plot_gamma(alignment, models = c(3,6,11)),
      nrow = 2, ncol = 1,
      align = "v"
    ),
    nrow = 2, ncol = 1,
    heights = c(1,4)
  )
}
