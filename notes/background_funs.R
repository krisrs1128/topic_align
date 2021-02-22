library("tidyr")
library("stringr")
library("dplyr")
library("MASS")

illustration_data <- function(n, m) {
  a_loc <- mvrnorm(n, c(0, 0), diag(2))
  b_loc <- mvrnorm(m, c(0, 0), diag(2))
  a <- cbind(a_loc, mass = rgamma(n, 1)) %>%
    as.data.frame() %>%
    mutate(id = as.character(row_number()))
  b <- cbind(b_loc, mass = rgamma(m, 1)) %>%
    as.data.frame() %>%
    mutate(id = as.character(row_number()))

  masses <- bind_rows(a, b, .id = "set") %>%
    group_by(set) %>%
    mutate(mass = mass / sum(mass))
}

sinkhorn <- function(a, b, C, epsilon = 1e-2, tol = 1e-5) {
  K <- exp(-C / epsilon)
  v <- rep(1, ncol(C))
  obj <- 1e10
  while (TRUE) {
    u <- a / (as.numeric(K %*% v))
    v <- b / (as.numeric(t(K) %*% u))

    # evaluate objective
    P <- diag(u) %*% K %*% diag(v)
    obj_new <- sum(diag(t(C) %*% P))
    if (obj - obj_new < tol) break
    obj <- obj_new
  }
  list(u = u, v = v, P = P)
}

merge_link_data <- function(P, masses) {
  links <- data.frame(P) %>%
    mutate(source = str_c("1-", row_number())) %>%
    pivot_longer(cols = starts_with("X"), names_to = "target") %>%
    mutate(target = str_replace(target, "X", "2-"))

  masses_ <- masses %>%
    unite("set_id", set, id, sep = "-") %>%
    select(set_id, V1, V2)

  links %>%
    left_join(masses_, c("source" = "set_id")) %>%
    left_join(masses_, c("target" = "set_id"))
}

softmin <- function(x, epsilon) {
  mx <- min(x)
  mx + epsilon * log(sum(exp(- 1 / epsilon * (x - mx))))
}

rowmin <- function(X, epsilon) {
  apply(X, 1, softmin, epsilon)
}

colmin <- function(X, epsilon) {
  apply(X, 2, softmin, epsilon)
}

log_sinkhorn <- function(a, b, C, s_eps = 1e-2, epsilon = 1e-2, tol = 1e-5, n_iter=10) {
  f <- rnorm(nrow(C))
  g <- rnorm(ncol(C))
  ones <- list(
    m = matrix(1, nrow = 1, ncol = ncol(C)),
    n = matrix(1, nrow = nrow(C), ncol = 1)
  )

  for (iter in seq_len(n_iter)) {
    f <- rowmin(C - f %*% ones$m - ones$n %*% g, s_eps) + f + epsilon * log(a)
    print(f)
    g <- colmin(C - f %*% ones$m - ones$n %*% g, s_eps) + g + epsilon * log(b)
  }

  u <- exp(f / epsilon)
  v <- exp(g / epsilon)
  P <- diag(u) %*% exp(-C / epsilon) %*% diag(v)
  list(u = u, v = v, f = f, g = g, P = P)
}
