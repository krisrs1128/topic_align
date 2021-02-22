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
