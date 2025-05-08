calculate_age_score <- function(CA, AAF) {
  age_groups <- list(
    young = list(min = 18, max = 34),
    middle = list(min = 34, max = 60),
    old = list(min = 60, max = Inf)
  )
  
  quantile_data <- data.frame(
    Quantile = c(0.25, 0.5, 0.75),
    young = c(0.962, 0.991, 1.051),
    middle = c(0.975, 1.007, 1.043),
    old = c(0.967, 1.002, 1.044)
  )
  
  min_values <- c(young = 0.83, middle = 0.84, old = 0.88)
  max_values <- c(young = 1.27, middle = 1.18, old = 1.16)
  
  scores <- numeric(length(CA))
  
  for (i in seq_along(CA)) {
    if (CA[i] >= age_groups$young$min & CA[i] < age_groups$young$max) {
      group <- "young"
    } else if (CA[i] >= age_groups$middle$min & CA[i] < age_groups$middle$max) {
      group <- "middle"
    } else if (CA[i] >= age_groups$old$min) {
      group <- "old"
    }
    
    r_min <- min_values[group]
    r_max <- max_values[group]
    Q25 <- quantile_data[1, group]
    Q75 <- quantile_data[3, group]
    
    if (AAF[i] < r_min) {
      scores[i] <- 97
    } else if (AAF[i] <= Q25) {
      scores[i] <- 98 - (98 - 91) / (Q25 - r_min) * (AAF[i] - r_min)
    } else if (AAF[i] <= Q75) {
      scores[i] <- 90 - (90 - 76) / (Q75 - Q25) * (AAF[i] - Q25)
    } else if (AAF[i] <= r_max) {
      scores[i] <- 75 - (75 - 62) / (r_max - Q75) * (AAF[i] - Q75)
    } else {
      scores[i] <- 63
    }
  }
  scores = round(scores)
  return(scores)
}