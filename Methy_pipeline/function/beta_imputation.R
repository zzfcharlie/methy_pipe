beta_imputation <- function(beta_matrix, ref_means) {
  na_rows <- rownames(beta_matrix)[rowSums(is.na(beta_matrix)) > 0]
  beta_matrix[na_rows, ] <- ref_means[na_rows, ]
  return(beta_matrix)
}