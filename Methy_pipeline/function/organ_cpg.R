process_oragn_cpg_data <- function(betas, organ_cpg_file, organ_mean_file) {
  organ_cpg = read.csv(organ_cpg_file)
  organ_mean = read.csv(organ_mean_file)
  
  combine_re <- NULL
  
  for (organ in colnames(organ_cpg)) {
    
    organ_cpg_mean_ref = matrix(organ_mean[, organ], nrow = ncol(betas), ncol = 5, byrow = TRUE)
    rownames(organ_cpg_mean_ref) = colnames(betas)
    colnames(organ_cpg_mean_ref) = paste(organ, '_', organ_cpg[, organ], '_mean', sep = '')
    
    organ_cpg_values = t(betas[organ_cpg[, organ], ])
    colnames(organ_cpg_values) = paste(organ, '_', organ_cpg[, organ], '_value', sep = '')
    
    if (is.null(combine_re)) {
      combine_re = cbind(organ_cpg_mean_ref, organ_cpg_values)
    } else {
      combine_re = cbind(combine_re, cbind(organ_cpg_mean_ref, organ_cpg_values))
    }
  }
  
  return(combine_re)
}




