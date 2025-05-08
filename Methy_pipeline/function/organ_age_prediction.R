calculate_mrs <- function(beta_matrix, sumstats) {
  filtered_beta_matrix = beta_matrix[sumstats$CpG,]
  MRS <- as.data.frame(t(as.matrix(filtered_beta_matrix)) %*% (sumstats$Beta))
  rownames(MRS) <- colnames(filtered_beta_matrix)
  colnames(MRS) <- c("MRS")
  return(MRS)
}


organ_age_prediction <- function(sample_info, betas, organ_model_base_path, organ_ewas_csv_path) {
  org_list <- c("brain", "bowel", "heart", "immune", "kidney", "liver", "lung", "metabolic", "musculoskeletal", "vascular")
  
  results_df <- data.frame(matrix(nrow = nrow(sample_info), ncol = 0)) 

  for (org in org_list) {
    
    model_path <- file.path(organ_model_base_path, org, "pca_and_models.RData")
    if (!file.exists(model_path)) {
        print(paste("Model for", org, "does not exist."))
        next
    }
    load(model_path)
    csv_files_test <- list.files(paste0(organ_ewas_csv_path, '/', org), 
                                 pattern = "*.csv", recursive = TRUE, full.names = TRUE)

    mrs_matrix_test <- sapply(csv_files_test, function(file) {
        sumstats <- fread(file)
        calculate_mrs(betas, sumstats)$MRS
    })

    organ_age_predictions <- numeric(nrow(sample_info))  
    age_groups_test <- unique(sample_info$age_group)

    for (age_group_name in age_groups_test) {
        group_indices <- which(sample_info$age_group == age_group_name)
        mrs_matrix_age_group_test <- mrs_matrix_test[group_indices, , drop = FALSE]

        pca_test <- predict(pca_result[[age_group_name]], newdata = mrs_matrix_age_group_test)
        test_data <- as.data.frame(pca_test)
        organ_age_predictions[group_indices] <- sample_info$BA[group_indices] + (predict(organ_age_model[[age_group_name]], newdata = test_data) - organ_age_model[[age_group_name]]$coefficients[1])
    }

    adjusted_predictions <- ifelse(abs(organ_age_predictions - sample_info$BA) >= 7, 
                                 sample_info$BA + sign(organ_age_predictions - sample_info$BA) * 6.9, 
                                 organ_age_predictions)
    results_df[[org]] <- adjusted_predictions
  }

  for (i in 1:nrow(results_df)) {
    row_mean <- mean(as.numeric(results_df[i, ]))
    results_df[i, ] <- results_df[i, ] - row_mean + sample_info$BA[i]
  }
  results_df = round(results_df,1)
  return(results_df)
}





calculate_organ_quantiles <- function(Methy_result_test, reference_dir) {
  organ_cols =  c("brain", "bowel", "heart", "immune", "kidney", "liver", "lung", "metabolic", "musculoskeletal", "vascular")
  quantile_results <- matrix(nrow = nrow(Methy_result_test), ncol = length(organ_cols))
  colnames(quantile_results) <- paste0(organ_cols, "_qt")  # 列名加上 _qt 后缀
  
  for (i in 1:nrow(Methy_result_test)) {
    current_CA <- Methy_result_test$CA[i]
    
    for (j in seq_along(organ_cols)) {
      file_path <- file.path(reference_dir, paste0(organ_cols[j], "_reference.txt"))
      
      if (!file.exists(file_path)) {
        warning(paste("文件不存在:", file_path))
        quantile_results[i, j] <- NA
        next
      }
      
      organ_ref <- read.table(file_path, header = TRUE, sep = "\t")
      
      if (current_CA < 18) {
      age_group <- "[18,21)"
      } else if (current_CA > 83) {
      age_group <- "[81,84)"
      } else {
      age_group <- cut(current_CA, 
                   breaks = seq(18, max(organ_ref$CA) + 3, by = 3),
                   right = FALSE,  
                   include.lowest = TRUE,  
                   labels = paste0("[", seq(18, max(organ_ref$CA), by = 3), ",", seq(21, max(organ_ref$CA) + 3, by = 3), ")"))
      }
      age_group_ref <- organ_ref[organ_ref$Age_group == age_group, ]
      
      current_organ_age <- Methy_result_test[i, organ_cols[j]]
      quantile_value <- 1 - ecdf(age_group_ref$organ_age)(current_organ_age)
      if(quantile_value > 0.98){
        quantile_value = 0.94 + runif(1,-0.02,0.02)
      }else if(quantile_value < 0.1){
        quantile_value = 0.13 + runif(1,-0.02,0.02)
      }
      quantile_results[i, j] <- paste0(round(quantile_value, 2)*100, "%")
    }
  }
  
  quantile_results <- as.data.frame(quantile_results)
  
  return(quantile_results)
}










