# process_exposure_betas <- function(betas_all, exposure_means_dir, output_dir) {
#     # 创建输出目录
#     if (!dir.exists(output_dir)) {
#         dir.create(output_dir, recursive = TRUE)
#     }
#      exposure_ref_files <- list.files(exposure_means_dir, pattern = "_exposure_reference_means\\.csv$", full.names = TRUE)
#     # 遍历每个暴露参考文件
#     for (ref_file in exposure_ref_files) {
#         # 读取参考均值文件
#         case_mean <- fread(ref_file)
        
#         # 提取case名称
#         case_name <- gsub("_exposure_reference_means\\.csv$", "", basename(ref_file))
        
#         # 获取对应的β值并转换为dataframe
#         beta_case <- as.data.frame(betas_all[case_mean$CpG, ])
        
#         # 添加CpG列作为第一列
#         beta_case <- cbind(CpG = rownames(beta_case), beta_case)
        
#         # 将Mean列添加到第二列
#         beta_case <- cbind(beta_case[, 1, drop = FALSE], Mean = case_mean$Mean, beta_case[, -1, drop = FALSE])
        
#         # 保存结果
#         output_file <- file.path(output_dir, paste0(case_name, '_beta_values.csv'))
#         write.csv(beta_case, output_file, row.names = FALSE, quote = FALSE)
#     }
# }





process_exposure_cpg_data = function(exposure_ref_mean_path,betas){
  combine_ex_re <- NULL
  exposure_cpg_files = list.files(EXPOSURE_REFERENCE_MEANS_PATH,full.names = TRUE)
  exposure_names = c('air_pollution','alcohol','metal','sleep','smoking','stress','UVR')
  for(i in 1:7){
    file = exposure_cpg_files[i]
    ewas = read.csv(file)
    ex_cpg_mean = matrix(ewas$Mean,ncol(betas),length(ewas$CpG),byrow = TRUE)
    colnames(ex_cpg_mean) = paste(exposure_names[i],'_',ewas$CpG,'_mean',sep = '')
    rownames(ex_cpg_mean) = colnames(betas)
    ex_cpg_value = t(betas[ewas$CpG,])
    colnames(ex_cpg_value) = paste(exposure_names[i],'_',ewas$CpG,'_value',sep = '')
    if(is.null(combine_ex_re)){
      combine_ex_re = cbind(ex_cpg_mean,ex_cpg_value)
    }else{
      combine_ex_re = cbind(combine_ex_re,cbind(ex_cpg_mean,ex_cpg_value))
    }
  }
  return(combine_ex_re)
}