library(optparse)
options(warn = 0)
option_list <- list(
  make_option(c("-c", "--config"), type = "character", default = NULL,
              help = "Path to the configuration file (config.R)", metavar = "FILE"),
  make_option("--idat-dir", type = "character", default = NULL,
              help = "Path to IDAT files", metavar = "DIR"),
  make_option("--sample-mapping-path", type = "character", default = NULL,
              help = "Path to sample mapping file", metavar = "FILE"),
  make_option("--sample-info-path", type = "character", default = NULL,
              help = "Path to sample info file", metavar = "FILE"),
  make_option("--output-dir", type = "character", default = NULL,
              help = "Output directory", metavar = "DIR")
)


opt_parser <- OptionParser(option_list = option_list)
args <- parse_args(opt_parser)


if (is.null(args$config) || is.null(args$`idat-dir`) ||
    is.null(args$`sample-mapping-path`) || is.null(args$`sample-info-path`) ||
    is.null(args$`output-dir`)) {
  print_help(opt_parser)
  stop("必须提供配置文件、IDAT目录、sample mapping 路径、sample info 路径和输出目录。", call. = FALSE)
}




IDAT_DIR <- args$`idat-dir`
SAMPLE_MAPPING_PATH <- args$`sample-mapping-path`
SAMPLE_INFO_PATH <- args$`sample-info-path`
OUTPUT_DIR <- args$`output-dir`
source(args$config)






suppressPackageStartupMessages({
  library(data.table)
  library(dnaMethyAge)
  library(dplyr)
  library(EpiDISH)
  library(tibble)
  library(sesame)
  library(readxl)
  library(BiocParallel)
  library(fmsb)
  library(scales)
  library(ggplot2)
  library(plotly)
  library(reticulate)
})

suppressWarnings(suppressMessages(library(meffonym)))






# 加载函数
{
source(CALCULATE_AGE_SCORE_PATH)
source(AGING_PLOT_PATH)
source(ORGAN_AGE_PREDICTION_PATH)
source(ORGAN_CPG_PATH)
source(TL_PREDICTION_PATH)
source(IMMUNE_SCORE_PATH)
source(EXPOSURE_PREDICTION_PATH)
source(BETA_IMPUTATION_PATH)
source(READ_IDAT_PATH)
source(EXPOSURE_CPG_PATH)

}



# 读取甲基化数据（读取-匹配-质控）
idat_re = process_methylation_data(
  idat_dir = IDAT_DIR, 
  sample_mapping_path = SAMPLE_MAPPING_PATH, 
  overlap_probe_path = OVERLAP_PROBE_PATH
)




ref_means <- read.csv(REF_MEANS_PATH, row.names = 1)
betas_all = data.frame(idat_re$betas_all, check.names = FALSE)
betas_all = beta_imputation(betas_all,ref_means)
betas_all_save = betas_all
if (ncol(betas_all) == 1) {
  col_name <- colnames(betas_all)[1]
  betas_all <- cbind(betas_all, dump = betas_all[, 1])
  colnames(betas_all)[1] <- col_name
}
betas = data.frame(idat_re$betas, check.names = FALSE)
betas = beta_imputation(betas,ref_means)
if (ncol(betas) == 1) {
  col_name <- colnames(betas)[1]
  betas <- cbind(betas, dump = betas[, 1])
  colnames(betas)[1] <- col_name
}
pval = data.frame(pval = idat_re$pval, check.names = FALSE)
colnames(pval) = 'pval'
detection_rates = data.frame(detection_rate = idat_re$detection_rates, check.names = FALSE)
qc_fail_samples = idat_re$qc_failed_samples





# 样本信息&出生年月日
sample_info = fread(SAMPLE_INFO_PATH)
sample_info = sample_info[which(sample_info$birth_day!=''),]
if(nrow(sample_info)>1){
  betas_all = betas_all[,as.character(sample_info$ID)]
  betas = betas[,as.character(sample_info$ID)]
}






cat('Step 1: Biological Age ...\n')
# 生物学年龄
Methy_result = data.frame(ID = colnames(betas))
Methy_result$detect_time = format(Sys.Date(), "%Y/%m/%d")
Methy_result$QC = rep('Y', nrow(Methy_result))
Methy_result$CA = round(as.numeric(difftime(sample_info$collect_day, sample_info$birth_day, units = "days")) / 365.25, 1)
suppressMessages(suppressWarnings(
  Methy_result$BA <- round(2.29 + 1.015 * methyAge(betas, clock = 'HorvathS2018')$mAge, 1)
))
Methy_result$AAF = round(Methy_result$BA / Methy_result$CA, 2)
Methy_result$HVI = calculate_age_score(Methy_result$CA, Methy_result$AAF)
cat('Done! ...\n')









cat('Step 2: Organ and System Age ...\n')
# 器官年龄
organ_age_by_age_group <- data.frame(BA = Methy_result$BA)
organ_age_by_age_group$age_group <- ifelse(organ_age_by_age_group$BA <= 18, "18-30", 
                                            as.character(cut(organ_age_by_age_group$BA, 
                                                              breaks = c(18, 30, 50, 70, Inf), 
                                                              labels = c("18-30", "30-50", "50-70", "70+"), 
                                                              right = FALSE)))
organ_age <- organ_age_prediction(organ_age_by_age_group, betas, ORGAN_PCA_MODEL_PATH, ORGAN_EWAS_PATH)
Methy_result <- cbind(Methy_result, organ_age)
organ_qt <- calculate_organ_quantiles(Methy_result, ORGAN_REFERENCE_PATH)
Methy_result <- cbind(Methy_result, organ_qt)
cat('Done! ...\n')






cat('Step 3: Telomere Length Evaluation ...\n')
# 端粒长度评估
suppressMessages(suppressWarnings(
  Methy_result$TL <- round(methyAge(betas, clock = 'LuA2019')$mAge, 2)
))
TL_group = fread(TL_GROUP_PATH)
Methy_result = cbind(Methy_result, calculate_TL_result(Methy_result$CA, Methy_result$TL, TL_group))
cat('Done! ...\n')






cat('Step 4: Immune Response Evaluation ...\n')
# 免疫功能评估
EPIC_reference = read.table(IMMUNE_REFERENCE_PATH)
EPIC_reference = as.matrix(EPIC_reference)
immune_result = calculate_immune_scores(betas_all, Methy_result$immune, Methy_result$BA, EPIC_reference)
Methy_result = cbind(Methy_result, immune_result)
cat('Done! ...\n')




cat('Step 5: Environmental Exposure Score ...\n')
# 环境暴露评估
Q_matrix = fread(EXPOSURE_Q_MATRIX_PATH) %>% column_to_rownames(var = "V1")
exposure_scores = suppressWarnings(predict_exposure(betas, exposure_ewas_dir = EXPOSURE_EWAS_DIR, Q_matrix = Q_matrix))
Methy_result = cbind(Methy_result, exposure_scores)












cat('Step 6: Get CPG value ...\n')
Methy_result = cbind(Methy_result, process_oragn_cpg_data(betas, ORGAN_CPG_TABLE_PATH, ORGAN_CPG_MEAN_PATH))
Methy_result = cbind(Methy_result,process_exposure_cpg_data(EXPOSURE_REFERENCE_MEANS_PATH,betas))
cat('Done! ...\n')







cat('Step 7: Combining Final result ...\n')
# 处理 QC 失败样本
if (exists("qc_fail_samples") && length(qc_fail_samples) > 0) {
  qc_fail_df <- data.frame(ID = qc_fail_samples, QC = "N")
  col_names <- colnames(Methy_result)
  na_cols <- setdiff(col_names, c("ID", "QC"))
  qc_fail_df[na_cols] <- NA
  Methy_result <- rbind(Methy_result, qc_fail_df)
}


if(nrow(sample_info)==1){
  Methy_result = Methy_result[1,]
}


cat('Generating Methy_Result ...\n')
write.csv(Methy_result,paste(OUTPUT_DIR,'Methy_result.csv',sep = '/'),fileEncoding = 'UTF-8',quote = FALSE,row.names = FALSE)
cat('Done!\n')

cat('Saving beta matrix ...\n')
write.csv(betas_all_save,paste(OUTPUT_DIR,'betas.csv',sep = '/'))
cat('Done!\n')

cat('Saving pval matrix ...\n')
write.csv(pval,paste(OUTPUT_DIR,'pval.csv',sep = '/'))
cat('Done!\n')

cat('Saving detection_rates ...\n')
write.csv(detection_rates, paste(OUTPUT_DIR,'detection_rates.csv',sep = '/'))
cat('Done!\n')


cat('Finish!\n')




