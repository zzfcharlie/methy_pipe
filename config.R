# config.R

# 工作目录
WORK_DIR <- file.path(getwd(), "Methy_pipeline")

# 函数文件路径
FUNCTION_DIR <- file.path(WORK_DIR, 'function')
CALCULATE_AGE_SCORE_PATH <- file.path(FUNCTION_DIR, 'calculate_age_score.R')
AGING_PLOT_PATH <- file.path(FUNCTION_DIR, 'Aging_plot.R')
ORGAN_AGE_PREDICTION_PATH <- file.path(FUNCTION_DIR, 'organ_age_prediction.R')
TL_PREDICTION_PATH <- file.path(FUNCTION_DIR, 'TL_prediction.R')
IMMUNE_SCORE_PATH <- file.path(FUNCTION_DIR, 'immune_score.R')
EXPOSURE_PREDICTION_PATH <- file.path(FUNCTION_DIR, 'Exposure_prediction.R')
BETA_IMPUTATION_PATH <- file.path(FUNCTION_DIR, 'beta_imputation.R')
READ_IDAT_PATH <- file.path(FUNCTION_DIR, 'Read_idat.R')
EXPOSURE_CPG_PATH <- file.path(FUNCTION_DIR, 'Exposure_CpG.R')
ORGAN_CPG_PATH <- file.path(FUNCTION_DIR, 'organ_cpg.R')







# 其他文件路径
OVERLAP_PROBE_PATH <- file.path(WORK_DIR, 'files/Overlap_probes/overlap_probe.txt')
REF_MEANS_PATH <- file.path(WORK_DIR,'files/Ref_Means/ref_means.csv')
ORGAN_PCA_MODEL_PATH <- file.path(WORK_DIR, 'files/Organ_age/pca_model')
ORGAN_EWAS_PATH <- file.path(WORK_DIR, 'files/Organ_age/organ_EWAS')
ORGAN_REFERENCE_PATH <- file.path(WORK_DIR, 'files/Organ_age/organ_reference')
ORGAN_CPG_TABLE_PATH <- file.path(WORK_DIR, 'files/Organ_cpg/cpg_table.csv')
ORGAN_CPG_MEAN_PATH <- file.path(WORK_DIR, 'files/Organ_cpg/cpg_ref_mean.csv')
TL_GROUP_PATH <- file.path(WORK_DIR, 'files/TL/TL_group_reference.txt')
IMMUNE_REFERENCE_PATH <- file.path(WORK_DIR, 'files/Immune_reference/IDOLOptimizedCpGs.txt')
EXPOSURE_Q_MATRIX_PATH <- file.path(WORK_DIR, 'files/Exposure_Q_matrix/Exposure_Q_matrix.txt')
EXPOSURE_EWAS_DIR <- file.path(WORK_DIR, 'files/Exposure_EWAS')
EXPOSURE_REFERENCE_MEANS_PATH <- file.path(WORK_DIR, 'files/Exposure_reference_means')

