process_methylation_data <- function(idat_dir, sample_mapping_path, overlap_probe_path, detection_threshold = 0.95, workers = 4) {
    idat_files <- list.files(idat_dir, pattern = "\\.idat$", full.names = TRUE)
    idat_files <- unique(gsub("_(Grn|Red).idat$", "", idat_files))    
    sample_mapping <- readxl::read_excel(sample_mapping_path)
    idat_samples <- basename(idat_files)
    match_idx <- match(paste(sample_mapping$芯片号, sample_mapping$文件名, sep = "_"), idat_samples)
    valid_idx <- !is.na(match_idx)
    matched_files <- idat_files[match_idx[valid_idx]]
    sample_bc <- sample_mapping$样本编号[valid_idx]
    param <- MulticoreParam(workers = workers, progressbar = TRUE)
    cat('Calculating QC statistics...\n')
    calcStats <- sesame::openSesame(matched_files, func = sesameQC_calcStats, BPPARAM = param)
    if (!is.list(calcStats)) {
        calcStats <- list(calcStats)
    }
    detection_rate <- sapply(calcStats, function(sample) sample@stat$frac_dt)
    detection_rate = data.frame(detection_rate)
    cat('Done!/n')
    param <- BiocParallel::MulticoreParam(workers = workers, progressbar = TRUE)
    cat('Calculating detection p-value...\n')
    sdfs <- openSesame(matched_files, func = NULL, BPPARAM = param, collapseToPfx = TRUE)
    bpstop(param)
    if (inherits(sdfs, "SigDF")) {
      detection_pvals <- list(pOOBAH(sdfs, return.pval = TRUE))
    } else if (is.list(sdfs) && all(sapply(sdfs, function(x) inherits(x, "SigDF")))) {
      detection_pvals <- lapply(sdfs, function(sdf) {
        pOOBAH(sdf, return.pval = TRUE)
      })
    } else {
      stop("openSesame did not return valid SigDF objects.")
    }
    pval_matrix <- data.frame(detection_pvals)
    cat('Done!\n')
    bpstop(param)
    qc_y_idx <- which(detection_rate >= detection_threshold)
    qc_n_idx <- which(detection_rate < detection_threshold)
    param <- MulticoreParam(workers = workers, progressbar = TRUE)
    cat('Calculating betas matrix...\n')
    betas_all <- openSesame(matched_files[qc_y_idx], collapseToPfx = TRUE, BPPARAM = param)
    bpstop(param)
    cat('Done!\n')
    overlap_probe <- fread(overlap_probe_path, header = FALSE)
    # betas_all <- beta_imputation(betas_all)
    betas_all = data.frame(betas = betas_all)
    colnames(betas_all) <- sample_bc[qc_y_idx]
    betas <- betas_all[overlap_probe$V1, ]
    betas = data.frame(betas = betas)
    colnames(betas) <- sample_bc[qc_y_idx]
    rownames(betas) = overlap_probe$V1
    qc_failed_samples <- sample_bc[qc_n_idx]    
    return(list(
        betas_all = betas_all, 
        betas = betas,  
        pval = pval_matrix,
        detection_rates = detection_rate,
        qc_failed_samples = qc_failed_samples 
    ))
}




read_idat <- function(idat_dir, overlap_probe_path, detection_threshold = 0.95, workers = 4) {
    idat_files <- list.files(idat_dir, pattern = "\\.idat$", full.names = TRUE)
    idat_files <- unique(gsub("_(Grn|Red).idat$", "", idat_files))
    param <- MulticoreParam(workers = workers, progressbar = TRUE)
    cat('Calculating QC statistics...\n')
    calcStats <- openSesame(idat_files, func = sesameQC_calcStats, BPPARAM = param)
    bpstop(param)
    cat('Done!\n')
    detection_rate <- sapply(calcStats, function(sample) {
        return(sample@stat$frac_dt)
    })

    param <- MulticoreParam(workers = workers, progressbar = TRUE)
    cat('Calculating detection p-value...\n')
    sdfs <- openSesame(idat_files, func = NULL, BPPARAM = param,collapseToPfx = TRUE)
    detection_pvals <- lapply(sdfs, function(sdf) {
        pOOBAH(sdf, return.pval = TRUE)
    })
    pval_matrix <- do.call(cbind, detection_pvals)
    cat('Done!\n')
    bpstop(param)
    qc_y_idx <- which(detection_rate >= detection_threshold)
    qc_n_idx <- which(detection_rate < detection_threshold)
    param <- MulticoreParam(workers = workers, progressbar = TRUE)
    cat('\nCalculating betas matrix...\n')
    betas_all <- openSesame(idat_files[qc_y_idx], collapseToPfx = TRUE, BPPARAM = param)
    bpstop(param)
    cat('Done!\n')
    overlap_probe <- fread(overlap_probe_path, header = FALSE)
    betas_all <- beta_imputation(betas_all)
    betas <- betas_all[overlap_probe$V1, ]
    qc_failed_samples <- colnames(betas)[qc_n_idx]
    return(list(
        betas_all = betas_all,  
        betas = betas,
        pval = pval_matrix,
        qc_failed_samples = qc_failed_samples  
    ))
}