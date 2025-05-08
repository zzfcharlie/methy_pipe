calculate_immune_scores <- function(betas_all, immune_age, BA, epic_reference) {
  pro = 0.3
  adjust_CD4T_CD8T_ratio <- function(cell_pro, pro) {
    CD4T_CD8T_ratio <- cell_pro$CD4T / cell_pro$CD8T
    adjust_idx <- which(CD4T_CD8T_ratio < 1.5)
    cell_pro$CD4T[adjust_idx] <- cell_pro$CD4T[adjust_idx] + cell_pro$CD8T[adjust_idx] * pro
    cell_pro$CD8T[adjust_idx] <- cell_pro$CD8T[adjust_idx] * (1 - pro)
    return(cell_pro)
  }

  calculate_immune_score <- function(cell_pro) {
    CD4T_CD8T_ratio <- cell_pro$CD4T / cell_pro$CD8T
    normal_ratio_min <- 1.5
    normal_ratio_max <- 4.0
    lower_ratio_min <- 0.5
    upper_ratio_max <- 8.0
    scores <- numeric(length(CD4T_CD8T_ratio))
    
    for (i in seq_along(CD4T_CD8T_ratio)) {
      ratio <- round(CD4T_CD8T_ratio[i],1)
      if (ratio < normal_ratio_min) {
        score <- -25 + ((ratio - lower_ratio_min) / (normal_ratio_min - lower_ratio_min)) * 20
      } else if (ratio <= normal_ratio_max) {
        score <- -5 + ((ratio - normal_ratio_min) / (normal_ratio_max - normal_ratio_min)) * 10
      } else {
        score <- 5 + ((ratio - normal_ratio_max) / (upper_ratio_max - normal_ratio_max)) * 20
      }
      scores[i] <- max(min(score, 25), -25)
    }
    
    conclusion <- cut(scores,
                      breaks = c(-Inf, -15, -10, -5, 5, 10, 15, Inf),
                      labels = c("免疫重度抑制", "免疫中度抑制", "免疫轻度抑制", 
                                "免疫正常", "免疫轻度激活", "免疫中度激活", 
                                "免疫重度激活"),
                      right = FALSE)
    
    return(data.frame(
      CD4T_CD8T_ratio = round(CD4T_CD8T_ratio,1),
      Immune_Score = round(scores)
    ))
  }

  calc_score <- function(P, L, U) {
    Median <- (L + U) / 2
    Range <- (U - L) / 2
    if (L <= P && P <= U) {
      return(1 - abs(P - Median) / Range)
    } else if (P > U) {
      return(exp(-0.3 + -abs(P - U) / Range))
    } else {
      return(exp(-0.3 + -abs(L - P) / Range))
    }
  }

  
  EPIC_reference <- as.matrix(epic_reference)
  cell_pro_o <- data.frame(epidish(beta.m = betas_all, ref.m = EPIC_reference)$estF)
  adjust_immune_cells <- function(immune_data) {
    target_values <- c(CD8T = 0.1, CD4T = 0.18, NK = 0.05, Bcell = 0.08, Mono = 0.05, Neu = 0.5)
    threshold_cd8t <- 0.07
    threshold_nk <- 0.04
    adjusted_data <- apply(immune_data, 1, function(sample) {
      if (sample["CD8T"] < threshold_cd8t || sample["NK"] < threshold_nk) {
        diff <- target_values - sample
        sample <- sample + diff * 0.5
        sample <- pmax(sample, 0)
        sample <- sample / sum(sample)
      }
      return(sample)
    })
      adjusted_data <- t(adjusted_data)
      colnames(adjusted_data) <- colnames(immune_data)
      return(as.data.frame(adjusted_data))
    }
  cell_pro_o = adjust_immune_cells(cell_pro_o)
  cell_pro <- adjust_CD4T_CD8T_ratio(cell_pro = cell_pro_o, pro = pro)
  immune_scores <- calculate_immune_score(cell_pro)
  rownames(immune_scores) <- colnames(betas_all)

  limits <- list(
    CD8T = c(L = 0.07, U = 0.12),
    CD4T = c(L = 0.15, U = 0.25),
    NK = c(L = 0.05, U = 0.1),
    Bcell = c(L = 0.05, U = 0.1),
    Mono = c(L = 0.05, U = 0.1),
    Neu = c(L = 0.4, U = 0.6)
  )

  scores <- t(sapply(1:nrow(cell_pro), function(i) {
    sapply(colnames(cell_pro), function(cell_type) {
      calc_score(cell_pro[i, cell_type], limits[[cell_type]]["L"], limits[[cell_type]]["U"])
    })
  }))
  colnames(scores) <- colnames(cell_pro)
  rownames(scores) <- rownames(cell_pro)

  weights <- list(
    Defense = c(CD8T = 0.2, NK = 0.2, Neu = 0.6),
    Surveillance = c(CD4T = 0.2, Mono = 0.3, Neu = 0.5),
    Homeostasis = c(CD4T = 0.3, Mono = 0.2, Bcell = 0.2, Neu = 0.3)
  )

  scores_defense <- apply(scores, 1, function(row) sum(row[names(weights$Defense)] * weights$Defense))
  scores_surveillance <- apply(scores, 1, function(row) sum(row[names(weights$Surveillance)] * weights$Surveillance))
  scores_homeostasis <- apply(scores, 1, function(row) sum(row[names(weights$Homeostasis)] * weights$Homeostasis))

  final_scores <- data.frame(
    Defense = scores_defense,
    Surveillance = scores_surveillance,
    Homeostasis = scores_homeostasis
  )

  immune_residual <- immune_age - BA
  score_residual <- (7 - immune_residual) / (2 * 7)
  final_scores$Residual = score_residual

  final_scores$Defense_Score = 0.7 * final_scores$Defense + 0.3 * final_scores$Residual
  final_scores$Surveillance_Score = 0.7 * final_scores$Surveillance + 0.3 * final_scores$Residual
  final_scores$Homeostasis_Score = 0.7 * final_scores$Homeostasis + 0.3 * final_scores$Residual

  # 计算Level
  final_scores$Defense_Level <- cut(
    final_scores$Defense_Score,
    breaks = c(-Inf, 0.25, 0.50, 0.75, Inf),
    labels = c(2, 3, 4, 5),
    right = TRUE
  )

  final_scores$Surveillance_Level <- cut(
    final_scores$Surveillance_Score,
    breaks = c(-Inf, 0.25, 0.50, 0.75, Inf),
    labels = c(2, 3, 4, 5),
    right = TRUE
  )

  final_scores$Homeostasis_Level <- cut(
    final_scores$Homeostasis_Score,
    breaks = c(-Inf, 0.25, 0.50, 0.75, Inf),
    labels = c(2, 3, 4, 5),
    right = TRUE
  )

  result_df <- cbind(
    (round(cell_pro,3)*100),
    immune_scores,
    final_scores[, c("Defense_Level", "Surveillance_Level", "Homeostasis_Level")]
  )
  return(result_df)
}