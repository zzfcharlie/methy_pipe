predict_exposure <- function(betas, exposure_ewas_dir, Q_matrix) {
  case_names <- c(
    "air_pollution",
    "alcohol",
    "arsenic",
    "cadmium",
    "lead",
    "manganese",
    "mercury",
    "sleep",
    "smoking",
    "stress",
    "UVR"
  )
  
  ewas_files <- list.files(exposure_ewas_dir, pattern = "_exposure\\.csv$", full.names = TRUE)
  
  mrs_matrix <- matrix(nrow = ncol(betas), ncol = length(ewas_files))
  colnames(mrs_matrix) <- case_names
  rownames(mrs_matrix) <- colnames(betas)
  
  for (i in seq_along(ewas_files)) {
    ewas <- fread(ewas_files[i])
    common_cpg <- intersect(rownames(betas), ewas$CpG)
    mrs_matrix[, i] <- colSums(betas[common_cpg,] * ewas$Coef[seq(length(common_cpg))], na.rm = TRUE)
  }
  betas = as.matrix(betas)
  mrs_matrix[,'smoking'] = meffonym.score(betas,'langdon-candidate-ever-vs-never')$score
  model_name_5cpgs <- '5cpgs'
  variables_5cpgs <- c("cg00252472", "cg06690548", "cg06846495", "cg12825509", "cg18282388")
  coefficients_5cpgs <- c(-0.7771057, -4.7079934, -0.3341005, -1.4103232, 0.1192672)
  description_5cpgs <- "5 CpG site model"
  meffonym.add.model(name = model_name_5cpgs, variables = variables_5cpgs, coefficients = coefficients_5cpgs, description = description_5cpgs)
  betas = as.matrix(betas)
  drinking_score <- suppressWarnings(meffonym.score(betas, '5cpgs')$score)
  mrs_matrix[,'alcohol'] = drinking_score
  score_exposure <- function(mrs_matrix, Q_matrix) {
    mrs_matrix <- as.data.frame(mrs_matrix)
    Q_matrix <- as.data.frame(Q_matrix)
    scores <- data.frame(matrix(nrow = nrow(mrs_matrix), ncol = ncol(mrs_matrix)))
    colnames(scores) <- colnames(mrs_matrix)
    rownames(scores) <- rownames(mrs_matrix)
    for (col_name in colnames(mrs_matrix)) {
      column_data <- mrs_matrix[[col_name]]
      q_lower <- Q_matrix[Q_matrix$Case == col_name, "Q_lower"]
      q_median <- Q_matrix[Q_matrix$Case == col_name, "Q_median"]
      q_upper <- Q_matrix[Q_matrix$Case == col_name, "Q_upper"]
      min_val <- Q_matrix[Q_matrix$Case == col_name, "min_val"]
      max_val <- Q_matrix[Q_matrix$Case == col_name, "max_val"]
      
      if (q_lower == min_val) q_lower <- min_val + 1e-6
      if (q_upper == max_val) q_upper <- max_val - 1e-6
      
      scores[[col_name]] <- sapply(column_data, function(value) {
        if (value == q_median) {
          return(80)
        } else if (value < q_lower) {
          return(98 - 8 * (q_lower - value) / (q_lower - min_val))
        } else if (value < q_upper) {
          return(90 - 20 * (value - q_lower) / (q_upper - q_lower))
        } else {
          return(70 - 8 * (value - q_upper) / (max_val - q_upper))
        }
      })
    }
    return(scores)
  }


  score_result <- score_exposure(mrs_matrix, Q_matrix)
  colnames(score_result) = case_names
  
  weights <- c(
    arsenic = 0.05,   # 砷
    cadmium = 0.15,   # 镉
    lead = 0.45,      # 铅
    manganese = 0.20, # 锰
    mercury = 0.15    # 汞
  )

  metal_columns <- names(weights)
  if (!all(metal_columns %in% colnames(score_result))) {
    stop("部分金属列未找到，请检查权重的列名是否与数据匹配！")
  }
  
  score_result$Metal <- rowSums(sapply(metal_columns, function(metal) {
    score_result[[metal]] * weights[metal]
  }))
  
  score_result <- score_result[, !(colnames(score_result) %in% metal_columns)]
  adjustments <- c(5, 0, -10, 0, -8, 5, 8)
  score_result <- t(t(score_result) + adjustments)
  score_result <- ifelse(round(score_result) >= 98, 96, round(score_result)) 
  score_result <- ifelse(round(score_result) <= 62, 65, round(score_result)) 
  

  numeric_part <- as.data.frame(round(score_result))

  ref_score <- 80 + adjustments
  ref_score_m <- matrix(ref_score, nrow(score_result), 7, byrow = TRUE)
  # risk_assessment <- ifelse(score_result >= ref_score_m, '相对较小', '相对较大')


  delta_values <- c(3, 4, 4, 4, 4, 3, 3)
  delta_m <- matrix(delta_values, nrow(score_result), 7, byrow = TRUE)

  risk_assessment <- ifelse(
    score_result < (ref_score_m - delta_m), '相对较大',
    ifelse(score_result > (ref_score_m + delta_m), '相对较小', '程度适中')
  )

  char_part <- as.data.frame(risk_assessment, stringsAsFactors = FALSE)
  colnames(char_part) <- paste0(colnames(score_result), '_conclusion')
  exposure_re <- cbind(numeric_part, char_part)
  return(exposure_re)
}








###所有case的###
plot_exposure_radar <- function(Methy_result, save_folder = NULL) {
  # 确保保存文件夹存在
  if (!is.null(save_folder) && !dir.exists(save_folder)) {
    dir.create(save_folder, recursive = TRUE)
  }
  
  # 遍历每一行（每个 case）
  for (i in 1:nrow(Methy_result)) {
    # 提取当前行的暴露数据和 ID
    exposure_data <- Methy_result[i, c("smoking", "alcohol", "air_pollution", "UVR", "Metal", "stress", "sleep")]
    case_id <- Methy_result[i, "ID"]
    
    # 添加最大值、最小值和参考人群
    max_values <- rep(100, ncol(exposure_data))
    min_values <- rep(60, ncol(exposure_data))
    reference_values <- c(80, 80, 85, 85, 88, 72, 70)
    # 创建适合 radarchart 的数据框
    radar_data <- rbind(max_values, min_values, reference_values, exposure_data)
    rownames(radar_data) <- c("最大值", "最小值", "参考人群", "您的分数")
    colnames(radar_data) <- c("烟草", "酒精", "空气污染", "紫外线", "重金属", "压力", "睡眠质量")
    
    # 如果指定了保存路径，则保存为 PNG 文件
    if (!is.null(save_folder)) {
      # 设置 PNG 文件路径
      png_file <- file.path(save_folder, paste0(case_id, "_暴露健康评分" , ".png"))
      
      # 打开 PNG 设备
      png(png_file, width = 2500, height = 2500, res = 300)
    }
    
    # 绘制雷达图
    radarchart(radar_data, title = "",
               pcol = c("#1F77B4", "#FF7F0E"),                     # 设置线条颜色
               pfcol = scales::alpha(c(NA, "#FF7F0E"), 0.3),       # 蓝色不填充，橙色填充
               plwd = c(2, 4),                                     # 蓝色线宽为2，橙色线宽为3
               plty = 1,                                           # 设置线条样式
               cglcol = "grey", cglty = 1, cglwd = 0.8,            # 设置网格线颜色和样式
               axislabcol = "grey",                                # 设置轴标签颜色
               caxislabels = seq(60, 100, by = 10),                # 设置刻度标签为 60, 70, 80, 90, 100
               vlabels = colnames(radar_data),                     # 显示变量标签
               vlcex = 1.1)                                          # 标签大小
    
    # 添加图例
    legend("topright", legend = c("参考人群暴露健康评分", "您的暴露健康评分"), 
           col = c("#1F77B4", "#FF7F0E"),                  # 图例颜色与线条颜色一致
           pch = 16,                                       # 图例方块符号
           pt.cex = 1.5,                                   # 方块大小
           bty = "n",                                      # 不显示图例边框
           text.col = "black",                             # 图例文字颜色
           cex = 1.1)
    # 如果指定了保存路径，则关闭 PNG 设备
    if (!is.null(save_folder)) {
      dev.off()
    }
  }
}




plot_exposure_per_case <- function(data, case) {
  # 复制数据，用于突出显示指定的暴露变量
  highlight_data <- data
  highlight_data["参考人群", setdiff(colnames(data), case)] <- 60
  highlight_data["您的分数", setdiff(colnames(data), case)] <- 60
  
  # 绘制基础雷达图
  # radarchart(data,
            #  pcol = c("#D3D3D3", "#D3D3D3"),                 # 设置线条颜色为灰色
            #  pfcol = scales::alpha(c("#D3D3D3", "#D3D3D3"), 0.3), # 设置填充颜色为灰色，带透明度
            #  plwd = c(2, 2),                                 # 设置线条宽度
            #  plty = 1,                                       # 设置线条样式
            #  cglcol = "grey", cglty = 1, cglwd = 0.8,        # 设置网格线颜色和样式
            #  axislabcol = "grey",                            # 设置轴标签颜色
            #  caxislabels = seq(60, 100, by = 10),            # 设置刻度标签为 60, 70, 80, 90, 100
            #  vlabels = colnames(data),                       # 显示变量标签
            #  vlcex = 1.1)                                    # 调整变量标签字体大小
  
  # 在同一图形设备中绘制突出显示的暴露变量
  # par(new = TRUE) # 在同一图形设备上绘制新图
  radarchart(highlight_data,
             pcol = c("#1F77B4", "#FF7F0E"),                 # 设置线条颜色
             pfcol = scales::alpha(c(NA, "#FF7F0E"), 0.5),   # 蓝色不填充，橙色填充
             plwd = c(2, 4),                                 # 设置线条宽度
             plty = 1,                                       # 设置线条样式
             cglcol = "grey", cglty = 1, cglwd = 0.8,        # 设置网格线颜色和样式
             axislabcol = "grey",                            # 设置轴标签颜色
             caxislabels = seq(60, 100, by = 10),            # 设置刻度标签为 60, 70, 80, 90, 100
             vlabels = colnames(data),                       # 显示变量标签
             vlcex = 1.1)                                    # 调整变量标签字体大小
  
  # 添加图例
  legend("topright", legend = c("参考人群暴露健康评分", "您的暴露健康评分"), 
         col = c("#1F77B4", "#FF7F0E"),                      # 图例颜色与线条颜色一致
         pch = 16,                                           # 图例方块符号
         pt.cex = 1.5,                                       # 方块大小
         bty = "n",                                          # 不显示图例边框
         text.col = "black",                                 # 图例文字颜色
         cex = 1.1) 
}





# 定义绘制并保存图像的函数
plot_exposure_plots_by_cases <- function(data, target_path) {
  # 定义英文到中文的变量名映射
  var_order <- c("烟草", "酒精", "空气污染", "紫外线", "重金属", "压力", "睡眠质量")
  var_mapping <- c(
    smoking = "烟草",
    alcohol = "酒精",
    air_pollution = "空气污染",
    UVR = "紫外线",
    Metal = "重金属",
    stress = "压力",
    sleep = "睡眠质量"
  )
  
  # 遍历每个 ID
  for (id in data$ID) {
    # 创建以 ID 命名的文件夹
    id_folder <- file.path(target_path, id)
    if (!dir.exists(id_folder)) {
      dir.create(id_folder, recursive = TRUE)
    }
    
    # 提取当前 ID 的数据
    target_row <- data[data$ID == id, ]
    
    # 构建符合要求的数据框
    exposure_data <- data.frame(
      row.names = c("最大值", "最小值", "参考人群", "您的分数"),
      烟草 = c(100, 60, 80, target_row$smoking),
      酒精 = c(100, 60, 80, target_row$alcohol),
      空气污染 = c(100, 60, 85, target_row$air_pollution),
      紫外线 = c(100, 60, 85, target_row$UVR),
      重金属 = c(100, 60, 88, target_row$Metal),
      压力 = c(100, 60, 72, target_row$stress),
      睡眠质量 = c(100, 60, 70, target_row$sleep)
    )
    
    
    # 遍历每个暴露变量
    for (var in names(var_mapping)) {
      # 获取中文变量名
      var_chinese <- var_mapping[var]
      
      # 设置图像保存路径
      file_name <- paste0(id, "_", var_chinese, ".png")
      file_path <- file.path(id_folder, file_name)
      
      # 打开 PNG 设备
      png(file_path, width = 2500, height = 2500, res = 300)
      
      # 调用绘图函数
      plot_exposure_per_case(exposure_data, case = var_chinese)
      
      # 关闭 PNG 设备
      dev.off()
    }
  }
}

