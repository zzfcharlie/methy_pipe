calculate_TL_result <- function(CA, TL, TL_group) {
  result_df <- data.frame(
    TL_age = numeric(length(CA)),
    quantile = character(length(CA)),
    Q25 = numeric(length(CA)),
    Q75 = numeric(length(CA)),
    stringsAsFactors = FALSE
  )
  age_groups <- seq(18, max(TL_group$Age) + 3, by = 3)
  group_stats <- TL_group %>%
    mutate(Age_group = cut(Age, breaks = age_groups, right = FALSE, include.lowest = TRUE)) %>%
    group_by(Age_group) %>%
    summarise(
      group_median_age = median(Age),
      TL_median = median(TL)
    )
  for (i in seq_along(CA)) {
    if (CA[i] < 18) {
      age_group <- factor("[18,21)", levels = levels(group_stats$Age_group))
    } else if (CA[i] > 78) {
        age_group <- factor("[75,78]", levels = levels(group_stats$Age_group))
    } else {
        age_group <- cut(CA[i], 
                        breaks = age_groups, 
                        right = FALSE, 
                        include.lowest = TRUE,
                        labels = levels(group_stats$Age_group))
    }
    group_info <- group_stats %>% filter(Age_group == age_group)
    TL_diff <- TL[i] - group_info$TL_median
    TL_age <- group_info$group_median_age - (TL_diff / 0.033)
    TL_group[which(TL_group$Age_group == '[75,78)'),'Age_group'] = '[75,78]'
    group_TL <- TL_group[Age_group == age_group, TL]
    quantile_num <- ecdf(group_TL)(TL[i])
    quantile_num = ifelse(quantile_num > 0.98, 0.95 + runif(1,-0.02,0.02), ifelse(quantile_num < 0.1, 0.14 + runif(1,-0.02,0.02), quantile_num))
    quantile <- paste(round(quantile_num * 100), "%", sep = "")
    if (quantile == "0%") {
      quantile <- paste(round(runif(1, 0.01, 0.05) * 100), "%", sep = "")
    }
    Q25 <- quantile(group_TL, 0.25)
    Q75 <- quantile(group_TL, 0.75)
    result_df$TL_age[i] <- round(TL_age)
    result_df$quantile[i] <- quantile
    result_df$Q25[i] <- round(Q25, 2)
    result_df$Q75[i] <- round(Q75, 2)
  }
  return(result_df)
}






# plot_TL_scatter <- function(TL_group, Methy_result, output_dir) {
#   # 检查并创建输出目录
#   if (!dir.exists(output_dir)) {
#     dir.create(output_dir, recursive = TRUE)
#   }

#   # 对每个样本绘制散点图
#   for (i in 1:nrow(Methy_result)) {
#     p <- ggplot(TL_group, aes(x = Age, y = TL)) +
#       # 蓝点：参考人群端粒长度
#       geom_point(aes(color = "参考人群端粒长度"), size = 1) +
#       # 红线：线性拟合（不显示在图例中）
#       geom_smooth(
#         method = "lm", formula = y ~ x,
#         color = "black", size = 1.5, linetype = "dashed", se = FALSE, show.legend = FALSE # 虚线
#       ) +
#       # 红点：您的端粒长度
#       geom_point(aes(x = Methy_result$CA[i], y = Methy_result$TL[i], color = "您的端粒长度"),
#         size = 2.5, show.legend = TRUE
#       ) + # 显示在图例中
#       theme_minimal(base_size = 15) +
#       theme(
#         panel.border = element_rect(color = "black", fill = NA, size = 1),
#         axis.line = element_line(color = "black"),
#         aspect.ratio = 0.5,
#         axis.text = element_text(size = 15),
#         axis.title = element_text(size = 15),
#         legend.position = "bottom", # 将图例放在底部
#         legend.title = element_blank() # 去掉图例标题
#       ) +
#       labs(
#         x = "日历年龄(岁)",
#         y = "端粒长度(kb)"
#       ) +
#       scale_color_manual(
#         name = "", # 图例标题为空
#         values = c("参考人群端粒长度" = "blue", "您的端粒长度" = "red"), # 指定颜色
#         guide = guide_legend(override.aes = list(size = 2)) # 统一图例中点的大小
#       ) +
#       xlim(0, 100) + # 设置 x 轴范围
#       ylim(6, 8.5) # 设置 y 轴范围

#     # 保存图像，使用ID作为文件名
#     ggsave(
#       p,
#       filename = file.path(output_dir, paste0(Methy_result$ID[i], "_端粒长度.png")),
#       width = 8, height = 6, units = "in", dpi = 300, bg = "white"
#     )
#   }
# }




