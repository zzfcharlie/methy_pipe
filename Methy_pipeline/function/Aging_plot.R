# 定义函数
# generate_sample_plots <- function(Methy_result, pheno, save_path) {
#   # 合并数据框，确保每个样本的 collect_day 可用
#   data <- merge(Methy_result, pheno[, c("ID", "collect_day")], by = "ID", all.x = TRUE)
  
#   # 创建虚拟数据框用于图例
#   line_legend <- data.frame(
#     slope = 1,
#     intercept = 0,
#     label = "日历年龄"
#   )
  
#   # 遍历每个样本
#   for (i in 1:nrow(data)) {
#     # 提取当前样本的数据
#     sample_data <- data[i, ]
    
#     # 以 CA 取整为中心，设置横坐标和纵坐标范围
#     ca_center <- round(sample_data$CA)  # CA 取整
#     age_range <- c(ca_center - 10, ca_center + 10)  # CA 加减 10
#     age_seq <- seq(age_range[1], age_range[2], by = 2)  # 间隔为 2
    
#     # 绘制图表
#     p <- ggplot(sample_data, aes(x = CA, y = BA)) +
#       geom_point(size = 3, color = "blue") +  # 点
#       geom_line(aes(color = "生物学年龄", linetype = "生物学年龄"), size = 1) +  # 生物学年龄折线
#       geom_abline(data = line_legend, aes(slope = slope, intercept = intercept, color = label, linetype = label),
#                   size = 1, show.legend = TRUE) +  # y=x 虚线
#       geom_hline(yintercept = age_seq, color = "grey", linetype = "dotted", size = 0.3) +  # 水平辅助线
#       geom_vline(xintercept = age_seq, color = "grey", linetype = "dotted", size = 0.3) +  # 垂直辅助线
#       # 注释
#       annotate("text", x = sample_data$CA, y = sample_data$BA - 2, label = sample_data$collect_day, hjust = 0.5, vjust = -0.5, 
#                size = 4, color = "darkblue", fontface = "bold") +  # 使用 collect_day 作为注释
#       labs(
#         x = "日历年龄 (岁)",
#         y = "生物学年龄 (岁)",
#         color = "",
#         linetype = ""
#       ) +
#       scale_color_manual(values = c("生物学年龄" = "blue", "日历年龄" = "red")) +
#       scale_linetype_manual(values = c("生物学年龄" = "solid", "日历年龄" = "dashed")) +
#       scale_x_continuous(limits = age_range, breaks = age_seq) +
#       scale_y_continuous(limits = age_range, breaks = age_seq) +
#       # 主题设置
#       theme_minimal(base_family = "PingFangSC") +
#       theme(
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         axis.line = element_blank(),
#         legend.position = "top",
#         legend.title = element_text(size = 12, face = "bold"),
#         legend.text = element_text(size = 12, face = "bold"),
#         legend.background = element_rect(fill = "white", color = "black", size = 0.8),
#         axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 15)),
#         axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 15)),
#         axis.text = element_text(size = 12, face = "bold"),
#         aspect.ratio = 1
#       ) +
#       # 添加手动轴线和箭头
#       annotate("segment", x = age_range[1], xend = age_range[2], y = age_range[1], yend = age_range[1], 
#                arrow = arrow(length = unit(0.2, "cm")), color = "black", size = 1) +
#       annotate("segment", x = age_range[1], xend = age_range[1], y = age_range[1], yend = age_range[2], 
#                arrow = arrow(length = unit(0.2, "cm")), color = "black", size = 1)
    
#     # 保存图表，文件名为 ID_动态生物学年龄.png
#     file_name <- file.path(save_path, paste0(sample_data$ID, "_动态生物学年龄.png"))
#     ggsave(file_name, plot = p, width = 6, height = 6, dpi = 300, bg = 'white')
#   }
# }


generate_sample_plots <- function(Methy_result, pheno, save_path) {

  if (!dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
  }
  # 合并数据框，确保每个样本的 collect_day 可用
  data <- merge(Methy_result, pheno[, c("ID", "collect_day")], by = "ID", all.x = TRUE)
  
  # 创建虚拟数据框用于图例
  line_legend <- data.frame(
    slope = 1,
    intercept = 0,
    label = "日历年龄"
  )
  
  # 遍历每个样本
  for (i in 1:nrow(data)) {
    # 提取当前样本的数据
    sample_data <- data[i, ]
    
    # 以 CA 取整为中心，设置横坐标和纵坐标范围
    ca_center <- round(sample_data$CA)  # CA 取整
    age_range <- c(ca_center - 10, ca_center + 10)  # CA 加减 5
    age_seq <- seq(age_range[1], age_range[2], by = 2)  # 间隔为 1
    
    # 绘制图表
    p <- ggplot(sample_data, aes(x = CA, y = BA)) +
      geom_point(size = 3, color = "blue") +  # 点
      geom_abline(data = line_legend, aes(slope = slope, intercept = intercept, color = label, linetype = label),
                  size = 1, show.legend = TRUE) +  # y=x 虚线
      geom_hline(yintercept = age_seq, color = "grey", linetype = "dotted", size = 0.3) +  # 水平辅助线
      geom_vline(xintercept = age_seq, color = "grey", linetype = "dotted", size = 0.3) +  # 垂直辅助线
      # 注释
      annotate("text", x = sample_data$CA, y = sample_data$BA - 2, label = sample_data$collect_day, hjust = 0.5, vjust = -0.5, 
               size = 4, color = "darkblue", fontface = "bold") +  # 使用 collect_day 作为注释
      labs(
        x = "日历年龄 (岁)",
        y = "生物学年龄 (岁)",
        color = "",
        linetype = ""
      ) +
      scale_color_manual(values = c("生物学年龄" = "blue", "日历年龄" = "red")) +
      scale_linetype_manual(values = c("生物学年龄" = "solid", "日历年龄" = "dashed")) +
      scale_x_continuous(limits = age_range, breaks = age_seq) +
      scale_y_continuous(limits = age_range, breaks = age_seq) +
      # 主题设置
      theme_minimal(base_family = "PingFangSC") +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_blank(),
        legend.position = "top",
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.background = element_rect(fill = "white", color = "black", size = 0.8),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 15)),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 15)),
        axis.text = element_text(size = 12, face = "bold"),
        aspect.ratio = 1
      ) +
      # 添加手动轴线和箭头
      annotate("segment", x = age_range[1], xend = age_range[2], y = age_range[1], yend = age_range[1], 
               arrow = arrow(length = unit(0.2, "cm")), color = "black", size = 1) +
      annotate("segment", x = age_range[1], xend = age_range[1], y = age_range[1], yend = age_range[2], 
               arrow = arrow(length = unit(0.2, "cm")), color = "black", size = 1)
    
    # 保存图表，文件名为 ID_动态生物学年龄.png
    file_name <- file.path(save_path, paste0(sample_data$ID, "_动态生物学年龄.png"))
    ggsave(file_name, plot = p, width = 6, height = 6, dpi = 300, bg = 'white')
  }
}





# 定义函数：为每个样本绘制并保存HVI仪表图
generate_HVI_plots <- function(Methy_result, output_dir) {
  # 确保输出目录存在
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # 循环处理每个样本
  for (i in 1:nrow(Methy_result)) {
    sample_data <- Methy_result[i, ]
    sample_id <- sample_data$ID
    hvi_value <- sample_data$HVI
    output_path <- file.path(output_dir, paste0(sample_id, "_HVI.png"))
    # 创建仪表图
    fig <- plot_ly(
      domain = list(x = c(0, 1), y = c(0, 1)),
      value = hvi_value, # 使用样本的HVI值
      title = list(text = ""),
      type = "indicator",
      mode = "gauge+number+delta",
      gauge = list(
        axis = list(
          range = list(0, 100), # 设置范围
          tickvals = seq(0, 100, by = 20), # 显示主要刻度
          ticktext = as.character(seq(0, 100, by = 20)), # 确保文字显示
          tickfont = list(size = 20) # 调整刻度字体大小
        ),
        bar = list(
          color = "blue"
        ),
        steps = list(
          list(range = c(0, 20), color = "lightblue"),
          list(range = c(20, 40), color = "skyblue"),
          list(range = c(40, 60), color = "deepskyblue"),
          list(range = c(60, 80), color = "dodgerblue"),
          list(range = c(80, 100), color = "royalblue")
        ),
        borderwidth = 5, # 调整外圈的宽度
        bordercolor = "gray" # 设置外圈颜色
      )
    ) %>%
    layout(
      margin = list(l = 50, r = 50, t = 50, b = 50)  # 调整边距
    )
    
    # 使用save_image保存图表为PNG
    save_image(fig, file = output_path, width = 800, height = 800, engine = "kaleido")
  }
}
