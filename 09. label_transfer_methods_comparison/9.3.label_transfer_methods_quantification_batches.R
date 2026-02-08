library(readxl)
library(ggplot2)
library(dplyr)
library(stringr)


setwd("D:/xueying-work/afterPHD/Liu-lab/project1-embryo benchmarking/2024April/manuscript/NCB_rebuttal_2025_june/code/20250812_label_transfer_method_quantification_v3")
outputDir = getwd()

#load label transfer benchmarking result
scib <- read.csv("dataset_scib_benchmarking.csv")

# Step 1: extract type
scib$type <- ifelse(str_detect(scib$Lineage_key, "lineage"), "lineage", "reanno")

# Step 2: extract method（去掉 _lineage_pred / _reanno_pred）
scib$method <- gsub("_lineage_pred|_reanno_pred", "", scib$Lineage_key)

# Step 3: convert  wide to long
scib_long <- scib %>%
  select(Dataset, method, type, Silhouette, NMI, ARI) %>%
  pivot_longer(cols = c(Silhouette, NMI, ARI),
               names_to = "Metric",
               values_to = "Value")


# Step 4: ensure factor type
scib_long$Dataset <- as.factor(scib_long$Dataset)
scib_long$method <- as.factor(scib_long$method)
scib_long$type <- as.factor(scib_long$type)
scib_long$Metric <- as.factor(scib_long$Metric)

# Step 5: define plot setting
theme_custom <- function() {
  theme_classic() +
    theme(
      axis.line = element_line(color = "black", size = 0.8),
      panel.border = element_rect(fill = NA, color = "black", size = 0.8),
      panel.background = element_blank(),
      legend.position = "right",
      plot.title = element_text(hjust = 0.5),
      strip.background = element_rect(fill = "white", color = "black", size = 0.8),
      strip.text = element_text(size = 10)
    )
}

plot_metrics_grouped_bar <- function(data, metric_title) {
  ggplot(
    data = subset(data, Metric == metric_title),
    aes(x = method, y = Value, fill = type)
  ) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    facet_wrap(~ Dataset, scales = "free_y") +  # 每个 Dataset 一个子图
    labs(
      x = "Prediction Method",
      y = metric_title,
      title = paste0("Benchmark: ", metric_title)
    ) +
    theme_custom() +
    scale_fill_brewer(palette = "Set2") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Step 6: generate plots
plot_silhouette <- plot_metrics_grouped_bar(scib_long, "Silhouette")
plot_nmi <- plot_metrics_grouped_bar(scib_long, "NMI")
plot_ari <- plot_metrics_grouped_bar(scib_long, "ARI")

print(plot_silhouette)
print(plot_nmi)
print(plot_ari)

# Step 7: save plots
ggsave("Silhouette_grouped_barplot.pdf", plot_silhouette, width = 10, height = 8)
ggsave("NMI_grouped_barplot.pdf", plot_nmi, width = 10, height = 8)
ggsave("ARI_grouped_barplot.pdf", plot_ari, width = 10, height = 8)


##########funky plot###################
library(funkyheatmap)
library(dplyr)
library(tibble)
library(ggplot2)
library(RColorBrewer)

# Step 1: reorganize data + pivot_wider
scib_wide <- scib_long %>%
  select(Dataset, method, type, Metric, Value) %>%
  pivot_wider(names_from = Metric, values_from = Value)

# Step 2: split by type 
scib_lineage <- scib_wide %>% filter(type == "lineage")
scib_reanno <- scib_wide %>% filter(type == "reanno")

# Step 3: define label_top_3 
label_top_3 <- function(scores) {
  ranks <- rank(-scores, ties.method = "min")
  labels <- ifelse(ranks <= 3, as.character(ranks), "")
  return(labels)
}

# Step 4: calculate mean + add label
process_data <- function(df) {
  df %>%
    group_by(method) %>%
    summarise(
      Silhouette = mean(Silhouette, na.rm = TRUE),
      NMI = mean(NMI, na.rm = TRUE),
      ARI = mean(ARI, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    ungroup() %>%
    mutate(id = as.character(1:n())) %>%
    mutate(
      label_Silhouette = label_top_3(Silhouette),
      label_NMI = label_top_3(NMI),
      label_ARI = label_top_3(ARI)
    ) %>%
    as.data.frame()
}

# Step 5: separate lineage and reanno
scib_lineage_plot <- process_data(scib_lineage)
scib_reanno_plot <- process_data(scib_reanno)

# Step 6: define column_info
scib_column_info <- tribble(
  ~id,                   ~id_color,     ~name,             ~geom,   ~group,               ~options,
  "id",                  NA,            "Rank",            "text",  "Method",             list(hjust = 0),
  "method",              NA,            "Method",          "text",  "Method",             list(hjust = 0, width = 5),
  "Silhouette",          "Silhouette",  "Silhouette",      "bar",   "Bio conservation",   list(palette = "olives", width = 1.5, draw_outline = FALSE),
  "label_Silhouette",    NA,            NA,                "text",  "Bio conservation",   list(hjust = 0.1, overlay = TRUE, palette = "reds"),
  "NMI",                 "NMI",         "NMI",             "bar",   "Bio conservation",   list(palette = "navys", width = 1.5, draw_outline = FALSE),
  "label_NMI",           NA,            NA,                "text",  "Bio conservation",   list(hjust = 0.1, overlay = TRUE, palette = "greens"),
  "ARI",                 "ARI",         "ARI",             "bar",   "Bio conservation",   list(palette = "oranges2", width = 1.5, draw_outline = FALSE),
  "label_ARI",           NA,            NA,                "text",  "Bio conservation",   list(hjust = 0.1, overlay = TRUE, palette = "blues")
)

# Step 7: update column_groups
scib_column_groups <- tribble(
  ~group,               ~palette,     ~level1,
  "Method",             "black",      "Method",
  "Bio conservation",   "greys",      "Bio conservation"
)

# Step 8: row_info
row_info <- function(df) {
  data.frame(id = df$id, group = character(nrow(df)))
}

# Step 9: define palettes
palettes <- list(
  navys = rev(colorRampPalette(c("#748EBB", "white"))(9)),
  olives = rev(colorRampPalette(c("#5F9856", "#F7FCF5"))(9)),
  oranges2 = RColorBrewer::brewer.pal(9, "Oranges"),
  black = c("black", "black"),
  greys = rev(colorRampPalette(c("#808080", "white"))(9))
)

# Step 10: define legends
legends <- list(
  list(
    title = "Silhouette rank",
    palette = "olives",
    geom = "rect",
    labels = c("0", " ", "0.5", " ", "1"),
    size = c(1, 1, 1,  1, 1)
  ),
  list(
    title = "NMI rank",
    palette = "navys",
    geom = "rect",
    labels = c("0", " ", "0.5", " ", "1"),
    size = c(1, 1, 1, 1, 1)
  ),
  list(
    title = "ARI rank",
    palette = "oranges2",
    geom = "rect",
    labels = c("0", " ", "0.5", " ", "1"),
    size = c(1, 1, 1, 1, 1)
  )
)

# Step 11: plot
generate_plot <- function(data, filename) {
  row_info_df <- row_info(data)
  
  J <- funky_heatmap(
    data = data,
    column_info = scib_column_info,
    column_groups = scib_column_groups,
    row_info = row_info_df,
    palettes = palettes,
    legends = legends,
    position_args = position_arguments(col_annot_offset = 4),
    scale_column = FALSE
  )
  
  ggsave(filename, plot = J, width = 8, height = 7, units = "in", dpi = 300)
  return(J)
}

# Step 12: generate plots
plot_lineage <- generate_plot(scib_lineage_plot, "lineage_funkyplot_mean.pdf")
plot_reanno <- generate_plot(scib_reanno_plot, "reanno_funkyplot_mean.pdf")

print(plot_lineage)
print(plot_reanno)