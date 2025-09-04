library(readxl)
library(devtools)
library(RColorBrewer) 
library(reshape2)
library(plotly)
library(Seurat)
library(dplyr)
library(Matrix)
library(gplots)
library(tidyverse)


setwd("D:/xueying-work/afterPHD/Liu-lab/project1-embryo benchmarking/2024April/manuscript/NCB_rebuttal_2025_june/code/20250806_embryo_model_benchmarking/cell_type_composition_analysis")
outputDir = getwd()

df <- read.csv("composition_correlation_matrix.csv")

reference_data <- c("E6_zhou","E6_IVC_zhou","E6_IVC_xiang", "E7_IVC_xiang","E8_IVC_xiang","E8_IVC_zhou",
                    "E9_IVC_mole","E9_IVC_xiang", "E10_IVC_Ai", "E10_IVC_xiang","E10_IVC_zhou","E11_IVC_Ai",
                    "E11_IVC_mole", "E12_IVC_Ai","E12_IVC_xiang","E12_IVC_zhou", "E13_IVC_Ai","E14_IVC_Ai",
                    "E14_IVC_xiang","E14_IVC_zhou", "CS7_tyser","CS9_Yuan_2025", "CS10_zeng")

embryo_model <- c( "D11_Liu", 
                   "D14_Rowan" , "D21_Rowan", 
                   "D2_Hislop","D3_Hislop","D4_Hislop","D5_Hislop", "D21_Hislop",  
                   "D4_Oldak", "D6_Oldak",  "D8_Oldak" ,
                   "D4_Pedroza" , "D8_Oldak" ,
                   "D4_Weatherbee","D6_Weatherbee" ,"D8_Weatherbee" ,
                   "D8_Ai_model",
                   "embryo_model_Zheng_2019",
                   "embryo_model_zheng_2022"  )

df <- df[which(df$stage_orig_ident %in% reference_data),]
df <- df[,-which(colnames(df) %in% reference_data)]



mat <- as.matrix(df[, -1])  
rownames(mat) <- df$stage_orig_ident 

# reorder rows
mat <- mat[reference_data, , drop = FALSE] 
# reorder columns
mat <- mat[, embryo_model, drop = FALSE] 


library(ggplot2)

# convert matrix to long format
bubble_data <- as.data.frame(as.table(mat))
colnames(bubble_data) <- c("Reference", "Model", "Correlation")

# create bubble plot
bubble_plot <- ggplot(bubble_data, aes(x = Model, y = Reference, 
                                       size = abs(Correlation),  # define bubble size to reflect correlation value
                                       color = Correlation)) +  # define color to reflect correlation value
  geom_point(alpha = 1) +  # transparency
  scale_size_continuous(range = c(1, 5),  # bubble size
                        name = "|Correlation|") +
  scale_color_gradient2(low = "#08519c", mid = "white", high = "#cb181d", 
                        midpoint = 0,  
                        limits = c(-0.25, 1), 
                        name = "Correlation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        panel.grid.major = element_line(color = "grey90"),
        legend.position = "right") +
  labs(x = "Embryo Models", y = "Reference Data") +
  coord_fixed(ratio = 0.8)  

# view plot
print(bubble_plot)

# save plot
ggsave(paste0(outputDir, "/stage_composition_correlation_bubble.pdf"), 
       plot = bubble_plot, 
       width = 6.5, height = 6.5)


# organize table

mat2 <- mat["CS7_tyser",, drop=FALSE]
mat2 <- as.matrix(mat2) 

# extract dataset names
extract_dataset_name <- function(colname) {
  parts <- unlist(strsplit(colname, "_"))
  if (length(parts) > 1) {
    paste(parts[-1], collapse = "_")
  } else {
    colname 
  }
}

dataset_names <- sapply(colnames(mat2), extract_dataset_name)

print(dataset_names)

# use the max value
result <- data.frame(
  CS7_correlation = sapply(unique(dataset_names), function(ds) {
    # summarizes the highest correlation value between each data and CS7_tyser
    cols <- which(dataset_names == ds)
    vals <- mat2["CS7_tyser", cols]
    if(length(vals) == 0 || all(is.na(vals))) NA else max(vals, na.rm = TRUE)
  }),
  row.names = unique(dataset_names)
)

print(result)

write.csv(result, file = "stage_composition_correlation.csv")



