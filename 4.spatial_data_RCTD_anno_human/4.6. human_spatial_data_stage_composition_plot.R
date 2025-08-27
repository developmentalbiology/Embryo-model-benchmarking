
setwd("D:/xueying-work/afterPHD/Liu-lab/project1-embryo benchmarking/2024April/manuscript/NCB_rebuttal_2025_june/code/20250724_human_refined_anno_RCTD_v3")
outputDir = getwd()

library(Seurat)
library(tidyverse)
library(Matrix)
library(dplyr)

cs8 <- readRDS("human_cs8_reanno.Rds")
cs7 <- readRDS("human_cs7_reanno.Rds")


#assign colors
godsnot_102 <- c(
  "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
  "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF",
  "#997D87", "#5A0007", "#809693", "#6A3A4C", "#1B4400", "#4FC601", "#3B5DFF",
  "#4A3B53", "#FF2F80", "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92",
  "#FF90C9", "#B903AA", "#D16100", "#DDEFFF", "#000035", "#7B4F4B", "#A1C299",
  "#300018", "#0AA6D8", "#013349", "#00846F", "#372101", "#FFB500", "#C2FFED",
  "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09", "#00489C", "#6F0062",
  "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66", "#885578",
  "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
  "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F",
  "#938A81", "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757",
  "#C8A1A1", "#1E6E00", "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C",
  "#772600", "#D790FF", "#9B9700", "#549E79", "#FFF69F", "#201625", "#72418F",
  "#BC23FF", "#99ADC0", "#3A2465", "#922329", "#5B4534", "#FDE8DC", "#404E55",
  "#0089A3", "#CB7E98", "#A4E804", "#324E72"
)


# Define the ordered labels
final_anno_labels = c('TE', 'CTB_1','CTB_2', 'STB_1', 'STB_2', 'STB_3', 'EVT_1', 'EVT_2',
                      'Epiblast_1','Epiblast_2','Epiblast_3','Ectoderm',
                      'Amniontic.epi','Amniontic.ectoderm',
                      'PGC',
                      'Primitive.streak',
                      'Neuromesodermal.progenitor',
                      'Neural.crest', 'Neural.ectoderm.forebrain', 'Neural.ectoderm.hindbrain', 'Neural.ectoderm.midbrain','Spinal.cord',
                      'Paraxial.mesoderm','Emergent.mesoderm','Pre-somatic.mesoderm','Somite', 'Rostral.mesoderm', 'Lateral.plate.mesoderm_1',
                      'Lateral.plate.mesoderm_2','Lateral.plate.mesoderm_3','Cardiac.mesoderm','Amniotic.mesoderm','Exe.meso.progenitor','YS.mesoderm_1', 'YS.mesoderm_2',
                      'Hypoblast_1', 'Hypoblast_2', 'AVE', 'VE', 'YS.endoderm',
                      'DE','Gut',
                      'Notochord',
                      'Hemogenic.endothelial.progenitor','Endothelium','Erythroid','Primitive.megakaryocyte','Myeloid.progenitor'
)


# Ensure the number of colors matches the number of clusters
color_palette <- setNames(godsnot_102[1:length(final_anno_labels)], final_anno_labels)

# ============================================================================
# 1. compare cs7 reanno
# ============================================================================

df <- as.data.frame(table(cs7$clusters, cs7$final_anno))
colnames(df)[1:2] <- c("orig_anno", "reanno")

df$reanno <- factor(df$reanno, levels = final_anno_labels)


df <- do.call(rbind, lapply(unique(df$orig_anno), function(x) {
  
  df <- df[which(df$orig_anno==x),]
  df$perct <- df$Freq/sum(df$Freq)
  
  return(df)
}))

p <- ggplot(data=df, aes(x=orig_anno, y=perct, fill=reanno)) +
  geom_bar(stat="identity") +
  ggtitle("reanno") +
  scale_fill_manual(values = color_palette) + 
  theme(text=element_text(size=10),
        axis.text.x = element_text(angle = 90, hjust = 1))  

# plot stage percentage for all
pdf(paste0(outputDir, "/human_cs7_spatial_reanno_proportion.pdf"), width = 12, height = 4)
p
dev.off()

# ============================================================================
# 2. compare cs8 reanno 
# ============================================================================

df <- as.data.frame(table(cs8$clusters, cs8$final_anno))
colnames(df)[1:2] <- c("orig_anno", "reanno")

df$reanno <- factor(df$reanno, levels = final_anno_labels)


df <- do.call(rbind, lapply(unique(df$orig_anno), function(x) {
  
  df <- df[which(df$orig_anno==x),]
  df$perct <- df$Freq/sum(df$Freq)
  
  return(df)
}))

p <- ggplot(data=df, aes(x=orig_anno, y=perct, fill=reanno)) +
  geom_bar(stat="identity") +
  ggtitle("reanno") +
  scale_fill_manual(values = color_palette) + 
  theme(text=element_text(size=10),
        axis.text.x = element_text(angle = 90, hjust = 1))  

# plot stage percentage for all
pdf(paste0(outputDir, "/human_cs8_spatial_reanno_proportion.pdf"), width = 12, height = 4)
p
dev.off()

# ============================================================================
# 3. cs7 reanno composition
# ============================================================================

cs7$stage <- "CS7"
df <- as.data.frame(table(cs7$stage,cs7$final_anno))
colnames(df)[1:2] <- c("stage", "reanno")

df$reanno <- factor(df$reanno, levels = final_anno_labels)


df <- do.call(rbind, lapply(unique(df$stage), function(x) {
  
  df <- df[which(df$stage==x),]
  df$perct <- df$Freq/sum(df$Freq)
  
  return(df)
}))

p <- ggplot(data=df, aes(x=stage, y=perct, fill=reanno)) +
  geom_bar(stat="identity") +
  ggtitle("reanno") +
  scale_fill_manual(values = color_palette) + 
  theme(text=element_text(size=10),
        axis.text.x = element_text(angle = 90, hjust = 1))  

# plot stage percentage for all
pdf(paste0(outputDir, "/human_cs7_spatial_reanno_dataset_proportion.pdf"), width = 6, height = 4)
p
dev.off()

# ============================================================================
# 4. cs8 reanno composition
# ============================================================================

cs8$stage <- "CS8"
df <- as.data.frame(table(cs8$stage,cs8$final_anno))
colnames(df)[1:2] <- c("stage", "reanno")

df$reanno <- factor(df$reanno, levels = final_anno_labels)


df <- do.call(rbind, lapply(unique(df$stage), function(x) {
  
  df <- df[which(df$stage==x),]
  df$perct <- df$Freq/sum(df$Freq)
  
  return(df)
}))

p <- ggplot(data=df, aes(x=stage, y=perct, fill=reanno)) +
  geom_bar(stat="identity") +
  ggtitle("reanno") +
  scale_fill_manual(values = color_palette) + 
  theme(text=element_text(size=10),
        axis.text.x = element_text(angle = 90, hjust = 1))  

# plot stage percentage for all
pdf(paste0(outputDir, "/human_cs8_spatial_reanno_dataset_proportion.pdf"), width = 6, height = 4)
p
dev.off()
