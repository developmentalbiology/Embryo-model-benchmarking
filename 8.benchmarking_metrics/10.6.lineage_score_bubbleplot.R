library(funkyheatmap)
library(dplyr)
library(tibble)
library(ggplot2)

##palettes
palettes <- list(
  features = c(FULL = "black", HVG = "black"),
  black = c("black", "black"),
  orange = rev(colorRampPalette(c("#FFA13F", "white"))(9)),
  purples = rev(colorRampPalette(c("#BA7FB5", "white"))(9)),
  navys = rev(colorRampPalette(c("#748EBB", "white"))(9)),
  olives = rev(colorRampPalette(c("#5F9856", "#F7FCF5"))(9)),
  oranges2 = RColorBrewer::brewer.pal(9, "Oranges")
)


# Get all _lineage.xlsx files in the current folder
# Define the folder containing the CSV files
folder_folder <- "./output_folder"

# Get all _lineage.csv files in the output folder
file_list <- list.files(path = folder_folder, pattern = "_lineage\\.csv$", full.names = TRUE)

metric_read <- list()

# Loop through the list of files and read them into the list
for (file in file_list) {
  # Read the Excel file into a data frame
  df <- read.csv(file, row.names = 1)
  
  # Convert the data frame to a matrix and store it in the list
  metric_read[[file]] <- as.matrix(df)
}

# Remove the "./output_folder/" and "_lineage.csv" from the file names
names(metric_read) <- gsub("./output_folder/", "", file_list)
names(metric_read) <- gsub("_lineage.csv$", "", names(metric_read))


# Loop through each dataset in metric_read
metrics_sum <- lapply(1:length(metric_read), function(i) {
  
  # Convert the matrix into a data frame with numeric values
  metrics <- as.data.frame(apply(metric_read[[i]], 2, as.numeric))
  
  # Calculate the sum of columns 2 to 7 (adjusted to the correct columns)
  metrics$sum <- rowSums(metrics[, 1:6])  # Adjusted to columns 2 to 7
  
  # Replace NA values with 0 in the entire dataframe
  metrics[is.na(metrics)] <- 0
  
  # Normalize the metrics$sum column to the range 0 to 1
  metrics$sum <- metrics$sum / 6
  
  # Create a new data frame with row names and the scaled sum
  metrics_sum <- data.frame(Row.names = rownames(metric_read[[i]]),
                            scaled.sum = metrics$sum)
  
  # Set column name as the lineage name (from the names of the metric_read list)
  colnames(metrics_sum)[2] <- names(metric_read)[i]
  
  return(metrics_sum)
})

# Combine results into a final data frame, merging by Row.names (to avoid duplicate columns)
metrics_sum <- Reduce(function(x, y) merge(x, y, by = "Row.names", all = TRUE), metrics_sum)

############################generate plot############################################################
# Reshape the matrix into long format using pivot_longer
metrics_long <- metrics_sum %>%
  pivot_longer(cols = -Row.names, 
               names_to = "Lineage", 
               values_to = "Value")

# Create a bubble plot heatmap using ggplot2
p <- ggplot(metrics_long, aes(x = Lineage, y = Row.names, size = Value, color = Value)) +
  geom_point() +  # Add points with varying size
  scale_size_continuous(range = c(1, 10)) +  # Control size of bubbles
  scale_color_gradientn(colors = c("blue", "cyan", "yellow", "orange", "red"),
              values = scales::rescale(c(0, 0.4,0.6,0.8, 1))) +  # Rescale color breaks so yellow is at 0.7
  #scale_color_gradientn(colors = c("blue", "yellow", "red")) +  # Use a custom color gradient
  #scale_color_viridis_c() +  # Use viridis color scale
  #scale_color_viridis_c(option = "inferno") +  # Use inferno color scale
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
  labs(x = "Lineage", y = "Row Names", size = "Value", title = "Bubble Plot Heatmap of Metrics") +
  theme(legend.position = "bottom")

p

# Save the plot as a PDF
ggsave(
  "lineage_benchmarking.pdf", 
  plot = p, 
  device = "pdf", 
  width = 5, 
  height = 5)
