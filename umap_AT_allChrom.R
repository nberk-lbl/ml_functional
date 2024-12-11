# Trying to replicate the umap plot with a poplar genome
library(umap)
library(ggplot2)

plot_umap_with_labels <- function(umap, labels, threshold = 1.5, margin_ratio = 0.05) {
  # Ensure data frames have the same number of rows
  if (nrow(umap) != nrow(labels)) stop("Data frames must have the same number of rows.")
  
  
  # Filter rows where labels[,3] is not "ambiguous"
  filtered_umap <- umap[labels[,3] != "ambiguous", ]
  filtered_labels <- labels[labels[,3] != "ambiguous", ]
  
  # Combine umap and labels into a single data frame for plotting
  combined <- data.frame(
    x = filtered_umap[,1],
    y = filtered_umap[,2],
    label = filtered_labels[,3]
  )
  
  # Calculate IQR and axis limits
  x_iqr <- IQR(combined$x, na.rm = TRUE)
  y_iqr <- IQR(combined$y, na.rm = TRUE)
  x_median <- median(combined$x, na.rm = TRUE)
  y_median <- median(combined$y, na.rm = TRUE)
  
  x_lower <- x_median - threshold * x_iqr
  x_upper <- x_median + threshold * x_iqr
  y_lower <- y_median - threshold * y_iqr
  y_upper <- y_median + threshold * y_iqr
  
  # Add a margin to ensure no points are clipped
  x_margin <- margin_ratio * (x_upper - x_lower)
  y_margin <- margin_ratio * (y_upper - y_lower)
  x_lower <- x_lower - x_margin
  x_upper <- x_upper + x_margin
  y_lower <- y_lower - y_margin
  y_upper <- y_upper + y_margin
  
  
  color_palette <- brewer.pal(n = 8, name = "Set2")
  
  # Create ggplot
ggplot(combined, aes(x = x, y = y, color = label)) +
  geom_point(size = 10e-6) +  # Small points in the plot
  scale_x_continuous(limits = c(x_lower, x_upper)) +
  scale_y_continuous(limits = c(y_lower, y_upper)) +
  scale_color_manual(values = color_palette) +
  labs(
    title = "TAIR10 chrom 4",
    x = "UMAP 1",
    y = "UMAP 2",
    color = "Label"
  ) +
  theme_minimal() +
  theme(
    legend.key.size = unit(1.5, "cm"),    # Larger key box
    legend.text = element_text(size = 12),  # Readable legend text
    legend.title = element_text(size = 14)  # Readable legend title
  ) +
  guides(
    color = guide_legend(
      override.aes = list(size = 5, shape = 15)  # Use squares (shape = 15), larger size
    )
  )
}


# Example Usage
set.seed(12)

setwd("/clusterfs/jgi/groups/science/homes/nberkowitz/gpn/ml_functional/")

#labels = read.csv("data/Arabidopsis_thaliana.TAIR10.60.allchrom.csv", head=FALSE)
#embeddings = read.csv(paste("data/TAIR10_all_avg_embeddings.tsv", sep =""), sep="\t", head=TRUE)

labels = read.csv("data/TAIR10.allchrom.anno.csv", head=FALSE)
embeddings = read.csv(paste("data/TAIR10.with_rc_all_avg_embeddings.tsv", sep =""), sep="\t", head=TRUE)

embeddings.umap = umap(embeddings)
plot_umap_with_labels(embeddings.umap, labels)
write.csv(embeddings.umap$layout, file = "data/TAIR10_all_4_UMAP.csv", row.names = FALSE, quote = TRUE)



