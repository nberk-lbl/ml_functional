# Trying to replicate the umap plot with a poplar genome

embeddings = ''
embeddings_umap = ''

find_bounds <- function(x) {
  # Sort the values to handle outliers
  x <- sort(x)
  
  # Define a cutoff for extreme outliers based on quantiles
  lower_quantile <- quantile(x, 0.05) # Lower 5th percentile
  upper_quantile <- quantile(x, 0.95) # Upper 95th percentile

  # Calculate the interquartile range for the main bulk of the data
  iqr <- IQR(x)
  
  # Define the main bounds by trimming extreme values based on IQR
  lower_bound <- max(min(x), lower_quantile - 1.5 * iqr)
  upper_bound <- min(max(x), upper_quantile + 1.5 * iqr)
  
  return(c(lower_bound, upper_bound))
}

plot(embeddings_umap $layout[,1], embeddings_umap $layout[,2], cex=.4, pch=19, 
    xlim=find_bounds(embeddings_umap $layout[,1]), 
    ylim=find_bounds(embeddings_umap $layout[,2]))
