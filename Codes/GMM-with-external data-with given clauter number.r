# Load necessary libraries
library(MASS)
library(mclust)
library(ggplot2)

# Load the saved data into a data.frame
data <- read.table("~/Desktop/dendrites.txt", header = TRUE)

# Ensure the data has exactly two columns
if (ncol(data) != 2) {
  stop("Error: The dataset must have exactly two columns.")
}

colnames(data) <- c("X1", "X2")

# Fit GMM with automatic selection of clusters
gmm_model <- Mclust(data)

# Alternatively, you can fit GMM with a specified number of clusters (e.g., G = 7)
gmm_model <- Mclust(data, G = 7)

# Add cluster labels
data$cluster <- as.factor(gmm_model$classification)

# Plot data with cluster colors
ggplot(data, aes(x = X1, y = X2, color = cluster)) +
  geom_point(size = 2, alpha = 0.7) +
  theme_minimal() +
  labs(title = "GMM Clustering", x = "X1", y = "X2") +
  theme(legend.title = element_blank())