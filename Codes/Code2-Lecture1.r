# Set seed for reproducibility
set.seed(123)

# Generate a large population (true distribution)
population <- rnorm(10000, mean = 50, sd = 10)  # Mean 50, SD 10

# Define a range of sample sizes
sample_sizes <- seq(2, 500, by = 5)  # From 2 to 500, step 5
num_trials <- 100  # Number of samples per sample size

# Store mean standard deviations for each sample size
mean_biased_sds <- numeric(length(sample_sizes))
mean_unbiased_sds <- numeric(length(sample_sizes))

# Loop through different sample sizes
for (i in seq_along(sample_sizes)) {
  n <- sample_sizes[i]
  biased_sds <- numeric(num_trials)
  unbiased_sds <- numeric(num_trials)
  for (j in 1:num_trials) {
    sample_data <- sample(population, n)
    # Biased standard deviation (dividing by n)
    biased_sds[j] <- sqrt(sum((sample_data - mean(sample_data))^2) / n)
    # Unbiased standard deviation (dividing by n-1) - R's default sd()
    unbiased_sds[j] <- sd(sample_data)
  }
  # Store the mean values for this sample size
  mean_biased_sds[i] <- mean(biased_sds)
  mean_unbiased_sds[i] <- mean(unbiased_sds)
}

# Plot the results
plot(sample_sizes, mean_biased_sds, type = "o", col = "red", pch = 16, lwd = 2,
     xlab = "Sample Size (n)", ylab = "Average Standard Deviation",
     main = "Convergence of Biased and Unbiased Standard Deviation")
lines(sample_sizes, mean_unbiased_sds, type = "o", col = "blue", pch = 16, lwd = 2)

# Add legend
legend("bottomright", legend = c("Biased SD (dividing by n)", "Unbiased SD (dividing by n-1)"),
       col = c("red", "blue"), pch = 16, lwd = 2)