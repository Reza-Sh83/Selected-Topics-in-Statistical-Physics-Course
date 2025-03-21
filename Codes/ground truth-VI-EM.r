# Load necessary libraries
library(MASS)
library(mixtools)

# Set ground truth parameters for the synthetic data
pi_true <- c(0.6, 0.4)
mu_true <- c(2, 5)
sigma_true <- c(1, 1.5)

# Generate synthetic data based on the ground truth
set.seed(42)  # For reproducibility
n <- 10000  # Number of data points
data <- c(
  rnorm(n * pi_true[1], mean = mu_true[1], sd = sigma_true[1]),
  rnorm(n * pi_true[2], mean = mu_true[2], sd = sigma_true[2])
)

# Number of components in the GMM
K <- 2

# Function for EM Algorithm
EM_GMM <- function(data, K) {
  # Initial parameters (random guess)
  pi <- rep(1/K, K)
  mu <- c(2.0, 5.0)
  sigma <- rep(1.0, K)
  log_likelihood <- function(data, pi, mu, sigma) {
    sum(log(sapply(data, function(x) sum(pi * dnorm(x, mu, sigma)))))
  }
  # E-step: Compute responsibilities
  e_step <- function(data, pi, mu, sigma) {
    gamma <- matrix(0, nrow = length(data), ncol = K)
    for (n in 1:length(data)) {
      for (k in 1:K) {
        gamma[n, k] <- pi[k] * dnorm(data[n], mu[k], sigma[k])
      }
      gamma[n, ] <- gamma[n, ] / sum(gamma[n, ])
    }
    return(gamma)
  }
  # M-step: Update parameters
  m_step <- function(gamma, data) {
    N_k <- colSums(gamma)
    pi <- N_k / length(data)
    mu <- colSums(gamma * data) / N_k
    sigma <- sqrt(colSums(gamma * (data - mu)^2) / N_k)
    return(list(pi = pi, mu = mu, sigma = sigma))
  }
  # EM Iteration
  log_likelihood_values <- c()
  for (i in 1:100) {
    gamma <- e_step(data, pi, mu, sigma)
    params <- m_step(gamma, data)
    pi <- params$pi
    mu <- params$mu
    sigma <- params$sigma
    log_likelihood_values[i] <- log_likelihood(data, pi, mu, sigma)
    # Check for convergence (log-likelihood stability)
    if (i > 1 && abs(log_likelihood_values[i] - log_likelihood_values[i-1]) < 1e-6) {
      break
    }
  }
  return(list(pi = pi, mu = mu, sigma = sigma, log_likelihood = log_likelihood_values))
}

# Function for Variational Inference (VI) using the EM algorithm approach
VI_GMM <- function(data, K) {
  # Using 'mixtools' package for a quick implementation of GMM with variational inference
  model <- normalmixEM(data, k = K)
  return(model)
}

# Run EM Algorithm on the data
em_result <- EM_GMM(data, K)

# Run Variational Inference (VI) on the data
vi_result <- VI_GMM(data, K)

# Print the results and compare with ground truth
cat("Ground Truth:\n")
cat("Mixing Coefficients (pi):", pi_true, "\n")
cat("Means (mu):", mu_true, "\n")
cat("Standard Deviations (sigma):", sigma_true, "\n\n")

cat("EM Algorithm Results:\n")
cat("Mixing Coefficients (pi):", em_result$pi, "\n")
cat("Means (mu):", em_result$mu, "\n")
cat("Standard Deviations (sigma):", em_result$sigma, "\n")

cat("\nVariational Inference Results:\n")
cat("Mixing Coefficients (pi):", vi_result$lambda, "\n")
cat("Means (mu):", vi_result$mu, "\n")
cat("Standard Deviations (sigma):", vi_result$sigma, "\n")

# Evaluate which method is closer to the ground truth by comparing the absolute differences
em_pi_diff <- abs(em_result$pi - pi_true)
vi_pi_diff <- abs(vi_result$lambda - pi_true)

em_mu_diff <- abs(em_result$mu - mu_true)
vi_mu_diff <- abs(vi_result$mu - mu_true)

em_sigma_diff <- abs(em_result$sigma - sigma_true)
vi_sigma_diff <- abs(vi_result$sigma - sigma_true)

cat("\nComparison of Differences from Ground Truth:\n")
cat("EM Algorithm - Pi Differences:", em_pi_diff, "\n")
cat("VI Algorithm - Pi Differences:", vi_pi_diff, "\n")

cat("EM Algorithm - Mu Differences:", em_mu_diff, "\n")
cat("VI Algorithm - Mu Differences:", vi_mu_diff, "\n")

cat("EM Algorithm - Sigma Differences:", em_sigma_diff, "\n")
cat("VI Algorithm - Sigma Differences:", vi_sigma_diff, "\n")