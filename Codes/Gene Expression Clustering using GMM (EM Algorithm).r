set.seed(42)
library(ggplot2)

# Generate synthetic gene expression data
generate_gene_expression <- function(n_genes, n_samples, means, sds) {
  data <- matrix(NA, nrow=n_samples, ncol=n_genes)
  condition_labels <- character(n_samples)  # Vector to store condition labels for all samples
  for (i in 1:n_genes) {
    # Randomly assign each sample to one of the three conditions
    condition_labels <- sample(1:3, n_samples, replace = TRUE)
    # Generate gene expression data for the gene based on the condition
    for (j in 1:n_samples) {
      data[j, i] <- rnorm(1, mean=means[condition_labels[j]], sd=sds[condition_labels[j]])
    }
  }
  return(list(data=data, conditions=condition_labels))
}

# Parameters for gene expression data
n_genes <- 100  # Number of genes
n_samples <- 1000  # Number of samples
means <- c(10, 5, 2)  # Mean gene expression levels for the three conditions
sds <- c(2, 1.5, 2.5)  # Standard deviations for the three conditions

# Generate gene expression data
gene_expression_results <- generate_gene_expression(n_genes, n_samples, means, sds)

# Create a data frame for easy visualization
gene_expression_data <- gene_expression_results$data
condition_labels <- gene_expression_results$conditions

gene_expression_df <- data.frame(
  gene = rep(1:n_genes, each = n_samples),
  sample = rep(1:n_samples, n_genes),
  expression = c(gene_expression_data),
  condition = rep(condition_labels, n_genes)
)

# GMM parameters and functions (EM Algorithm)
initialize_params <- function(k) {
  mu <- rnorm(k, mean=0, sd=5)  # Initial means
  sigma2 <- rep(1, k)  # Initial variances
  pi <- rep(1/k, k)  # Initial mixture weights
  return(list(mu=mu, sigma2=sigma2, pi=pi))
}

# Gaussian PDF
gaussian_pdf <- function(x, mu, sigma2) {
  return(dnorm(x, mean=mu, sd=sqrt(sigma2)))
}

# Log likelihood
log_likelihood <- function(x, phi, pi) {
  return(sum(log(rowSums(phi))))
}

# E-Step
e_step <- function(x, mu, sigma2, pi, k, n) {
  phi <- matrix(0, nrow=n, ncol=k)
  for (i in 1:k) {
    phi[, i] <- pi[i] * gaussian_pdf(x, mu[i], sigma2[i])
  }
  row_sums <- rowSums(phi)
  phi <- phi / row_sums
  return(phi)
}

# M-Step
m_step <- function(x, phi, k, n) {
  N_k <- colSums(phi)
  mu <- colSums(phi * x) / N_k
  sigma2 <- colSums(phi * (x - mu[apply(phi, 1, which.max)])^2) / N_k
  pi <- N_k / n
  return(list(mu=mu, sigma2=sigma2, pi=pi))
}

# EM Algorithm
em_algorithm <- function(x, k, max_iter) {
  n <- length(x)
  params <- initialize_params(k)
  mu <- params$mu
  sigma2 <- params$sigma2
  pi <- params$pi
  log_likelihoods <- numeric(max_iter)
  for (iter in 1:max_iter) {
    phi <- e_step(x, mu, sigma2, pi, k, n)
    params <- m_step(x, phi, k, n)
    mu <- params$mu
    sigma2 <- params$sigma2
    pi <- params$pi
    log_likelihoods[iter] <- log_likelihood(x, phi, pi)
    if (iter %% 10 == 0) {
      cat("Iteration:", iter, "Log Likelihood:", log_likelihoods[iter], "\n")
    }
  }
  return(list(mu=mu, sigma2=sigma2, pi=pi, phi=phi, log_likelihoods=log_likelihoods))
}

# Prepare the data for EM (flatten the gene expression data into a single vector)
gene_expression_flattened <- as.vector(gene_expression_df$expression)

# Run EM Algorithm on gene expression data
k <- 3  # Number of clusters (for example, the three conditions)
max_iter <- 100
result <- em_algorithm(gene_expression_flattened, k, max_iter)

# Assign clusters based on highest responsibility
cluster_assignments <- apply(result$phi, 1, which.max)

# Add cluster assignments to the data frame
gene_expression_df$cluster <- as.factor(cluster_assignments)

# **Improved Histogram Plot with Transparency**
ggplot(gene_expression_df, aes(x=expression, fill=cluster)) +
  geom_histogram(bins=51, alpha=0.4, position="identity") +  # Increased transparency for better overlap
  geom_vline(xintercept=result$mu, linetype="dashed", color="red", size=1) +
  scale_fill_manual(values=c("#FF9999", "#99CCFF", "#99FF99")) +  # Custom colors for clarity
  labs(title="Gene Expression Clustering using GMM (EM Algorithm)", 
       x="Gene Expression", y="Frequency") +
  theme_minimal()

# **Alternative Density Plot for Smooth Visualization**
ggplot(gene_expression_df, aes(x=expression, fill=cluster, color=cluster)) +
  geom_density(alpha=0.4) +  # Density plots with transparency
  geom_vline(xintercept=result$mu, linetype="dashed", color="red", size=1) +
  scale_fill_manual(values=c("#FF9999", "#99CCFF", "#99FF99")) +
  scale_color_manual(values=c("#FF4444", "#4466FF", "#44FF44")) +
  labs(title="Gene Expression Clustering using GMM (EM Algorithm)", 
       x="Gene Expression", y="Density") +
  theme_minimal()