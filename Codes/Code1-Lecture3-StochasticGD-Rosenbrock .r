 # Load required packages
library(parallel)

# Define the Rosenbrock function and its derivative
f_rosenbrock <- function(x) {
  return(100 * (x[2] - x[1]^2)^2 + (1 - x[1])^2)  # Rosenbrock function
}

df_rosenbrock <- function(x) {
  df1 <- -400 * x[1] * (x[2] - x[1]^2) - 2 * (1 - x[1])
  df2 <- 200 * (x[2] - x[1]^2)
  return(c(df1, df2))  # Gradient of Rosenbrock function
}

# Stochastic Gradient Descent with momentum for Rosenbrock in parallel
stochastic_gradient_descent_rosenbrock_parallel <- function(epochs = 10000, lr = 0.001, momentum = 0.9, restarts = 500, x_range = c(-2, 2)) {
  # Detect the number of cores
  num_cores <- detectCores(logical = TRUE) - 1  # Leave one core free
  cl <- makeCluster(num_cores)
  # Ensure functions and variables are available to workers
  clusterExport(cl, varlist = c("f_rosenbrock", "df_rosenbrock", "lr", "momentum", "epochs", "x_range"), envir = environment())

  # Run SGD in parallel
  results <- parLapply(cl, 1:restarts, function(i) {
    set.seed(i)  # Ensure reproducibility per worker
    x <- runif(2, x_range[1], x_range[2])  # Random initialization in 2D
    v <- c(0, 0)  # Momentum term
    best_x_local <- x
    best_f_local <- f_rosenbrock(x)
    for (j in 1:epochs) {
      grad <- df_rosenbrock(x) + rnorm(2, mean = 0, sd = 0.05)  # Add small noise
      if (any(is.nan(grad)) || any(is.infinite(grad))) break  # Handle invalid gradients
      v <- momentum * v - lr * grad  # Update velocity
      x <- x + v  # Update position
      x <- pmin(pmax(x, x_range[1]), x_range[2])  # Clamp within range
      # Track best solution found
      current_f <- f_rosenbrock(x)
      if (current_f < best_f_local) {
        best_x_local <- x
        best_f_local <- current_f
      }
    }
    return(list(best_x = best_x_local, best_f = best_f_local))
  })

  # Stop parallel execution
  stopCluster(cl)

  # Find the overall best result
  best_x <- NULL
  best_f <- Inf
  for (result in results) {
    if (!is.null(result$best_x) && result$best_f < best_f) {
      best_x <- result$best_x
      best_f <- result$best_f
    }
  }

  return(list(best_x = best_x, best_f = best_f))
}

# Run the parallelized SGD
result_rosenbrock_parallel <- stochastic_gradient_descent_rosenbrock_parallel()

# Check if result is valid before printing
if (!is.null(result_rosenbrock_parallel$best_x)) {
  cat(sprintf("Best solution found: x = (%.5f, %.5f), f(x) = %.5f\n", 
              result_rosenbrock_parallel$best_x[1], 
              result_rosenbrock_parallel$best_x[2], 
              result_rosenbrock_parallel$best_f))
} else {
  cat("SGD failed to find a valid solution.\n")
}

# Plot the function
library(rgl)
x_vals <- seq(-2, 2, length.out = 100)
y_vals <- seq(-2, 2, length.out = 100)
z_vals <- outer(x_vals, y_vals, Vectorize(function(x, y) f_rosenbrock(c(x, y))))

persp3d(x_vals, y_vals, z_vals, col = "blue", alpha = 0.6, main = "SGD Optimization of Rosenbrock Function")
if (!is.null(result_rosenbrock_parallel$best_x)) {
  points3d(result_rosenbrock_parallel$best_x[1], result_rosenbrock_parallel$best_x[2], result_rosenbrock_parallel$best_f, col = "red", size = 10)
}