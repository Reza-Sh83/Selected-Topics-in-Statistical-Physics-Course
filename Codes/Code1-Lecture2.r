# Stochastic Gradient Descent in R to find the global minimum of f(x) = x^3 + sin(10*x)

# Define the function and its derivative
f <- function(x) {
  return(x^3 + sin(10 * x))
}

df <- function(x) {
  return(3 * x^2 + 10 * cos(10 * x))
}

# Stochastic Gradient Descent with momentum
stochastic_gradient_descent <- function(epochs = 1000, lr = 0.005, momentum = 0.9, restarts = 5, x_range = c(-1, 0)) {
  best_x <- NULL
  best_f <- Inf
  for (i in 1:restarts) {
    x <- runif(1, x_range[1], x_range[2])  # Random initialization
    v <- 0  # Momentum term
    for (j in 1:epochs) {
      grad <- df(x) + rnorm(1, mean = 0, sd = 0.05)  # Add small noise for stochasticity
      if (is.nan(grad) || is.infinite(grad)) break  # Handle invalid gradients
      v <- momentum * v - lr * grad  # Update velocity
      x <- x + v  # Update position
      # Clamp x within [-1, 0] to avoid instability
      x <- max(min(x, x_range[2]), x_range[1])
      # Track best solution found
      if (!is.nan(f(x)) && f(x) < best_f) {
        best_x <- x
        best_f <- f(x)
      }
    }
  }
  return(list(best_x = best_x, best_f = best_f))
}

# Run SGD to find the global minimum
result <- stochastic_gradient_descent()
cat(sprintf("Best solution found: x = %.5f, f(x) = %.5f\n", result$best_x, result$best_f))

# Plot function and found minimum
x_vals <- seq(-1, 0, length.out = 1000)
y_vals <- sapply(x_vals, f)
plot(x_vals, y_vals, type = "l", col = "blue", main = "SGD Optimization of f(x)", xlab = "x", ylab = "f(x)")
points(result$best_x, result$best_f, col = "red", pch = 19)
legend("topright", legend = c("f(x)", "SGD Minimum"), col = c("blue", "red"), lty = 1, pch = c(NA, 19))