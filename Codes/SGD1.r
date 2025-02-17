# Define function and gradient
L <- function(x) { x^4 - x^2 }       # Function to minimize
grad_L <- function(x) { -2*x + 4*x^3 }  # Gradient of L(x)

# Gradient Descent Parameters
alpha <- 0.01      # Learning rate
tol <- 1e-6       # Tolerance for stopping criterion
max_iter <- 100000  # Maximum number of iterations
x <- 5            # Initial guess

# Gradient Descent Loop
for (i in 1:max_iter) {
  grad <- grad_L(x)  # Compute gradient
  x_new <- x - alpha * grad  # Update x
  # Check for convergence
  if (abs(x_new - x) < tol) {
    cat("Converged at iteration:", i, "\n")
    break
  }
  x <- x_new  # Update x
}

# Print results
cat("Minimum x:", x, "\n")
cat("Minimum L(x):", L(x), "\n")