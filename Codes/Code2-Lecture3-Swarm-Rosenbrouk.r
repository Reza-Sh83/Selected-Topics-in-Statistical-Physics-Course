# Load necessary library
library(ggplot2)

# Define the Rosenbrock function
rosenbrock_function <- function(x, y) {
  return((1 - x)^2 + 100 * (y - x^2)^2)
}

# PSO parameters
num_particles <- 2000
num_dimensions <- 2
max_iters <- 100
w <- 0.7
c1 <- 1.5
c2 <- 1.5

# Initialize particle positions and velocities randomly
set.seed(42)
positions <- matrix(runif(num_particles * num_dimensions, min = -2, max = 2), 
                    nrow = num_particles, ncol = num_dimensions)
velocities <- matrix(runif(num_particles * num_dimensions, min = -1, max = 1), 
                     nrow = num_particles, ncol = num_dimensions)

# Personal best positions and values
p_best_positions <- positions
p_best_values <- apply(p_best_positions, 1, function(p) rosenbrock_function(p[1], p[2]))

# Global best position and value
g_best_index <- which.min(p_best_values)
g_best_position <- p_best_positions[g_best_index, ]
g_best_value <- p_best_values[g_best_index]

# Store data for static visualization
pso_data <- list()

# PSO loop
for (iter in 1:max_iters) {
  # Generate random numbers
  r1 <- matrix(runif(num_particles * num_dimensions), nrow = num_particles, ncol = num_dimensions)
  r2 <- matrix(runif(num_particles * num_dimensions), nrow = num_particles, ncol = num_dimensions)
  # Update velocities
  velocities <- w * velocities + 
                c1 * r1 * (p_best_positions - positions) + 
                c2 * r2 * (matrix(rep(g_best_position, num_particles), nrow = num_particles, byrow = TRUE) - positions)
  # Update positions
  positions <- positions + velocities
  # Evaluate new positions
  new_values <- apply(positions, 1, function(p) rosenbrock_function(p[1], p[2]))
  # Update personal bests
  improved <- new_values < p_best_values
  p_best_positions[improved, ] <- positions[improved, ]
  p_best_values[improved] <- new_values[improved]
  # Update global best
  g_best_index <- which.min(p_best_values)
  g_best_position <- p_best_positions[g_best_index, ]
  g_best_value <- p_best_values[g_best_index]
  # Store every 10 iterations
  if (iter %% 10 == 0) {
    pso_data[[as.character(iter)]] <- data.frame(
      iter = iter,
      x = positions[, 1], 
      y = positions[, 2], 
      vx = velocities[, 1], 
      vy = velocities[, 2]
    )
  }
}

# Calculate mean values of x and y at the last step
mean_x <- mean(positions[, 1])
mean_y <- mean(positions[, 2])

# Print the mean values
cat("Mean of x at last step:", mean_x, "\n")
cat("Mean of y at last step:", mean_y, "\n")

# Combine data for plotting
pso_df <- do.call(rbind, pso_data)

# Plot particle positions and velocity vectors
ggplot(pso_df, aes(x, y)) +
  geom_point(aes(color = factor(iter)), size = 3, alpha = 0.7) +  # Particles
  geom_segment(aes(xend = x + vx * 0.1, yend = y + vy * 0.1), 
               arrow = arrow(type = "closed", length = unit(0.1, "inches")), 
               alpha = 0.5) +  # Velocity vectors
  facet_wrap(~ iter, ncol = 3) +  # Show different iterations
  labs(title = "Particle Swarm Optimization - Rosenbrock Function", 
       x = "X Position", y = "Y Position", color = "Iteration") +
  theme_minimal()