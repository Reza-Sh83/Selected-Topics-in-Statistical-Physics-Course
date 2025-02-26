# Load necessary library
library(GA)

# Define the Rosenbrock function (to be minimized)
rosenbrock <- function(x) {
  return(100 * (x[2] - x[1]^2)^2 + (1 - x[1])^2)
}

# Run the Genetic Algorithm (without specifying crossover)
ga_result <- ga(type = "real-valued", 
                fitness = function(x) -rosenbrock(x),  # GA maximizes, so negate function
                lower = c(-2, -2), 
                upper = c(2, 2), 
                popSize = 50,      # Population size
                maxiter = 20000,     # Number of generations
                pmutation = 0.2)   # Mutation probability

# Print best solution found
best_solution <- ga_result@solution
best_value <- -ga_result@fitnessValue  # Convert back to minimization
cat(sprintf("Best solution found: x = (%.5f, %.5f), f(x) = %.5f\n", 
            best_solution[1], best_solution[2], best_value))

# Plot optimization progress
plot(ga_result)