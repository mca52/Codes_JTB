# Load the inputs from TRANSFIL model - costs, DALYs averted, prob of elimination

# Load necessary libraries
library(lpSolve)
library(MASS)
library(Matrix)
library(caret)

# Function to compute the Linear Optimal Transport (LOT)
compute_lot <- function(mu_i, sigma) {
  # Placeholder for the actual transport map computation
  # For simplicity, using mu_i directly
  T_star <- mu_i
  projection <- T_star - diag(nrow(T_star))
  return(projection)
}

# Function to perform PCA
perform_pca <- function(P) {
  P_centered <- scale(P, center = TRUE, scale = FALSE)
  covariance_matrix <- cov(P_centered)
  if (nrow(covariance_matrix) == 0 || ncol(covariance_matrix) == 0) {
    stop("Covariance matrix is empty or incorrectly computed.")
  }
  eigen_decomp <- eigen(covariance_matrix)
  pca_projection_matrix <- eigen_decomp$vectors[, 1:2] # Top 2 components
  projected_data <- P_centered %*% pca_projection_matrix
  return(list(projected_data = projected_data, pca_matrix = pca_projection_matrix))
}

# Function to train LDA classifier and evaluate
train_evaluate_lda <- function(X_train, y_train, X_test, y_test) {
  lda_model <- lda(X_train, grouping = y_train)
  predictions <- predict(lda_model, X_test)$class
  accuracy <- mean(predictions == y_test)
  return(accuracy)
}

# Function to compute the utility function
compute_utility <- function(predictions, true_labels) {
  utility <- mean(predictions == true_labels)
  return(utility)
}

# Function to simulate the dataset
simulate_dataset <- function(N, n_features) {
  costs <- matrix(rnorm(N * n_features), nrow = N, ncol = n_features)
  dalys <- matrix(rnorm(N * n_features), nrow = N, ncol = n_features)
  labels <- sample(c(0, 1), N, replace = TRUE)
  return(list(costs = costs, dalys = dalys, labels = labels))
}

# Main function to compute EVSI
compute_evsi <- function(N, M, p) {
  
  # Step 1: Compute LOT for both costs and DALYs
  projections_costs <- t(apply(costs, 1, function(mu_i) compute_lot(mu_i, costs)))
  projections_dalys <- t(apply(dalys, 1, function(mu_i) compute_lot(mu_i, dalys)))
  
  # Combine projections for PCA
  projections <- cbind(projections_costs, projections_dalys)
  
  # Step 3: PCA on projections
  pca_results <- perform_pca(projections)
  X <- pca_results$projected_data
  
  # Split dataset into training (80%) and testing (20%)
  set.seed(123)
  train_indices <- createDataPartition(true_labels, p = 0.8, list = FALSE)
  X_train <- X[train_indices, ]
  y_train <- true_labels[train_indices]
  X_test <- X[-train_indices, ]
  y_test <- true_labels[-train_indices]
  
  # Train LDA classifier
  EU_current <- train_evaluate_lda(X_train, y_train, X_test, y_test)
  
  # Step 5: Simulate and augment dataset
  EU_new_list <- numeric(M)
  for (k in 1:M) {
    data_sim <- simulate_dataset(N, n_features)
    costs_sim <- data_sim$costs
    dalys_sim <- data_sim$dalys
    true_labels_sim <- data_sim$labels
    
    # Compute LOT for simulated costs and DALYs
    projections_costs_sim <- t(apply(costs_sim, 1, function(mu_i) compute_lot(mu_i, costs_sim)))
    projections_dalys_sim <- t(apply(dalys_sim, 1, function(mu_i) compute_lot(mu_i, dalys_sim)))
    
    # Combine projections for PCA
    projections_sim <- cbind(projections_costs_sim, projections_dalys_sim)
    
    # Augment training data
    augmented_projections <- rbind(projections, projections_sim)
    augmented_labels <- c(true_labels, true_labels_sim)
    pca_results_aug <- perform_pca(augmented_projections)
    X_train_aug <- pca_results_aug$projected_data[1:N, ]
    y_train_aug <- augmented_labels[1:N]
    
    # Train new LDA model
    EU_new <- train_evaluate_lda(X_train_aug, y_train_aug, X_test, y_test)
    EU_new_list[k] <- EU_new
  }
  
  # Step 6: Compute EVSI
  EUsample <- mean(EU_new_list)
  EVSI <- EUsample - EU_current
  
  return(EVSI)
}

# Example usage
N <- 1000 # Number of data points
M <- 10  # Number of simulated datasets
p <- 10   # Number of PCA components
evsi <- compute_evsi(N, M, p)
print(paste("Estimated EVSI:", evsi))
