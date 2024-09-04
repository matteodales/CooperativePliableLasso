




calc_mse <- function(actual, predicted) {
  return(mean((actual - predicted)^2))
}

# this function, given a source X, modifying variables Z and interaction coefficients theta,
# computes the interaction terms of the pliable model




compute_pliable<-function(X, Z, theta){
  p=ncol(X)
  N=nrow(X)
  K=ncol(Z)
  
  xz_theta <- lapply(seq_len(p),
                     function(j) (matrix(X[, j], nrow = N, ncol = K) * Z) %*% t(theta)[, j])
  xz_term<- (Reduce(f = '+', x = xz_theta))
  
  return(xz_term)
  
  
}





coop_pliable <- function(data_x1, data_x2, data_z, response_y, coop_alpha=0, fold_indices, num_folds=10, iterations=3, tt_param=NULL, penalty_factors=NA, linear_z=FALSE, pliable_alpha=0.5) {
  
  # Combine data_x1 and data_x2 for the initial setup
  combined_data_x = rbind(cbind(data_x1, data_x2),
                          cbind(-sqrt(coop_alpha) * data_x1, sqrt(coop_alpha) * data_x2))
  
  # Expand the response vector with zeros for the augmented part
  augmented_y = c(response_y, rep(0, nrow(data_x1)))
  
  # Combine data_z with a matrix of zeros for the augmented part
  combined_z = rbind(data_z, matrix(0, nrow(data_x1), ncol(data_z)))
  
  # Set penalty factors if not provided
  if (length(penalty_factors) == 1) {
    penalty_factors = rep(1, ncol(combined_data_x))
  }
  
  # Fit the initial pliable lasso model to obtain the lambda sequence
  pliable_fit = pliable(combined_data_x, combined_z, augmented_y, nlambda = 50, zlinear = linear_z, tt = tt_param, penalty.factor = penalty_factors, alpha = pliable_alpha)
  lambda_seq = pliable_fit$lambda
  
  # Initialize storage for models and errors
  fold_models = vector("list", num_folds)
  error_matrix = matrix(NA, num_folds, length(lambda_seq))
  
  # Loop through each fold
  for (fold in seq_len(num_folds)) {
    test_indices = (fold_indices == fold)
    train_indices = !test_indices
    
    # Subset the data for training
    x1_train = data_x1[train_indices, , drop = FALSE]
    x2_train = data_x2[train_indices, , drop = FALSE]
    z_train = data_z[train_indices, , drop = FALSE]
    y_train = response_y[train_indices]
    
    # Center the training data (no scaling)
    x1_train = scale(x1_train, center = TRUE, scale = FALSE)
    x2_train = scale(x2_train, center = TRUE, scale = FALSE)
    z_train = scale(z_train, center = TRUE, scale = FALSE)
    
    # Prepare the augmented data for pliable lasso
    combined_x_train = rbind(cbind(x1_train, x2_train),
                             cbind(-sqrt(coop_alpha) * x1_train, sqrt(coop_alpha) * x2_train))
    combined_z_train = rbind(z_train, matrix(0, nrow(x1_train), ncol(z_train)))
    augmented_y_train = c(y_train, rep(0, nrow(x1_train)))
    
    # Fit the pliable lasso model on the training data
    pliable_model = pliable(combined_x_train, combined_z_train, augmented_y_train, lambda = lambda_seq, zlinear = linear_z, tt = tt_param, penalty.factor = penalty_factors, alpha = pliable_alpha)
    fold_models[[fold]] = pliable_model
    
    # Predict on the test set
    x_test = cbind(data_x1[test_indices, , drop = FALSE], data_x2[test_indices, , drop = FALSE])
    z_test = data_z[test_indices, , drop = FALSE]
    predictions = predict(pliable_model, x_test, z_test, lambda = lambda_seq)
    
    # Calculate the mean squared error for the test set
    true_y = response_y[test_indices]
    err_mse = (predictions - replicate(ncol(predictions), true_y))^2
    err_cv = apply(err_mse, 2, mean, na.rm = TRUE)
    error_matrix[fold, ] = err_cv

  }
  
  # Calculate the cross-validated mean and standard deviation of the errors
  cv_mean = colMeans(error_matrix, na.rm = TRUE)
  cv_sd = sqrt(colMeans(scale(error_matrix, center = cv_mean, scale = FALSE)^2, na.rm = TRUE) / (num_folds - 1))
  
  # Find the index of the minimum mean error
  min_error_index = which.min(cv_mean)
  
  # Return the results
  return(list(cvm = cv_mean, cvsd = cv_sd, lambda = lambda_seq, lambda_min = lambda_seq[min_error_index], best_fit = pliable_fit, min_index = min_error_index))
}