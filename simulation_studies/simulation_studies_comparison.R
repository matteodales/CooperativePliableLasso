


# This script performs the analysis described in the simulation studies section.
# By adjusting parameters, we can explore different settings.
# Each simulation is repeated `num_simulations` times. We study test MSE, as well as sensitivity and specificity.

# Load the necessary library
library(pliable)
source('functions.R')

# Define rho values to explore
rho_values = c(0, 0.2, 0.4, 0.6, 0.8, 1, 3, 5, 9)

# Simulation parameters
simulation_settings <- list(
  list(sigma = 33, factor_strength = 4, tx1 = 2, tx2 = 2),
  list(sigma = 40, factor_strength = 4, tx1 = 4, tx2 = 4),
  list(sigma = 22, factor_strength = 4, tx1 = 0, tx2 = 0),
  list(sigma = 24, factor_strength = 8, tx1 = 0, tx2 = 0)
)



# Loop over simulations
for(sim_idx in 1:4){
  
  
  # Simulation setup
  seed = 100
  num_simulations = 10
  num_folds = 5
  
  # Extract parameters for this simulation
  sigma_val = simulation_settings[[sim_idx]]$sigma
  factor_strength = simulation_settings[[sim_idx]]$factor_strength
  tx1_effect = simulation_settings[[sim_idx]]$tx1
  tx2_effect = simulation_settings[[sim_idx]]$tx2
  
  
  # Define additional parameters
  n_train = 200
  n_test = 9800
  p1 = 500
  p2 = 500
  signal_predictors = 30
  u_std_dev = 1
  tz_effect = 1
  num_factors = 4
  snr_avg = 0
  
  # Initialize storage for results
  mse_x1 = numeric(num_simulations)
  mse_x2 = numeric(num_simulations)
  mse_early = numeric(num_simulations)
  mse_late = numeric(num_simulations)
  mse_coop = numeric(num_simulations)
  mse_coop_adap = numeric(num_simulations)
  
  beta_count_x1 = numeric(num_simulations)
  beta_count_x2 = numeric(num_simulations)
  beta_count_early = numeric(num_simulations)
  beta_count_late = numeric(num_simulations)
  beta_count_coop = numeric(num_simulations)
  beta_count_coop_adap = numeric(num_simulations)
  
  theta_count_x1 = numeric(num_simulations)
  theta_count_x2 = numeric(num_simulations)
  theta_count_early = numeric(num_simulations)
  theta_count_late = numeric(num_simulations)
  theta_count_coop = numeric(num_simulations)
  theta_count_coop_adap = numeric(num_simulations)
  
  beta_presence_x1 = numeric(p1 + p2)
  beta_presence_x2 = numeric(p1 + p2)
  beta_presence_early = numeric(p1 + p2)
  beta_presence_late = numeric(p1 + p2)
  beta_presence_coop = numeric(p1 + p2)
  beta_presence_coop_adap = numeric(p1 + p2)
  
  selected_rhos = numeric(num_simulations)
  selected_rhos_adap = numeric(num_simulations)
  
  
  # Define the real beta and theta coefficients
  beta_U = c(rep(factor_strength, signal_predictors))
  
  beta1_real = rep(0, p1)
  beta2_real = rep(0, p2)
  beta1_real[1:signal_predictors] = factor_strength
  if (sim_idx != 4) beta2_real[1:signal_predictors] = factor_strength
  if (sim_idx == 3) beta2_real[1:signal_predictors] = -factor_strength
  real_main_effects = c(beta1_real, beta2_real) != 0
  
  theta1_real = Matrix(0, p1, num_factors, sparse = TRUE)
  theta1_real[3, 1] = factor_strength
  theta1_real[4, 2] = -factor_strength
  theta1_real[1, 3] = factor_strength
  
  theta2_real = Matrix(0, p2, num_factors, sparse = TRUE)
  if (sim_idx != 4) {
    theta2_real[3, 1] = -factor_strength
    theta2_real[4, 2] = factor_strength
    theta2_real[1, 3] = -factor_strength
  }
  
  real_interactions = rbind(theta1_real, theta2_real) != 0

  # Run simulations
  for (j in seq_len(num_simulations)) {
    
    print(j)
    
    set.seed(j + seed)
    
    # Generate training data
    X1_train = matrix(rnorm(n_train * p1), nrow = n_train, ncol = p1)
    X2_train = matrix(rnorm(n_train * p2), nrow = n_train, ncol = p2)
    Z_train = matrix(rnorm(n_train * num_factors), nrow = n_train, ncol = num_factors)
    U_train = matrix(0, n_train, signal_predictors)
    
    for (i in seq_len(signal_predictors)) {
      u_sample = rnorm(n_train, sd = u_std_dev)
      X1_train[, i] = X1_train[, i] + tx1_effect * u_sample
      X2_train[, i] = X2_train[, i] + tx2_effect * u_sample
      U_train[, i] = U_train[, i] + u_sample
      if (i <= num_factors) {
        Z_train[, i] = Z_train[, i] + tz_effect * u_sample
      }
    }
    
    X1_train = scale(X1_train, center = TRUE, scale = TRUE)
    X2_train = scale(X2_train, center = TRUE, scale = TRUE)
    Z_train = scale(Z_train, center = TRUE, scale = TRUE)
    
    # Compute pliable components
    pliable1_train = compute_pliable(X1_train, Z_train, theta1_real)
    pliable2_train = compute_pliable(X2_train, Z_train, theta2_real)
    
    # Generate response variable
    mu_train = X1_train %*% beta1_real + pliable1_train + X2_train %*% beta2_real + pliable2_train
    y_train = mu_train + sigma_val * rnorm(n_train)
    y_train = y_train - mean(y_train)
    
    fold_ids = sample(rep_len(1:num_folds, nrow(X1_train)))
    
    snr = var(mu_train) / var(y_train - mu_train)
    snr_avg = snr_avg + snr / num_simulations
    cat("", fill = TRUE)
    cat(c("SNR =", snr), fill = TRUE)
    cat("", fill = TRUE)
    
    # Generate test data
    set.seed(j + seed)
    
    X1_test = matrix(rnorm(n_test * p1), nrow = n_test, ncol = p1)
    X2_test = matrix(rnorm(n_test * p2), nrow = n_test, ncol = p2)
    Z_test = matrix(rnorm(n_test * num_factors), nrow = n_test, ncol = num_factors)
    U_test = matrix(0, n_test, signal_predictors)
    
    for (i in seq_len(signal_predictors)) {
      u_sample = rnorm(n_test, sd = u_std_dev)
      X1_test[, i] = X1_test[, i] + tx1_effect * u_sample
      X2_test[, i] = X2_test[, i] + tx2_effect * u_sample
      U_test[, i] = U_test[, i] + u_sample
      if (i <= num_factors) {
        Z_test[, i] = Z_test[, i] + tz_effect * u_sample
      }
    }
    
    X1_test = scale(X1_test, center = TRUE, scale = TRUE)
    X2_test = scale(X2_test, center = TRUE, scale = TRUE)
    Z_test = scale(Z_test, center = TRUE, scale = TRUE)
    
    pliable1_test = compute_pliable(X1_test, Z_test, theta1_real)
    pliable2_test = compute_pliable(X2_test, Z_test, theta2_real)
    
    mu_test = X1_test %*% beta1_real + pliable1_test + X2_test %*% beta2_real + pliable2_test
    y_test = mu_test + sigma_val * rnorm(n_test)
    y_test = y_test - mean(y_test)
    
    
    # Model: Only X1
    print('Only X1')
    fit_x1 = pliable(x = X1_train, z = Z_train, y = y_train, nlambda = 50, zlinear = FALSE)
    cv_x1 = cv.pliable(fit_x1, x = X1_train, z = Z_train, y = y_train, foldid = fold_ids, verbose = FALSE)
    pred_x1 = predict(fit_x1, x = X1_test, z = Z_test, lambda = cv_x1$lambda.min)
    mse_x1[j] = calc_mse(mu_test, pred_x1)
    beta_count_x1[j] = fit_x1$nbeta[which(cv_x1$lambda == cv_x1$lambda.min)]
    theta_count_x1[j] = fit_x1$ntheta[which(cv_x1$lambda == cv_x1$lambda.min)]
    beta_presence_x1 = beta_presence_x1 + (c(fit_x1$beta[, which(cv_x1$lambda == cv_x1$lambda.min)], rep(0, p2)) != 0)
    lambda1 = cv_x1$lambda.min
    print(mse_x1[j])
    
    # Model: Only X2
    print('Only X2')
    fit_x2 = pliable(x = X2_train, z = Z_train, y = y_train, nlambda = 50, zlinear = FALSE)
    cv_x2 = cv.pliable(fit_x2, x = X2_train, z = Z_train, y = y_train, foldid = fold_ids, verbose = FALSE)
    pred_x2 = predict(fit_x2, x = X2_test, z = Z_test, lambda = cv_x2$lambda.min)
    mse_x2[j] = calc_mse(mu_test, pred_x2)
    beta_count_x2[j] = fit_x2$nbeta[which(cv_x2$lambda == cv_x2$lambda.min)]
    theta_count_x2[j] = fit_x2$ntheta[which(cv_x2$lambda == cv_x2$lambda.min)]
    beta_presence_x2 = beta_presence_x2 + (c(rep(0, p1), fit_x2$beta[, which(cv_x2$lambda == cv_x2$lambda.min)]) != 0)
    lambda2 = cv_x2$lambda.min
    print(mse_x2[j])
    
    # Model: Early Fusion
    print('Early fusion')
    X_train_combined = cbind(X1_train, X2_train)
    X_test_combined = cbind(X1_test, X2_test)
    fit_early = pliable(x = X_train_combined, z = Z_train, y = y_train, nlambda = 50, zlinear = FALSE)
    cv_early = cv.pliable(fit_early, x = X_train_combined, z = Z_train, y = y_train, foldid = fold_ids, verbose = FALSE)
    pred_early = predict(fit_early, x = X_test_combined, z = Z_test, lambda = cv_early$lambda.min)
    mse_early[j] = calc_mse(mu_test, pred_early)
    beta_count_early[j] = fit_early$nbeta[which(cv_early$lambda == cv_early$lambda.min)]
    theta_count_early[j] = fit_early$ntheta[which(cv_early$lambda == cv_early$lambda.min)]
    beta_presence_early = beta_presence_early + (fit_early$beta[, which(cv_early$lambda == cv_early$lambda.min)] != 0)
    print(mse_early[j])
    
    # Model: Late Fusion
    print('Late fusion')
    validation_fraction = 0.2
    num_val_samples = floor(validation_fraction * nrow(X1_train))
    val_indices = sort(sample(seq_len(nrow(X1_train)), size = num_val_samples))
    train_indices_late = setdiff(seq_len(nrow(X1_train)), val_indices)
    X1_val = X1_train[val_indices, ]
    X1_train_late = X1_train[train_indices_late, ]
    X2_val = X2_train[val_indices, ]
    X2_train_late = X2_train[train_indices_late, ]
    Z_val = Z_train[val_indices, ]
    Z_train_late = Z_train[train_indices_late, ]
    y_val = y_train[val_indices]
    y_train_late = y_train[train_indices_late]
    
    fit_x1_late = pliable(X1_train_late, Z_train_late, y_train_late, zlinear = FALSE)
    cv_x1_late = cv.pliable(fit_x1_late, X1_train_late, Z_train_late, y_train_late, nfolds = 5, verbose = FALSE)
    
    fit_x2_late = pliable(X2_train_late, Z_train_late, y_train_late, zlinear = FALSE)
    cv_x2_late = cv.pliable(fit_x2_late, X2_train_late, Z_train_late, y_train_late, nfolds = 5, verbose = FALSE)
    
    X1_pred_val = predict(fit_x1_late, x = X1_val, z = Z_val, lambda = cv_x1_late$lambda.min)
    X1_pred_test = predict(fit_x1_late, x = X1_test, z = Z_test, lambda = cv_x1_late$lambda.min)
    X2_pred_val = predict(fit_x2_late, x = X2_val, z = Z_val, lambda = cv_x2_late$lambda.min)
    X2_pred_test = predict(fit_x2_late, x = X2_test, z = Z_test, lambda = cv_x2_late$lambda.min)
    
    fused_data = data.frame(y = y_val, X1_pred = as.vector(X1_pred_val), X2_pred = as.vector(X2_pred_val))
    fit_fusion = lm(y ~ X1_pred + X2_pred, data = fused_data)
    pred_late = predict(fit_fusion, data.frame(X1_pred = as.vector(X1_pred_test), X2_pred = as.vector(X2_pred_test)))
    
    mse_late[j] = calc_mse(mu_test, pred_late)
    beta_count_late[j] = fit_x1_late$nbeta[which(cv_x1_late$lambda == cv_x1_late$lambda.min)] +
      fit_x2_late$nbeta[which(cv_x2_late$lambda == cv_x2_late$lambda.min)]
    theta_count_late[j] = fit_x1_late$ntheta[which(cv_x1_late$lambda == cv_x1_late$lambda.min)] +
      fit_x2_late$ntheta[which(cv_x2_late$lambda == cv_x2_late$lambda.min)]
    beta_presence_late = beta_presence_late + (c(fit_x1_late$beta[, which(cv_x1_late$lambda == cv_x1_late$lambda.min)],
                                                 fit_x2_late$beta[, which(cv_x2_late$lambda == cv_x2_late$lambda.min)]) != 0)
    
    print(mse_late[j])
    
    # Model: Cooperative Learning
    print('CoopLearn')
    cvm_values = numeric(length(rho_values))
    mse_test_no_pf = numeric(length(rho_values))
    beta_count_coop_rho = numeric(length(rho_values))
    theta_count_coop_rho = numeric(length(rho_values))
    beta_coop_rho = Matrix(0, p1 + p2, length(rho_values))
    
    for (rho_idx in seq_along(rho_values)) {
      current_rho = rho_values[rho_idx]
      print(current_rho)
      fit_coop_pliable = coop_pliable(X1_train, X2_train, Z_train, y_train,
                                           coop_alpha = current_rho, fold_indices = fold_ids,
                                           num_folds = max(fold_ids))
      
      cvm_values[rho_idx] = min(fit_coop_pliable$cvm)
      beta_count_coop_rho[rho_idx] = fit_coop_pliable$best_fit$nbeta[fit_coop_pliable$min_ind]
      theta_count_coop_rho[rho_idx] = fit_coop_pliable$best_fit$ntheta[fit_coop_pliable$min_ind]
      mse_test_no_pf[rho_idx] = calc_mse(mu_test, predict(fit_coop_pliable$best_fit, cbind(X1_test, X2_test), Z_test)[, fit_coop_pliable$min_ind])
      print(mse_test_no_pf[rho_idx])
      
      beta_coop_rho[, rho_idx] = fit_coop_pliable$best_fit$beta[, fit_coop_pliable$min_ind]
    }
    
    selected_rho = rho_values[which.min(cvm_values)]
    print(selected_rho)
    
    beta_count_coop[j] = beta_count_coop_rho[which.min(cvm_values)]
    theta_count_coop[j] = theta_count_coop_rho[which.min(cvm_values)]
    beta_coop_selected = beta_coop_rho[, which.min(cvm_values)]
    beta_presence_coop = beta_presence_coop + (beta_coop_selected != 0)
    mse_coop[j] = mse_test_no_pf[which.min(cvm_values)]
    selected_rhos[j] = selected_rho
    
    # Model: Adaptive Cooperative Learning
    print('AdapLearn')
    cvm_values_adap = numeric(length(rho_values))
    mse_test_no_pf_adap = numeric(length(rho_values))
    beta_count_coop_rho_adap = numeric(length(rho_values))
    theta_count_coop_rho_adap = numeric(length(rho_values))
    beta_coop_rho_adap = Matrix(0, p1 + p2, length(rho_values))
    
    for (rho_idx in seq_along(rho_values)) {
      current_rho = rho_values[rho_idx]
      print(current_rho)
      
      fit_coop_pliable_adap = coop_pliable(X1_train, X2_train, Z_train, y_train,
                                           coop_alpha = current_rho, fold_indices = fold_ids,
                                           num_folds = max(fold_ids),
                                           penalty_factors = c(rep(1, p1), rep(lambda2 / lambda1, p2)))
      
      cvm_values_adap[rho_idx] = min(fit_coop_pliable_adap$cvm)
      beta_count_coop_rho_adap[rho_idx] = fit_coop_pliable_adap$best_fit$nbeta[fit_coop_pliable_adap$min_ind]
      theta_count_coop_rho_adap[rho_idx] = fit_coop_pliable_adap$best_fit$ntheta[fit_coop_pliable_adap$min_ind]
      mse_test_no_pf_adap[rho_idx] = calc_mse(mu_test, predict(fit_coop_pliable_adap$best_fit, cbind(X1_test, X2_test), Z_test)[, fit_coop_pliable_adap$min_ind])
      print(mse_test_no_pf_adap[rho_idx])
      
      beta_coop_rho_adap[, rho_idx] = fit_coop_pliable_adap$best_fit$beta[, fit_coop_pliable_adap$min_ind]
    }
    
    selected_rho_adap = rho_values[which.min(cvm_values_adap)]
    print(selected_rho_adap)
    
    beta_count_coop_adap[j] = beta_count_coop_rho_adap[which.min(cvm_values_adap)]
    theta_count_coop_adap[j] = theta_count_coop_rho_adap[which.min(cvm_values_adap)]
    beta_coop_adap_selected = beta_coop_rho_adap[, which.min(cvm_values_adap)]
    beta_presence_coop_adap = beta_presence_coop_adap + (beta_coop_adap_selected != 0)
    mse_coop_adap[j] = mse_test_no_pf_adap[which.min(cvm_values_adap)]
    selected_rhos_adap[j] = selected_rho_adap
  }
  
  # Save results to file
  sim_filename = paste0(seed, "_ntrain", n_train, "_ntest", n_test,
                        "_pimp", signal_predictors, "_px1", p1, "_px2", p2,
                        "_tx1", tx1_effect, "_tx2", tx2_effect,
                        "_sigma", sigma_val, "_factstr", factor_strength,
                        "_SNR", as.integer(snr_avg))
  
  saveRDS(beta_presence_x1, file = paste0("beta_x1_pres_", sim_filename, ".RDS"))
  saveRDS(beta_presence_x2, file = paste0("beta_x2_pres_", sim_filename, ".RDS"))
  saveRDS(beta_presence_early, file = paste0("beta_early_pres_", sim_filename, ".RDS"))
  saveRDS(beta_presence_late, file = paste0("beta_late_pres_", sim_filename, ".RDS"))
  saveRDS(beta_presence_coop, file = paste0("beta_coop_pres_", sim_filename, ".RDS"))
  saveRDS(beta_presence_coop_adap, file = paste0("beta_coop_adap_pres_", sim_filename, ".RDS"))
  
  data_df <- data.frame(cbind(mse_x1, mse_x2, mse_early, mse_late, mse_coop, mse_coop_adap,
                              beta_count_x1, beta_count_x2, beta_count_early, beta_count_late, beta_count_coop, beta_count_coop_adap,
                              theta_count_x1, theta_count_x2, theta_count_early, theta_count_late, theta_count_coop, theta_count_coop_adap,
                              selected_rhos, selected_rhos_adap))
  
  write.csv(data_df, file = paste0("C:/Users/dalma/Desktop/Matteo/uni/magistrale/erasmus/thesis/", sim_filename, "_with_adap.csv"))
}
  
  
  
  
  
  