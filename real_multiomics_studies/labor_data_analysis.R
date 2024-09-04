


# Load necessary libraries and scripts
library(caret)
library(glmnet)
library(pliable)
source('functions.R')

# Load the labor onset data
load("labor_onset_data.rda")

# Set the seed
sim_seed = 100
set.seed(sim_seed)

# Number of simulations to run
num_simulations = 10

# Initialize variables to store results
mse_null = numeric(num_simulations)
mse_x1 = numeric(num_simulations)
mse_x2 = numeric(num_simulations)
mse_early = numeric(num_simulations)
mse_late = numeric(num_simulations)
mse_coop = numeric(num_simulations)
mse_coop_adap = numeric(num_simulations)

beta_selected_x1 = numeric(num_simulations)
beta_selected_x2 = numeric(num_simulations)
beta_selected_early = numeric(num_simulations)
beta_selected_late = numeric(num_simulations)
beta_selected_coop = numeric(num_simulations)
beta_selected_coop_adap = numeric(num_simulations)

theta_selected_x1 = numeric(num_simulations)
theta_selected_x2 = numeric(num_simulations)
theta_selected_early = numeric(num_simulations)
theta_selected_late = numeric(num_simulations)
theta_selected_coop = numeric(num_simulations)
theta_selected_coop_adap = numeric(num_simulations)


# Store beta and theta matrices for each simulation
beta_matrices = list()
theta_matrices = list()

selected_rhos = numeric(num_simulations)
selected_rhos_adap = numeric(num_simulations)

# Define response variable and predictors
response_y = DOS
predictor_X1 = Proteomics
predictor_X2 = Metabolomics
modifier_Z = cbind(as.integer(Timepoint == 'G1'), as.integer(Timepoint == 'G3'))
colnames(modifier_Z) <- c('T1', 'T3')

num_features_x1 = ncol(predictor_X1)
num_features_x2 = ncol(predictor_X2)
num_samples = nrow(predictor_X1)
num_modifiers = ncol(modifier_Z)

# Define cross-validation and split parameters
num_folds = 10
train_fraction = 0.8
validation_fraction = 0.4

# Define rho values to be explored
rho_values = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 3)


# Loop over each simulation
for (sim in seq_len(num_simulations)) {
  
  cat("Simulation:", sim, "\n")
  set.seed(sim_seed + sim)
  
  # Split data based on unique IDs
  unique_ids = unique(Id)
  num_train_samples = floor(train_fraction * length(unique_ids))
  train_ids = sort(sample(unique_ids, size = num_train_samples))
  train_indices = which(Id %in% train_ids)
  test_indices = setdiff(seq_len(num_samples), train_indices)
  
  # Subset the data into training and testing sets
  train_X1_raw = predictor_X1[train_indices, ]
  test_X1_raw = predictor_X1[test_indices, ]
  train_X2_raw = predictor_X2[train_indices, ]
  test_X2_raw = predictor_X2[test_indices, ]
  train_Z_raw = modifier_Z[train_indices, ]
  test_Z_raw = modifier_Z[test_indices, ]
  
  # Preprocess the training data (center and scale)
  preprocess_X1 = preProcess(train_X1_raw, method = c("center", "scale"))
  X1_train = predict(preprocess_X1, train_X1_raw)
  X1_test = predict(preprocess_X1, test_X1_raw)
  
  preprocess_X2 = preProcess(train_X2_raw, method = c("center", "scale"))
  X2_train = predict(preprocess_X2, train_X2_raw)
  X2_test = predict(preprocess_X2, test_X2_raw)
  
  # Use raw Z data
  Z_train = train_Z_raw
  Z_test = test_Z_raw
  
  # Center the response variable based on the training data
  y_train = response_y[train_indices]
  y_test = response_y[test_indices]
  
  y_test = y_test - mean(y_train)
  y_train = y_train - mean(y_train)
  
  
  # Get variable names for use in beta and theta matrices
  var_names_x1 = colnames(X1_train)
  var_names_x2 = colnames(X2_train)
  var_names_total = c(var_names_x1, var_names_x2)
  
  # Create folds based on unique IDs
  fold_ids = sample(rep_len(1:num_folds, length(unique_ids)))
  fold_assignments = fold_ids[Id[train_indices]]
  
  # Null model (mean of the training response)
  cat("Null model\n")
  mse_null[sim] = calc_mse(mean(y_train), y_test)
  cat("MSE Null:", mse_null[sim], "\n")

  # Model with only X1
  cat("Only X1\n")

  fit_X1 = pliable(x = X1_train, z = Z_train, y = y_train, nlambda = 50, zlinear = FALSE, alpha = 0.1)
  cv_X1 = cv.pliable(fit_X1, x = X1_train, z = Z_train, y = y_train, foldid = fold_assignments, verbose = FALSE)
  predictions_X1 = predict(fit_X1, x = X1_test, z = Z_test, lambda = cv_X1$lambda.min)

  mse_x1[sim] = calc_mse(y_test, predictions_X1)
  beta_selected_x1[sim] = fit_X1$nbeta[which(cv_X1$lambda == cv_X1$lambda.min)]
  theta_selected_x1[sim] = fit_X1$ntheta[which(cv_X1$lambda == cv_X1$lambda.min)]
  cat("MSE X1:", mse_x1[sim], "\n")

  beta_matrices[[1 + (sim - 1) * 6]] = setNames(fit_X1$beta[, which(cv_X1$lambda == cv_X1$lambda.min)], var_names_x1)
  theta_matrices[[1 + (sim - 1) * 6]] = fit_X1$theta[, , which(cv_X1$lambda == cv_X1$lambda.min)]
  rownames(theta_matrices[[1 + (sim - 1) * 6]]) = var_names_x1

  lambda1 = cv_X1$lambda.min
  cat("Lambda X1:", lambda1, "\n")

  # Model with only X2
  cat("Only X2\n")

  fit_X2 = pliable(x = X2_train, z = Z_train, y = y_train, nlambda = 50, zlinear = FALSE, alpha = 0.1)
  cv_X2 = cv.pliable(fit_X2, x = X2_train, z = Z_train, y = y_train, foldid = fold_assignments, verbose = FALSE)
  predictions_X2 = predict(fit_X2, x = X2_test, z = Z_test, lambda = cv_X2$lambda.min)

  mse_x2[sim] = calc_mse(y_test, predictions_X2)
  beta_selected_x2[sim] = fit_X2$nbeta[which(cv_X2$lambda == cv_X2$lambda.min)]
  theta_selected_x2[sim] = fit_X2$ntheta[which(cv_X2$lambda == cv_X2$lambda.min)]
  cat("MSE X2:", mse_x2[sim], "\n")

  beta_matrices[[2 + (sim - 1) * 6]] = setNames(fit_X2$beta[, which(cv_X2$lambda == cv_X2$lambda.min)], var_names_x2)
  theta_matrices[[2 + (sim - 1) * 6]] = fit_X2$theta[, , which(cv_X2$lambda == cv_X2$lambda.min)]
  rownames(theta_matrices[[2 + (sim - 1) * 6]]) = var_names_x2

  lambda2 = cv_X2$lambda.min
  cat("Lambda X2:", lambda2, "\n")

  # Early Fusion Model (combining X1 and X2)
  cat("Early Fusion\n")

  combined_train_X = cbind(X1_train, X2_train)
  combined_test_X = cbind(X1_test, X2_test)

  fit_early = pliable(x = combined_train_X, z = Z_train, y = y_train, nlambda = 50, zlinear = FALSE, alpha = 0.1)
  cv_early = cv.pliable(fit_early, x = combined_train_X, z = Z_train, y = y_train, foldid = fold_assignments, verbose = FALSE)

  predictions_early = predict(fit_early, x = combined_test_X, z = Z_test, lambda = cv_early$lambda.min)

  mse_early[sim] = calc_mse(y_test, predictions_early)
  beta_selected_early[sim] = fit_early$nbeta[which(cv_early$lambda == cv_early$lambda.min)]
  theta_selected_early[sim] = fit_early$ntheta[which(cv_early$lambda == cv_early$lambda.min)]
  cat("MSE Early:", mse_early[sim], "\n")

  beta_matrices[[3 + (sim - 1) * 6]] = setNames(fit_early$beta[, which(cv_early$lambda == cv_early$lambda.min)], var_names_total)
  theta_matrices[[3 + (sim - 1) * 6]] = fit_early$theta[, , which(cv_early$lambda == cv_early$lambda.min)]
  rownames(theta_matrices[[3 + (sim - 1) * 6]]) = var_names_total

  # Late Fusion Model
  cat("Late Fusion\n")
  
  second_stage_samples = floor(validation_fraction * num_train_samples)
  val_ids = sort(sample(train_ids, size = second_stage_samples))
  val_indices = which(Id %in% val_ids)
  train_late_indices = setdiff(seq_len(nrow(X1)), c(test_indices,val_indices))
  
  X1_val_raw = predictor_X1[val_indices,]
  X1_train_late_raw = predictor_X1[train_late_indices,]
  X2_val_raw = predictor_X2[val_indices,]
  X2_train_late_raw = predictor_X2[train_late_indices,]
  Z_val = modifier_Z[val_indices,]
  Z_train_late = modifier_Z[train_late_indices,]
  y_val = response_y[val_indices]
  y_train_late = response_y[train_late_indices]
  
  y_val <- y_val - mean(y_train_late)
  y_train_late <- y_train_late - mean(y_train_late)
  
  preprocess_values_train_X1_late = preProcess(X1_train_late_raw, method = c("center", "scale"))
  X1_train_late = predict(preprocess_values_train_X1_late, X1_train_late_raw)
  X1_val = predict(preprocess_values_train_X1_late, X1_val_raw)
  X1_test_late = predict(preprocess_values_train_X1_late, test_X1_raw)
  
  preprocess_values_train_X2_late = preProcess(X2_train_late_raw, method = c("center", "scale"))
  X2_train_late = predict(preprocess_values_train_X2_late, X2_train_late_raw)
  X2_val = predict(preprocess_values_train_X2_late, X2_val_raw)
  X2_test_late = predict(preprocess_values_train_X2_late, test_X2_raw)

  fit_X1_late = pliable(X1_train_late, Z_train_late, y_train_late, nlambda = 50, zlinear = FALSE, alpha = 0.1)
  cv_X1_late = cv.pliable(fit_X1_late, X1_train_late, Z_train_late, y_train_late, nfolds = 5, verbose = FALSE)

  fit_X2_late = pliable(X2_train_late, Z_train_late, y_train_late, nlambda = 50, zlinear = FALSE, alpha = 0.1)
  cv_X2_late = cv.pliable(fit_X2_late, X2_train_late, Z_train_late, y_train_late, nfolds = 5, verbose = FALSE)

  X1_pred_val = predict(fit_X1_late, x = X1_val, z = Z_val, lambda = cv_X1_late$lambda.min)
  X1_pred_test = predict(fit_X1_late, x = X1_test, z = Z_test, lambda = cv_X1_late$lambda.min)
  X2_pred_val = predict(fit_X2_late, x = X2_val, z = Z_val, lambda = cv_X2_late$lambda.min)
  X2_pred_test = predict(fit_X2_late, x = X2_test, z = Z_test, lambda = cv_X2_late$lambda.min)

  fusion_data = data.frame(y = y_val, X1_pred = as.vector(X1_pred_val), X2_pred = as.vector(X2_pred_val))
  fusion_fit = lm(y ~ X1_pred + X2_pred, data = fusion_data)
  predictions_late = predict(fusion_fit, data.frame(X1_pred = as.vector(X1_pred_test), X2_pred = as.vector(X2_pred_test)))

  mse_late[sim] = calc_mse(y_test, predictions_late)
  beta_selected_late[sim] = fit_X1_late$nbeta[which(cv_X1_late$lambda == cv_X1_late$lambda.min)] +
    fit_X2_late$nbeta[which(cv_X2_late$lambda == cv_X2_late$lambda.min)]
  theta_selected_late[sim] = fit_X1_late$ntheta[which(cv_X1_late$lambda == cv_X1_late$lambda.min)] +
    fit_X2_late$ntheta[which(cv_X2_late$lambda == cv_X2_late$lambda.min)]
  cat("MSE Late:", mse_late[sim], "\n")

  beta_matrices[[4 + (sim - 1) * 6]] = setNames(c(fit_X1_late$beta[, which(cv_X1_late$lambda == cv_X1_late$lambda.min)],
                                                  fit_X2_late$beta[, which(cv_X2_late$lambda == cv_X2_late$lambda.min)]),
                                                var_names_total)

  theta_matrices[[4 + (sim - 1) * 6]] = rbind(fit_X1_late$theta[, , which(cv_X1_late$lambda == cv_X1_late$lambda.min)],
                                              fit_X2_late$theta[, , which(cv_X2_late$lambda == cv_X2_late$lambda.min)])
  rownames(theta_matrices[[4 + (sim - 1) * 6]]) = var_names_total


  # Cooperative Learning Model
  cat("CoopLearn\n")

  cv_mse = numeric(length(rho_values))
  mse_test = numeric(length(rho_values))

  beta_selected_coop_rho = numeric(length(rho_values))
  theta_selected_coop_rho = numeric(length(rho_values))

  beta_coop_rho = list()
  theta_coop_rho = list()

  for (rho_idx in seq_along(rho_values)) {
    current_rho = rho_values[rho_idx]
    cat("Rho:", current_rho, "\n")

    fit_coop_pliable = coop_pliable(X1_train, X2_train, Z_train, y_train,
                                    coop_alpha = current_rho, fold_indices = fold_assignments,
                                    num_folds = max(fold_ids), pliable_alpha = 0.1)


    cv_mse[rho_idx] = min(fit_coop_pliable$cvm)
    beta_selected_coop_rho[rho_idx] = fit_coop_pliable$best_fit$nbeta[fit_coop_pliable$min_ind]
    theta_selected_coop_rho[rho_idx] = fit_coop_pliable$best_fit$ntheta[fit_coop_pliable$min_ind]
    mse_test[rho_idx] = calc_mse(y_test, predict(fit_coop_pliable$best_fit, cbind(X1_test, X2_test), Z_test)[, fit_coop_pliable$min_ind])
    cat("MSE:", mse_test[rho_idx], "\n")

    beta_coop_rho[[rho_idx]] = setNames(fit_coop_pliable$best_fit$beta[, fit_coop_pliable$min_ind], var_names_total)
    theta_coop_rho[[rho_idx]] = fit_coop_pliable$best_fit$theta[, , fit_coop_pliable$min_ind]
    rownames(theta_coop_rho[[rho_idx]]) = var_names_total
  }

  selected_rho = rho_values[which.min(cv_mse)]
  cat("Selected Rho:", selected_rho, "\n")

  beta_selected_coop[sim] = beta_selected_coop_rho[which.min(cv_mse)]
  theta_selected_coop[sim] = theta_selected_coop_rho[which.min(cv_mse)]

  mse_coop[sim] = mse_test[which.min(cv_mse)]
  selected_rhos[sim] = selected_rho

  beta_matrices[[5 + (sim - 1) * 6]] = beta_coop_rho[[which.min(cv_mse)]]
  theta_matrices[[5 + (sim - 1) * 6]] = theta_coop_rho[[which.min(cv_mse)]]

  # Adaptive Cooperative Learning Model
  cat("AdapLearn\n")

  cv_mse_adap = numeric(length(rho_values))
  mse_test_adap = numeric(length(rho_values))

  beta_selected_coop_rho_adap = numeric(length(rho_values))
  theta_selected_coop_rho_adap = numeric(length(rho_values))

  beta_coop_rho_adap = list()
  theta_coop_rho_adap = list()

  for (rho_idx in seq_along(rho_values)) {
    current_rho = rho_values[rho_idx]
    cat("Rho:", current_rho, "\n")

    fit_coop_pliable_adap = coop_pliable(X1_train, X2_train, Z_train, y_train,
                                         coop_alpha = current_rho, fold_indices = fold_assignments,
                                         num_folds = max(fold_ids), pliable_alpha = 0.1,
                                         penalty_factors = c(rep(lambda1, num_features_x1), rep(lambda2, num_features_x2)))





    cv_mse_adap[rho_idx] = min(fit_coop_pliable_adap$cvm)
    beta_selected_coop_rho_adap[rho_idx] = fit_coop_pliable_adap$best_fit$nbeta[fit_coop_pliable_adap$min_ind]
    theta_selected_coop_rho_adap[rho_idx] = fit_coop_pliable_adap$best_fit$ntheta[fit_coop_pliable_adap$min_ind]
    mse_test_adap[rho_idx] = calc_mse(y_test, predict(fit_coop_pliable_adap$best_fit, cbind(X1_test, X2_test), Z_test)[, fit_coop_pliable_adap$min_ind])
    cat("MSE:", mse_test_adap[rho_idx], "\n")

    beta_coop_rho_adap[[rho_idx]] = setNames(fit_coop_pliable_adap$best_fit$beta[, fit_coop_pliable_adap$min_ind], var_names_total)
    theta_coop_rho_adap[[rho_idx]] = fit_coop_pliable_adap$best_fit$theta[, , fit_coop_pliable_adap$min_ind]
    rownames(theta_coop_rho_adap[[rho_idx]]) = var_names_total
  }

  selected_rho_adap = rho_values[which.min(cv_mse_adap)]
  cat("Selected Rho Adap:", selected_rho_adap, "\n")

  beta_selected_coop_adap[sim] = beta_selected_coop_rho_adap[which.min(cv_mse_adap)]
  theta_selected_coop_adap[sim] = theta_selected_coop_rho_adap[which.min(cv_mse_adap)]

  mse_coop_adap[sim] = mse_test_adap[which.min(cv_mse_adap)]
  selected_rhos_adap[sim] = selected_rho_adap

  beta_matrices[[6 + (sim - 1) * 6]] = beta_coop_rho_adap[[which.min(cv_mse_adap)]]
  theta_matrices[[6 + (sim - 1) * 6]] = theta_coop_rho_adap[[which.min(cv_mse_adap)]]
}

print(mean(mse_late))
print(sd(mse_late)/sqrt(10))

# Compile results into matrices
test_errors = cbind(mse_x1, mse_x2, mse_early, mse_late, mse_coop, mse_coop_adap)
beta_selection = cbind(beta_selected_x1, beta_selected_x2, beta_selected_early, beta_selected_late, beta_selected_coop, beta_selected_coop_adap)
theta_selection = cbind(theta_selected_x1, theta_selected_x2, theta_selected_early, theta_selected_late, theta_selected_coop, theta_selected_coop_adap)

labor_onset_results = rbind(colMeans(test_errors), apply(test_errors, 2, sd) / sqrt(num_simulations),
                            colMeans(beta_selection), apply(beta_selection, 2, sd) / sqrt(num_simulations),
                            colMeans(theta_selection), apply(theta_selection, 2, sd) / sqrt(num_simulations))

# Save results to files
saveRDS(labor_onset_results, paste0("labor_onset_results.rds"))
write(c(mean(selected_rhos), mean(selected_rhos_adap)), file = paste0("labor_onset_chosen_rhos.rds"))
saveRDS(test_errors, paste0("labor_onset_testerr.rds"))
saveRDS(beta_matrices, paste0("labor_onset_beta_matrices.RDS"))
saveRDS(theta_matrices, paste0("labor_onset_theta_matrices.RDS"))

