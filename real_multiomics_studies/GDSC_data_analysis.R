
library(caret)
library(pliable)
source('functions.R')


## this script performs the analysis of GDSC data

## we study ten different train-test splits and get results on test MSE and features selected


# Load the GDSC data
GDSC<-readRDS("GDSC_data_processed.rds")

# Select the specific drug for analysis
selected_drug <- "Nilotinib"

# Extract the drug's pharmacological profiling and tissue dummy variables
col_filter <- colnames(GDSC$y) %in% paste("IC50.", selected_drug, sep = "")
response_and_tissue <- cbind(
  GDSC$y[, col_filter],
  GDSC$z[, 1:12]
)

# Prepare feature matrices (X1: expression data, X2: mutation data) and response variable (y)
X1 <- log2(GDSC$x1)  # Log2 transformation for expression levels
X2 <- GDSC$x2
y <- GDSC$y[, col_filter]
Z <- GDSC$z[, -6]
response_name <- stringr::str_remove(colnames(y), pattern = "IC50.")

# Set row names for feature matrices and response variable
rownames(X1) <- GDSC$names.cell.line
rownames(X2) <- GDSC$names.cell.line
rownames(Z) <- GDSC$names.cell.line

# Define the number of features and samples
num.nonpen <- GDSC$num.nonpen
p1 = ncol(X1)
p2 = ncol(X2)
n = nrow(X1)
num_modifiers = ncol(Z)
num_folds = 10
num_simulations = 10
train_fraction = 0.75
sim_seed = 0
validation_fraction = 0.4
set.seed(sim_seed)

# Define the sequence of rho values to explore
rho_values = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 3)

# Initialize vectors to store results
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

chosen_rhos = numeric(num_simulations)
chosen_rhos_adap = numeric(num_simulations)

# Initialize lists to store beta and theta matrices
beta_matrices = list()
theta_matrices = list()

# Combine Z (modifiers), X1 (expression), and X2 (mutation) into a single matrix
GDSC$x <- cbind(Z, X1, X2)

# Loop over each simulation
for (sim in seq_len(num_simulations)) {
  
  cat("Simulation:", sim, "\n")
  
  # Select train dataset of 80% cell lines of each tissue
  trainID <- NULL
  for(i.tissue in 1:num.nonpen) trainID <- c(trainID, sample(which(GDSC$x[,i.tissue]==1),round(sum(GDSC$x[,i.tissue]==1)*0.8)))
  
  # make sure at least one mutated cell lines of each cancer tissue for every mutation features
  repeat{
    if(min(colSums(GDSC$x[trainID, num.nonpen+p1:(p1+p2-1)]))>=1) break
    trainID <- NULL
    for(i.tissue in 1:num.nonpen) trainID <- c(trainID, sample(which(GDSC$x[,i.tissue]==1),round(sum(GDSC$x[,i.tissue]==1)*0.8)))
  }
  
  # Split the data into training and testing sets
  y_train = y[trainID]
  Z_train = Z[trainID, ]
  X1_train = X1[trainID, ]
  X2_train = X2[trainID, ]
  
  y_test = y[-trainID]
  Z_test = Z[-trainID, ]
  X1_test = X1[-trainID, ]
  X2_test = X2[-trainID, ]
  
  # Preprocess (center and scale) the training data and apply the same transformation to the test data
  preprocess_X1 = preProcess(X1_train, method = c("center", "scale"))
  X1_train = predict(preprocess_X1, X1_train)
  X1_test = predict(preprocess_X1, X1_test)
  
  preprocess_X2 = preProcess(X2_train, method = c("center", "scale"))
  X2_train = predict(preprocess_X2, X2_train)
  X2_test = predict(preprocess_X2, X2_test)
  
  # Center the response variable based on the training data
  y_test = y_test - mean(y_train)
  y_train = y_train - mean(y_train)
  
  # Variable names for later use
  var_names_x1 = colnames(X1_train)
  var_names_x2 = colnames(X2_train)
  var_names_total = c(var_names_x1, var_names_x2)
  
  # Create folds for cross-validation
  fold_ids = sample(rep_len(1:num_folds, nrow(X1_train)))
  
  # Null model (mean of the training response)
  cat("Null model\n")
  mse_null[sim] = calc_mse(mean(y_train), y_test)
  cat("MSE Null:", mse_null[sim], "\n")
  
  # Model with only X1
  cat("Only X1\n")
  
  fit_X1 = pliable(x = X1_train, z = Z_train, y = y_train, nlambda = 50, zlinear = FALSE, alpha = 0.1)
  cv_X1 = cv.pliable(fit_X1, x = X1_train, z = Z_train, y = y_train, foldid = fold_ids, verbose = FALSE)
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
  cv_X2 = cv.pliable(fit_X2, x = X2_train, z = Z_train, y = y_train, foldid = fold_ids, verbose = FALSE)
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
  cv_early = cv.pliable(fit_early, x = combined_train_X, z = Z_train, y = y_train, foldid = fold_ids, verbose = FALSE)
  
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
  
  # Select validation set and remaining training set
  num_val_samples = floor(validation_fraction * nrow(X1_train))
  val_indices = sort(sample(seq_len(nrow(X1_train)), size = num_val_samples))
  train_late_indices = setdiff(seq_len(nrow(X1_train)), val_indices)
  
  X1_val = X1_train[val_indices, ]
  X1_train_late = X1_train[train_late_indices, ]
  X2_val = X2_train[val_indices, ]
  X2_train_late = X2_train[train_late_indices, ]
  Z_val = Z_train[val_indices, ]
  Z_train_late = Z_train[train_late_indices, ]
  y_val = y_train[val_indices]
  y_train_late = y_train[train_late_indices]
  
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
  
  beta_coops_rho = list()
  theta_coops_rho = list()
  
  for (rho_idx in seq_along(rho_values)) {
    current_rho = rho_values[rho_idx]
    cat("Rho:", current_rho, "\n")
    
    coop_fit = coop_pliable(X1_train, X2_train, Z_train, y_train,
                            coop_alpha = current_rho, fold_indices = fold_ids,
                            num_folds = max(fold_ids), pliable_alpha = 0.1)
    
    cv_mse[rho_idx] = min(coop_fit$cvm)
    beta_selected_coop_rho[rho_idx] = coop_fit$best_fit$nbeta[coop_fit$min_ind]
    theta_selected_coop_rho[rho_idx] = coop_fit$best_fit$ntheta[coop_fit$min_ind]
    mse_test[rho_idx] = calc_mse(y_test, predict(coop_fit$best_fit, cbind(X1_test, X2_test), Z_test)[, coop_fit$min_ind])
    cat("MSE:", mse_test[rho_idx], "\n")
    
    beta_coops_rho[[rho_idx]] = setNames(coop_fit$best_fit$beta[, coop_fit$min_ind], var_names_total)
    theta_coops_rho[[rho_idx]] = coop_fit$best_fit$theta[, , coop_fit$min_ind]
    rownames(theta_coops_rho[[rho_idx]]) = var_names_total
  }
  
  selected_rho = rho_values[which.min(cv_mse)]
  cat("Selected Rho:", selected_rho, "\n")
  
  beta_selected_coop[sim] = beta_selected_coop_rho[which.min(cv_mse)]
  theta_selected_coop[sim] = theta_selected_coop_rho[which.min(cv_mse)]
  
  mse_coop[sim] = mse_test[which.min(cv_mse)]
  chosen_rhos[sim] = selected_rho
  
  beta_matrices[[5 + (sim - 1) * 6]] = beta_coops_rho[[which.min(cv_mse)]]
  theta_matrices[[5 + (sim - 1) * 6]] = theta_coops_rho[[which.min(cv_mse)]]
  
  
  
  
  # Adaptive Cooperative Learning Model
  cat("AdapCoopLearn\n")
  
  cv_mse_adap = numeric(length(rho_values))
  mse_test_adap = numeric(length(rho_values))
  
  beta_selected_coop_adap_rho = numeric(length(rho_values))
  theta_selected_coop_adap_rho = numeric(length(rho_values))
  
  beta_coops_adap_rho = list()
  theta_coops_adap_rho = list()
  
  for (rho_idx in seq_along(rho_values)) {
    current_rho = rho_values[rho_idx]
    cat("Rho:", current_rho, "\n")
    
    coop_adap_fit = coop_pliable(X1_train, X2_train, Z_train, y_train,
                                 coop_alpha = current_rho, fold_indices = fold_ids,
                                 num_folds = max(fold_ids), pliable_alpha = 0.1,
                                 penalty_factors = c(rep(lambda1, p1), rep(lambda2, p1)))
    
    
    cv_mse_adap[rho_idx] = min(coop_adap_fit$cvm)
    beta_selected_coop_adap_rho[rho_idx] = coop_adap_fit$best_fit$nbeta[coop_adap_fit$min_ind]
    theta_selected_coop_adap_rho[rho_idx] = coop_adap_fit$best_fit$ntheta[coop_adap_fit$min_ind]
    mse_test_adap[rho_idx] = calc_mse(y_test, predict(coop_adap_fit$best_fit, cbind(X1_test, X2_test), Z_test)[, coop_adap_fit$min_ind])
    cat("MSE:", mse_test_adap[rho_idx], "\n")
    
    beta_coops_adap_rho[[rho_idx]] = setNames(coop_adap_fit$best_fit$beta[, coop_adap_fit$min_ind], var_names_total)
    theta_coops_adap_rho[[rho_idx]] = coop_adap_fit$best_fit$theta[, , coop_adap_fit$min_ind]
    rownames(theta_coops_adap_rho[[rho_idx]]) = var_names_total
  }
  
  selected_rho_adap = rho_values[which.min(cv_mse_adap)]
  cat("Selected Rho:", selected_rho_adap, "\n")
  
  beta_selected_coop_adap[sim] = beta_selected_coop_adap_rho[which.min(cv_mse_adap)]
  theta_selected_coop_adap[sim] = theta_selected_coop_adap_rho[which.min(cv_mse_adap)]
  
  mse_coop_adap[sim] = mse_test_adap[which.min(cv_mse_adap)]
  chosen_rhos_adap[sim] = selected_rho_adap
  
  beta_matrices[[6 + (sim - 1) * 6]] = beta_coops_adap_rho[[which.min(cv_mse_adap)]]
  theta_matrices[[6 + (sim - 1) * 6]] = theta_coops_adap_rho[[which.min(cv_mse_adap)]]  
 }


test_err = cbind(mse_x1, mse_x2, mse_early, mse_late, mse_coop, mse_coop_adap)
beta_selection = cbind(beta_selected_x1, beta_selected_x2, beta_selected_early, beta_selected_late, beta_selected_coop, beta_selected_coop_adap)
theta_selection = cbind(theta_selected_x1, theta_selected_x2, theta_selected_early, theta_selected_late, theta_selected_coop, theta_selected_coop_adap)


GDSC_data_mat = rbind(colMeans(test_err), apply(test_err,2,sd)/sqrt(10),
                      colMeans(beta_selection), apply(beta_selection,2,sd)/sqrt(10),
                      colMeans(theta_selection), apply(theta_selection,2,sd)/sqrt(10))


saveRDS(test_err, paste0("GDSC_testerr.rds"))
saveRDS(beta_matrices, "GDSC_beta_matrices.RDS")
saveRDS(theta_matrices, "GDSC_theta_matrices.RDS")

saveRDS(GDSC_data_mat, "GDSC_data_mat.RDS")

write(c(chosen_rhos), file = paste0("GDSC_chosen_rhos.txt"))
write(c(chosen_rhos_adap), file = paste0("GDSC_chosen_rhos_adap.txt"))
