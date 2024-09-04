

library(caret)

cutoff = 3 # features selected more than - cutoff - times

simN <- 10

theta_matrices <- readRDS(paste0("labor_onset_theta_matrices.RDS"))
beta_matrices <- readRDS(paste0("labor_onset_beta_matrices.RDS"))
test_err <- readRDS(paste0("labor_onset_testerr.rds"))
labor_onset_mat <- readRDS(paste0("labor_onset_results.rds"))


load("labor_onset_data.rda")


# Function to extract the second character from a string
extract_second_char <- function(string) {
  return(substr(string, 2, 2))
}

y = DOS
X1 = Proteomics
X2 = Metabolomics
Z = cbind(as.integer(Timepoint=='G1'),as.integer(Timepoint=='G3'))

all_rows <- c(colnames(X1),colnames(X2))

# selected betas
{
  beta_matrix_full <- rep(0,length(all_rows))
  names(beta_matrix_full) <- all_rows
  
avg_abs_beta_x1 = list()
nonzero_abs_beta_x1 = list()

for(j in 1:simN){
  avg_abs_beta_x1[[j]] <- beta_matrix_full
  avg_abs_beta_x1[[j]][names(beta_matrices[[1 + (j-1)*6]])] <- beta_matrices[[1 + (j-1)*6]]
  
  nonzero_abs_beta_x1[[j]] <- beta_matrix_full
  nonzero_abs_beta_x1[[j]][names(beta_matrices[[1 + (j-1)*6]])] <- beta_matrices[[1 + (j-1)*6]]!=0
}

coefs_beta_x1 = Reduce('+', avg_abs_beta_x1)/10
counts_beta_x1 = Reduce('+', nonzero_abs_beta_x1)

avg_abs_beta_x2 = list()
nonzero_abs_beta_x2 = list()

for(j in 1:simN){
  avg_abs_beta_x2[[j]] <- beta_matrix_full
  avg_abs_beta_x2[[j]][names(beta_matrices[[2 + (j-1)*6]])] <- beta_matrices[[2 + (j-1)*6]]
  
  nonzero_abs_beta_x2[[j]] <- beta_matrix_full
  nonzero_abs_beta_x2[[j]][names(beta_matrices[[2 + (j-1)*6]])] <- beta_matrices[[2 + (j-1)*6]]!=0
}

coefs_beta_x2 = Reduce('+', avg_abs_beta_x2)/10
counts_beta_x2 = Reduce('+', nonzero_abs_beta_x2)

avg_abs_beta_early = list()
nonzero_abs_beta_early = list()

for(j in 1:simN){
  avg_abs_beta_early[[j]] <- beta_matrix_full
  avg_abs_beta_early[[j]][names(beta_matrices[[3 + (j-1)*6]])] <- beta_matrices[[3 + (j-1)*6]]
  
  nonzero_abs_beta_early[[j]] <- beta_matrix_full
  nonzero_abs_beta_early[[j]][names(beta_matrices[[3 + (j-1)*6]])] <- beta_matrices[[3 + (j-1)*6]]!=0
}

coefs_beta_early = Reduce('+', avg_abs_beta_early)/10
counts_beta_early = Reduce('+', nonzero_abs_beta_early)

avg_abs_beta_late = list()
nonzero_abs_beta_late = list()

for(j in 1:simN){
  avg_abs_beta_late[[j]] <- beta_matrix_full
  avg_abs_beta_late[[j]][names(beta_matrices[[4 + (j-1)*6]])] <- beta_matrices[[4 + (j-1)*6]]
  
  nonzero_abs_beta_late[[j]] <- beta_matrix_full
  nonzero_abs_beta_late[[j]][names(beta_matrices[[4 + (j-1)*6]])] <- beta_matrices[[4 + (j-1)*6]]!=0
}

coefs_beta_late = Reduce('+', avg_abs_beta_late)/10
counts_beta_late = Reduce('+', nonzero_abs_beta_late)

avg_abs_beta_coop = list()
nonzero_abs_beta_coop = list()

for(j in 1:simN){
  avg_abs_beta_coop[[j]] <- beta_matrix_full
  avg_abs_beta_coop[[j]][names(beta_matrices[[5 + (j-1)*6]])] <- beta_matrices[[5 + (j-1)*6]]
  
  nonzero_abs_beta_coop[[j]] <- beta_matrix_full
  nonzero_abs_beta_coop[[j]][names(beta_matrices[[5 + (j-1)*6]])] <- beta_matrices[[5 + (j-1)*6]]!=0
}

coefs_beta_coop = Reduce('+', avg_abs_beta_coop)/10
counts_beta_coop = Reduce('+', nonzero_abs_beta_coop)

avg_abs_beta_coop_adap = list()
nonzero_abs_beta_coop_adap = list()

for(j in 1:simN){
  avg_abs_beta_coop_adap[[j]] <- beta_matrix_full
  avg_abs_beta_coop_adap[[j]][names(beta_matrices[[6 + (j-1)*6]])] <- beta_matrices[[6 + (j-1)*6]]
  
  nonzero_abs_beta_coop_adap[[j]] <- beta_matrix_full
  nonzero_abs_beta_coop_adap[[j]][names(beta_matrices[[6 + (j-1)*6]])] <- as.numeric(beta_matrices[[6 + (j-1)*6]]!=0)
}

coefs_beta_coop_adap = Reduce('+', avg_abs_beta_coop_adap)/10
counts_beta_coop_adap = Reduce('+', nonzero_abs_beta_coop_adap)
}


{
## selected thetas

theta_matrix_full <- matrix(0, nrow = length(all_rows), ncol = 2, dimnames = list(all_rows, c(1,2)))

avg_abs_theta_x1 = list()
nonzero_abs_theta_x1 = list()

for(j in 1:simN){
  avg_abs_theta_x1[[j]] <- theta_matrix_full
  avg_abs_theta_x1[[j]][rownames(theta_matrices[[1 + (j-1)*6]]),] <- theta_matrices[[1 + (j-1)*6]]
  
  nonzero_abs_theta_x1[[j]] <- theta_matrix_full
  nonzero_abs_theta_x1[[j]][rownames(theta_matrices[[1 + (j-1)*6]]),] <- theta_matrices[[1 + (j-1)*6]]!=0
}

coefs_theta_x1 = Reduce('+',avg_abs_theta_x1)/10
counts_theta_x1 = Reduce('+',nonzero_abs_theta_x1)

  
avg_abs_theta_x2 = list()
nonzero_abs_theta_x2 = list()

for(j in 1:simN){
  avg_abs_theta_x2[[j]] <- theta_matrix_full
  avg_abs_theta_x2[[j]][rownames(theta_matrices[[2 + (j-1)*6]]),] <- theta_matrices[[2 + (j-1)*6]]
  
  nonzero_abs_theta_x2[[j]] <- theta_matrix_full
  nonzero_abs_theta_x2[[j]][rownames(theta_matrices[[2 + (j-1)*6]]),] <- theta_matrices[[2 + (j-1)*6]]!=0
}

coefs_theta_x2 = Reduce('+',avg_abs_theta_x2)/10
counts_theta_x2 = Reduce('+',nonzero_abs_theta_x2)


avg_abs_theta_early = list()
nonzero_abs_theta_early = list()

for(j in 1:simN){
  avg_abs_theta_early[[j]] <- theta_matrix_full
  avg_abs_theta_early[[j]][rownames(theta_matrices[[3 + (j-1)*6]]),] <- theta_matrices[[3 + (j-1)*6]]
  
  nonzero_abs_theta_early[[j]] <- theta_matrix_full
  nonzero_abs_theta_early[[j]][rownames(theta_matrices[[3 + (j-1)*6]]),] <- theta_matrices[[3 + (j-1)*6]]!=0
}

coefs_theta_early = Reduce('+',avg_abs_theta_early)/10
counts_theta_early = Reduce('+',nonzero_abs_theta_early)


avg_abs_theta_late = list()
nonzero_abs_theta_late = list()

for(j in 1:simN){
  avg_abs_theta_late[[j]] <- theta_matrix_full
  avg_abs_theta_late[[j]][rownames(theta_matrices[[4 + (j-1)*6]]),] <- theta_matrices[[4 + (j-1)*6]]
  
  nonzero_abs_theta_late[[j]] <- theta_matrix_full
  nonzero_abs_theta_late[[j]][rownames(theta_matrices[[4 + (j-1)*6]]),] <- theta_matrices[[4 + (j-1)*6]]!=0
}

coefs_theta_late = Reduce('+',avg_abs_theta_late)/10
counts_theta_late = Reduce('+',nonzero_abs_theta_late)


avg_abs_theta_coop = list()
nonzero_abs_theta_coop = list()

for(j in 1:simN){
  avg_abs_theta_coop[[j]] <- theta_matrix_full
  avg_abs_theta_coop[[j]][rownames(theta_matrices[[5 + (j-1)*6]]),] <- theta_matrices[[5 + (j-1)*6]]
  
  nonzero_abs_theta_coop[[j]] <- theta_matrix_full
  nonzero_abs_theta_coop[[j]][rownames(theta_matrices[[5 + (j-1)*6]]),] <- theta_matrices[[5 + (j-1)*6]]!=0
}

coefs_theta_coop = Reduce('+',avg_abs_theta_coop)/10
counts_theta_coop = Reduce('+',nonzero_abs_theta_coop)


avg_abs_theta_coop_adap = list()
nonzero_abs_theta_coop_adap = list()

for(j in 1:simN){
  avg_abs_theta_coop_adap[[j]] <- theta_matrix_full
  avg_abs_theta_coop_adap[[j]][rownames(theta_matrices[[6 + (j-1)*6]]),] <- theta_matrices[[6 + (j-1)*6]]
  
  nonzero_abs_theta_coop_adap[[j]] <- theta_matrix_full
  nonzero_abs_theta_coop_adap[[j]][rownames(theta_matrices[[6 + (j-1)*6]]),] <- theta_matrices[[6 + (j-1)*6]]!=0
}

coefs_theta_coop_adap = Reduce('+',avg_abs_theta_coop_adap)/10
counts_theta_coop_adap = Reduce('+',nonzero_abs_theta_coop_adap)

}



beta_sel <- c(sum(counts_beta_x1>cutoff), sum(counts_beta_x2>cutoff), sum(counts_beta_early>cutoff), sum(counts_beta_late>cutoff), sum(counts_beta_coop>cutoff), sum(counts_beta_coop_adap>cutoff))
theta_sel <- c(sum(counts_theta_x1>cutoff), sum(counts_theta_x2>cutoff), sum(counts_theta_early>cutoff), sum(counts_theta_late>cutoff), sum(counts_theta_coop>cutoff), sum(counts_theta_coop_adap>cutoff))

results_table = data.frame(colMeans(test_err), apply(test_err,2,sd)/sqrt(10), beta_sel,theta_sel)
colnames(results_table) <- c('Test MSE', 'SD', 'Main effects', 'Interactions')
rownames(results_table) <- c('Only X1', 'Only X2', 'Early fusion', 'Late fusion', 'CoopPliable','AdapCoopPliable')
print(results_table)
