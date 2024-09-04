simN = 10
# features selected more than - cutoff - times


cutoff = 3

theta_matrices_nilo <- readRDS('GDSC_theta_matrices.RDS')
beta_matrices_nilo <- readRDS('GDSC_beta_matrices.RDS')
test_err <- readRDS('GDSC_testerr.RDS')
GDSC_data_mat <- readRDS('GDSC_data_mat.RDS')

# selected betas

avg_abs_beta_x1 = list()

for(j in 1:simN){
  avg_abs_beta_x1[[j]] = abs(beta_matrices_nilo[[1 + (j-1)*6]])!=0
}

sum_beta_x1 = tapply(unlist(avg_abs_beta_x1), names(unlist(avg_abs_beta_x1)), sum)
#selected_betas_x1 = sort(sum_beta_x1, decreasing = TRUE)

avg_abs_beta_x2 = list()

for(j in 1:simN){
  avg_abs_beta_x2[[j]] = abs(beta_matrices_nilo[[2 + (j-1)*6]])!=0
}

sum_beta_x2 = tapply(unlist(avg_abs_beta_x2), names(unlist(avg_abs_beta_x2)), sum)
#selected_betas_x2 = sort(sum_beta_x2, decreasing = TRUE)

avg_abs_beta_early = list()

for(j in 1:simN){
  avg_abs_beta_early[[j]] = abs(beta_matrices_nilo[[3 + (j-1)*6]])!=0
}

sum_beta_early = tapply(unlist(avg_abs_beta_early), names(unlist(avg_abs_beta_early)), sum)
#selected_betas_early = sort(sum_beta_early, decreasing = TRUE)

avg_abs_beta_late = list()

for(j in 1:simN){
  avg_abs_beta_late[[j]] = abs(beta_matrices_nilo[[4 + (j-1)*6]])!=0
}

sum_beta_late = tapply(unlist(avg_abs_beta_late), names(unlist(avg_abs_beta_late)), sum)
#selected_betas_late = sort(sum_beta_late, decreasing = TRUE)

avg_abs_beta_coop_nilo = list()

for(j in 1:simN){
  avg_abs_beta_coop_nilo[[j]] = beta_matrices_nilo[[5 + (j-1)*6]]!=0
}

sum_beta_coop_nilo = tapply(unlist(avg_abs_beta_coop_nilo), names(unlist(avg_abs_beta_coop_nilo)), sum)


avg_abs_beta_coop_adap_nilo = list()

for(j in 1:simN){
  avg_abs_beta_coop_adap_nilo[[j]] = beta_matrices_nilo[[6 + (j-1)*6]]!=0
}

sum_beta_coop_adap_nilo = tapply(unlist(avg_abs_beta_coop_adap_nilo), names(unlist(avg_abs_beta_coop_adap_nilo)), sum)



## selected thetas

rownames(theta_matrices_nilo[[1]]) <- names(beta_matrices_nilo[[1]])
sum_theta_x1 <- as.numeric(theta_matrices_nilo[[1]]!=0)

for(j in 2:length(avg_abs_beta_x1)){
  
  rownames(theta_matrices_nilo[[1 + (j-1)*6]]) <- names(beta_matrices_nilo[[1 + (j-1)*6]])
  
  sum_theta_x1 <- sum_theta_x1 + as.numeric(theta_matrices_nilo[[1 + (j-1)*6]]!=0)
}

rownames(theta_matrices_nilo[[2]]) <- names(beta_matrices_nilo[[2]])
sum_theta_x2 <- as.numeric(theta_matrices_nilo[[2]]!=0)

for(j in 2:length(avg_abs_beta_x2)){
  
  rownames(theta_matrices_nilo[[2 + (j-1)*6]]) <- names(beta_matrices_nilo[[2 + (j-1)*6]])
  
  sum_theta_x2 <- sum_theta_x2 + as.numeric(theta_matrices_nilo[[2 + (j-1)*6]]!=0)
}

rownames(theta_matrices_nilo[[3]]) <- names(beta_matrices_nilo[[3]])
sum_theta_early <- as.numeric(theta_matrices_nilo[[3]]!=0)

for(j in 2:length(avg_abs_beta_early)){
  
  rownames(theta_matrices_nilo[[3 + (j-1)*6]]) <- names(beta_matrices_nilo[[3 + (j-1)*6]])
  
  sum_theta_early <- sum_theta_early + as.numeric(theta_matrices_nilo[[3 + (j-1)*6]]!=0)
}

rownames(theta_matrices_nilo[[4]]) <- names(beta_matrices_nilo[[4]])
sum_theta_late <- as.numeric(theta_matrices_nilo[[4]]!=0)

for(j in 2:length(avg_abs_beta_late)){
  
  rownames(theta_matrices_nilo[[4 + (j-1)*6]]) <- names(beta_matrices_nilo[[4 + (j-1)*6]])
  
  sum_theta_late <- sum_theta_late + as.numeric(theta_matrices_nilo[[4 + (j-1)*6]]!=0)
}


rownames(theta_matrices_nilo[[5]]) <- names(beta_matrices_nilo[[5]])
sum_theta_coop <- as.numeric(theta_matrices_nilo[[5]]!=0)

for(j in 2:length(avg_abs_beta_coop_nilo)){

  rownames(theta_matrices_nilo[[5 + (j-1)*6]]) <- names(beta_matrices_nilo[[5 + (j-1)*6]])

  sum_theta_coop <- sum_theta_coop + as.numeric(theta_matrices_nilo[[5 + (j-1)*6]]!=0)
}

rownames(theta_matrices_nilo[[6]]) <- names(beta_matrices_nilo[[6]])
sum_theta_coop_adap <- as.numeric(theta_matrices_nilo[[6]]!=0)

for(j in 2:length(avg_abs_beta_coop_adap_nilo)){

  rownames(theta_matrices_nilo[[6 + (j-1)*6]]) <- names(beta_matrices_nilo[[6 + (j-1)*6]])

  sum_theta_coop_adap <- sum_theta_coop_adap + as.numeric(theta_matrices_nilo[[6 + (j-1)*6]]!=0)
}

beta_sel <- c(sum(sum_beta_x1>cutoff), sum(sum_beta_x2>cutoff), sum(sum_beta_early>cutoff), sum(sum_beta_late>cutoff), sum(sum_beta_coop_nilo>cutoff), sum(sum_beta_coop_adap_nilo>cutoff))
theta_sel <- c(sum(sum_theta_x1>cutoff), sum(sum_theta_x2>cutoff), sum(sum_theta_early>cutoff), sum(sum_theta_late>cutoff), sum(sum_theta_coop>cutoff), sum(sum_theta_coop_adap>cutoff))

results_table = data.frame(colMeans(test_err), apply(test_err,2,sd)/sqrt(10), beta_sel,theta_sel)
colnames(results_table) <- c('Test MSE', 'SD', 'Main effects', 'Interactions')
rownames(results_table) <- c('Only X1', 'Only X2', 'Early fusion', 'Late fusion', 'CoopPliable','AdapCoopPliable')
results_table




### chord diagram

library(circlize)
library(networkD3)
library(igraph)
library(reshape2)
library(sysfonts)
library(showtext)
library(caret)
library(pliable)

beta_counts_list = list()
beta_values_list = list()

for(j in 1:simN){
  beta_counts_list[[j]] = beta_matrices_nilo[[5 + (j-1)*6]]!=0
  beta_values_list[[j]] = beta_matrices_nilo[[5 + (j-1)*6]]
}

beta_counts = tapply(unlist(beta_counts_list), names(unlist(beta_counts_list)), sum)
beta_values = tapply(unlist(beta_values_list), names(unlist(beta_values_list)), sum)

#rownames(theta_matrices_nilo[[5]]) <- names(beta_matrices_nilo[[5]])
theta_counts_list = list()
theta_values_list = list()
#theta_counts <- as.numeric(theta_matrices_nilo[[5]]!=0)
#theta_values <- theta_matrices_nilo[[5]]

for(j in 1:simN){
  
  rownames(theta_matrices_nilo[[5 + (j-1)*6]]) <- names(beta_matrices_nilo[[5 + (j-1)*6]])
  
  theta_counts_list[[j]] <- theta_matrices_nilo[[5 + (j-1)*6]]!=0
  theta_values_list[[j]] <- theta_matrices_nilo[[5 + (j-1)*6]]
}

theta_counts = Reduce('+', theta_counts_list)
theta_values = Reduce('+', theta_values_list)

tsc1_blood_mut <- c()
for(j in 1:simN){
  tsc1_blood_mut <- append(tsc1_blood_mut, beta_values_list[[j]]['TSC1.MUT'] + theta_values_list[[j]]['TSC1.MUT',3])
}

tsc1_lung_mut <- c()
for(j in 1:simN){
  tsc1_lung_mut <- append(tsc1_lung_mut, beta_values_list[[j]]['TSC1.MUT'] + theta_values_list[[j]]['TSC1.MUT',9])
}

bcrabl_blood_mut <- c()
for(j in 1:simN){
  bcrabl_blood_mut <- append(bcrabl_blood_mut, beta_values_list[[j]]['BCR_ABL.MUT'] + theta_values_list[[j]]['BCR_ABL.MUT',3])
}

boxplot(as.numeric(tsc1_blood_mut), as.numeric(lapply(beta_values_list, FUN=function(x){x['TSC1.MUT']})), as.numeric(tsc1_lung_mut))
boxplot(as.numeric(bcrabl_blood_mut), as.numeric(lapply(beta_values_list, FUN=function(x){x['BCR_ABL.MUT']})))


plot.circ<-function(beta,theta,X,Z,adj = c(.15, .8),cex = .5){
  
  data_mat = merge(as.data.frame(beta), as.data.frame(theta), by='row.names')
  rownames(data_mat) <- data_mat$Row.names
  data_mat = subset(data_mat, select = -c(Row.names) )
  colnames(data_mat) <- c("main", colnames(Z) )
  new_data<-data_mat[order( data_mat$main,decreasing = T),]
  new_data <- new_data[new_data$main!=0,]
  new_data<-as.matrix(new_data)
  new_data <- new_data[,-1]
  
  df = data.frame(from = rep(rownames(new_data), times = ncol(new_data)),
                  to = rep(colnames(new_data), each = nrow(new_data)),
                  value = as.vector(new_data),
                  stringsAsFactors = FALSE)
  df<-df[df$value!=0,]
  
  
  chordDiagram(df,big.gap = 10,transparency = 0,scale = F, annotationTrack = "grid",# grid.col = c("red", "blue", "purple","pink", "orange", "green"),
               preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(df))))))
  
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
                facing = "clockwise", niceFacing = TRUE, adj = adj)#,cex = cex)
  }, bg.border = NA) # here set bg.border to NA is important
  
}





# get X and Z

GDSC<-readRDS("C:/Users/dalma/Desktop/Matteo/uni/magistrale/erasmus/thesis/GDSC_data_matteo.rds")

name_drug <- c("Nilotinib")

# extract the drugs' pharmacological profiling and tissue dummy
col_filter <- colnames(GDSC$y) %in% paste("IC50.", name_drug,sep="")
YX0 <- cbind(
  GDSC$y[, col_filter],
  GDSC$z[,c(1:12)]
)





X23 <- GDSC$x1
X1 <- log2(X23)

X2=GDSC$x2;y=GDSC$y[, col_filter];Z=GDSC$z[,-c(6)]
y_name<-stringr::str_remove(colnames(y), pattern = "IC50.")

drug_names <- stringr::str_remove(colnames(GDSC$y), pattern = "IC50.")

X2=GDSC$x2;Z=GDSC$z[,-c(6)]

rownames(X1)<-GDSC$names.cell.line
rownames(X2)<-GDSC$names.cell.line
rownames(Z)<-GDSC$names.cell.line

X = cbind(X1,X2)

# plotting

# beta_values and theta_values contain the average value
# of beta and theta coefficients from the GDSC analysis

# beta_counts and theta_counts contain the number of times
# a certain coefficients was included in the model

beta = abs(beta_values) * (beta_counts>cutoff)
theta = abs(theta_values)* (theta_counts>cutoff)



#plot.circ(beta,theta,X,Z)


## interaction scatterplot

library(RColorBrewer)
library(latex2exp)
library(ggplot2)
library(caret)
library(gridExtra)

tsc1mut <- X2[,'TSC1.MUT']
bcrablmut <- X2[,'BCR_ABL.MUT']
myvar1 = as.factor(bcrablmut + 2*tsc1mut)
levels(myvar1) <- c('Other', 'BCR-ABL', 'TSC1')
myvar2 = X1[,'DDEF2']

cancertype_name <- as.factor(unlist(lapply(1:nrow(Z), FUN = function(x){names(GDSC$z[x,][GDSC$z[x,]==1])})))
cancertype_blood <- as.factor(
  ifelse(cancertype_name == 'blood',
         'Blood','Other'))


#ggplot(mydata1, aes(x=cancertype_name, y=GDSC$y[,i]))

mydata1 = data.frame(myvar1,y,cancertype_blood)
mydata2 = data.frame(myvar2,y,cancertype_blood)

p1 = ggplot(mydata1, aes(x=cancertype_blood, y=y, col=myvar1))+ xlab('') + ylab('IC50 Nilotinib') +
   theme(plot.title = element_text(face='bold', size='22'),
         axis.text=element_text(size=26),
         axis.title=element_text(size=24),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.background = element_blank(),
         legend.title=element_text(size='22'),
         legend.text=element_text(size='20'),
         axis.line = element_line(colour = "black"),
         panel.border = element_rect(colour = "black", fill=NA),
         plot.margin = unit(c(0.5, 1, 0.25, 0.5), "cm")
         ) +
  geom_jitter(size=4, width=0.1) +
  labs(colour="Mutation") + 
  geom_point(x=cancertype_blood[321],y=y[321], col='#619CFF', size=4)

# Assuming tsc1_blood_mut, beta_values_list, and tsc1_lung_mut are your data sources
tsc1_blood_mut <- as.numeric(tsc1_blood_mut)
tsc1_lung_mut <- as.numeric(tsc1_lung_mut)

# Extract the 'TSC1.MUT' values from beta_values_list
beta_values_tsc1_mut <- as.numeric(sapply(beta_values_list, function(x) x['TSC1.MUT']))

# Combine the data into a data frame
data <- data.frame(
  Group = factor(rep(c("Blood", "Lung"), c(length(tsc1_blood_mut), length(tsc1_lung_mut)))),
  Value = c(tsc1_blood_mut, tsc1_lung_mut)
)

# Create the boxplot using ggplot2
p2 <- ggplot(data, aes(x = Group, y = Value)) +
  geom_boxplot(fill = '#619CFF') +
  labs(x = NULL, y = "TSC1 Coefficients") +
  theme_bw() +
theme(plot.title = element_text(face='bold', size='22'),
      axis.text=element_text(size=26),
      axis.title=element_text(size=24),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      legend.position = "none",
      axis.line = element_line(colour = "black"),
      panel.border = element_rect(colour = "black", fill=NA),
      plot.margin = unit(c(0.5, 0.5, 1, 1), "cm")
)


grid_plot <- grid.arrange(p1, p2, widths=c(1.2,1), nrow = 1)

ggsave("GDSC_scatterplots.eps", grid_plot, width = 13, height = 8, dpi=800)
ggsave("GDSC_scatterplots.pdf", grid_plot, width = 13, height = 8)


print(results_table)