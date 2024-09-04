library(tidyr)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(grid)
library(png)


## Script to generate the results for low and high dimensional simulation studies

## Figures 6 and 7


num_sims = 4

filenames = c(
  "100_ntrain500_ntest9800_pimp30_px1100_px2100_tx12_tx22_sigma33_factstr4_SNR1_with_adap.csv",
  "100_ntrain500_ntest9800_pimp30_px1100_px2100_tx14_tx24_sigma40_factstr4_SNR1_with_adap.csv",
  "100_ntrain500_ntest9800_pimp30_px1100_px2100_tx10_tx20_sigma22_factstr4_SNR2_with_adap.csv",
  "100_ntrain500_ntest9800_pimp30_px1100_px2100_tx10_tx20_sigma24_factstr8_SNR3_with_adap.csv"
)

data_df <- read.csv(filename)[,-1]
fills <- c('#00BFC4', '#F8766D',rep('white', 4))
rholist = c(0,0.2,0.4,0.6,0.8,1,3,5,9)
limits <- list(c(NA,NA),c(NA,NA),c(NA,NA),c(NA,225))

all_plots <- list()
row_labels <- list()

for(i in 1:length(filenames)){
  
  filename = filenames[i]
  data_df <- read.csv(filename)[,-1]
  
  data_mse_long <- gather(data_df[,1:6], method, MSE, obj_x1s:obj_coops_adap, factor_key=TRUE)
  levels(data_mse_long$method) <- c('Only X1', 'Only X2', 'Early fusion', 'Late fusion', 'CoopPliable', 'Adap CoopPliable')
  data_mse_long$method <- factor(data_mse_long$method, levels = rev(unique(data_mse_long$method)))
  
  
  plot_mse <- ggplot(data_mse_long, aes(x=method, y=MSE, fill=method)) + 
    geom_boxplot(show.legend = FALSE) + scale_x_discrete(name=NULL)+ scale_y_continuous(limits=limits[[i]])+
    coord_flip() + ylab('Test MSE') + scale_fill_manual(values=fills)+
    theme_bw(base_size = 15)+ theme(axis.text.y = element_text(face = "bold", size=15), axis.title.x = element_text(face = "bold"), panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank(), plot.margin=unit(c(6.5,5.5,11,5.5), 'pt'))
  
  
  data_sparsity1_long <- gather(data_df[,7:12], method, feature_selected, beta_sel_x1s:beta_sel_coops_adap, factor_key=TRUE)
  levels(data_sparsity1_long$method) <- c('Only X1', 'Only X2', 'Early fusion', 'Late fusion', 'CoopPliable', 'Adap CoopPliable')
  data_sparsity1_long$method <- factor(data_sparsity1_long$method, levels = rev(unique(data_sparsity1_long$method)))
  
  plot_sparsity1 <- ggplot(data_sparsity1_long, aes(x=method, y=feature_selected, fill=method)) + 
    geom_boxplot(show.legend = FALSE) + scale_x_discrete(name=NULL, labels = NULL)+
    coord_flip() + ylab('Number of features') +scale_fill_manual(values=fills)+
    theme_bw(base_size = 15)+ theme(plot.margin=unit(c(6.5,5.5,11,5.5), 'pt'), panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank(), axis.title.x = element_text(face = "bold"))
  
  data_sparsity2_long <- gather(data_df[,13:18], method, feature_selected, theta_sel_x1s:theta_sel_coops_adap, factor_key=TRUE)
  levels(data_sparsity2_long$method) <- c('Only X1', 'Only X2', 'Early fusion', 'Late fusion', 'CoopPliable', 'Adap CoopPliable')
  data_sparsity2_long$method <- factor(data_sparsity2_long$method, levels = rev(unique(data_sparsity2_long$method)))
  
  plot_sparsity2 <- ggplot(data_sparsity2_long, aes(x=method, y=feature_selected, fill=method)) + 
    geom_boxplot(show.legend = FALSE) + scale_x_discrete(name=NULL, labels = NULL)+
    coord_flip() + ylab('Number of interactions') +scale_fill_manual(values=fills)+
    theme_bw(base_size = 15)+ theme(plot.margin=unit(c(6.5,5.5,11,5.5), 'pt'), panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank(), axis.title.x = element_text(face = "bold"))
  
  
  counts = sapply(rholist, function(x) sum(data_df$chosen_rhos == x))
  counts_adap = sapply(rholist, function(x) sum(data_df$chosen_rhos_adap == x))
  rho_df <- data.frame(rholist, counts, counts_adap)
  rho_df <- gather(rho_df, method, counts, counts:counts_adap, factor_key=TRUE)
  levels(rho_df$method) <- c('CoopPliable', 'Adap CoopPliable')
  colnames(rho_df)[1] <- 'rho'
  rho_df$rho <- as.factor(rho_df$rho)
  plot_rho <- ggplot(rho_df) + 
    geom_col(aes(x=rho, y=counts, fill=method), colour="black", position="dodge", show.legend = FALSE) + labs(x=as.expression(bquote(bold(paste('Selected ', rho))))) + ylim(c(0,10))+ 
    theme_bw(base_size = 15) + theme(axis.title.y = element_blank(), axis.ticks = element_blank() , panel.grid.major = element_blank(),
                                     panel.grid.minor = element_blank(),legend.title = element_blank(), 
                                     legend.position=c(.65,.8), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  all_plots <- c(all_plots, list(plot_mse, plot_sparsity1, plot_sparsity2, plot_rho))
  row_labels <- c(row_labels, list(textGrob(as.character(i), gp=gpar(fontsize=20, fontface="bold"),vjust = -6.5)))
  
}

grid_plot <- grid.arrange(
  arrangeGrob(textGrob(" "), grobs = row_labels, ncol = 1),  # Empty space for top left corner
  arrangeGrob(grobs = all_plots, ncol = 4, nrow = 4, widths=c(1.07,0.7, 0.7,0.7)),
  widths = c(0.3, 15),  # Adjust width to fit labels
  nrow = 1
)

ggsave("lowdim_wholegrid_withlabels.eps", grid_plot, width = 15, height = 13, dpi=800)
ggsave("lowdim_wholegrid_withlabels.pdf", grid_plot, width = 15, height = 13)

