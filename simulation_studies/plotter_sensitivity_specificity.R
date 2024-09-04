library(ggplot2)
library(gridExtra)
library(cowplot)
library(reshape2)


setwd("C:/Users/dalma/Desktop/Matteo/uni/magistrale/erasmus/thesis")

p_imp = 30
p1 = 500
p2 = 500


# Function to generate sensitivity and specificity data
generate_data <- function(simseed, sigma, tx1, tx2, factor_strength, snr_avg) {
  
  # Initialize the parameters
  real_main_effects <- c(rep(TRUE,p_imp),rep(FALSE,p1-p_imp),rep(TRUE,p_imp),rep(FALSE,p2-p_imp))
  if(factor_strength==8) real_main_effects <- c(rep(TRUE,p_imp),rep(FALSE,p1-p_imp),rep(FALSE,p2))

  
  sim_filename = paste(simseed, paste0("ntrain", n_train), paste0("ntest", n_test),
                       paste0("pimp", p_imp),
                       paste0("px1", p1), paste0("px2", p2),
                       paste0("tx1", tx1), paste0("tx2", tx2),
                       paste0("sigma", sigma),
                       paste0("factstr", factor_strength),
                       paste0("SNR", as.integer(snr_avg)),
                       sep = "_")
  
  beta_x1_pres <- readRDS(paste0("beta_x1_pres_",sim_filename,".RDS"))
  beta_x2_pres <- readRDS(paste0("beta_x2_pres_",sim_filename,".RDS"))
  beta_early_pres <- readRDS(paste0("beta_early_pres_",sim_filename,".RDS"))
  beta_late_pres <- readRDS(paste0("beta_late_pres_",sim_filename,".RDS"))
  beta_coop_pres <- readRDS(paste0("beta_coop_pres_",sim_filename,".RDS"))
  beta_coop_adap_pres <- readRDS(paste0("beta_coop_adap_pres_",sim_filename,".RDS"))
  
  sens_x1 <- rep(0, 10)
  sens_x2 <- rep(0, 10)
  sens_late <- rep(0, 10)
  sens_early <- rep(0, 10)
  sens_coop <- rep(0, 10)
  sens_coop_adap <- rep(0, 10)
  spec_x1 <- rep(0, 10)
  spec_x2 <- rep(0, 10)
  spec_early <- rep(0, 10)
  spec_late <- rep(0, 10)
  spec_coop <- rep(0, 10)
  spec_coop_adap <- rep(0, 10)
  
  for (cutoff in 1:10) {
    sens_x1[cutoff] <- sum(real_main_effects != 0 & beta_x1_pres >= cutoff) / sum(real_main_effects != 0)
    sens_x2[cutoff] <- sum(real_main_effects != 0 & beta_x2_pres >= cutoff) / sum(real_main_effects != 0)
    sens_early[cutoff] <- sum(real_main_effects != 0 & beta_early_pres >= cutoff) / sum(real_main_effects != 0)
    sens_late[cutoff] <- sum(real_main_effects != 0 & beta_late_pres >= cutoff) / sum(real_main_effects != 0)
    sens_coop[cutoff] <- sum(real_main_effects != 0 & beta_coop_pres >= cutoff) / sum(real_main_effects != 0)
    sens_coop_adap[cutoff] <- sum(real_main_effects != 0 & beta_coop_adap_pres >= cutoff) / sum(real_main_effects != 0)
    
    spec_x1[cutoff] <- sum(real_main_effects == 0 & beta_x1_pres < cutoff) / sum(real_main_effects == 0)
    spec_x2[cutoff] <- sum(real_main_effects == 0 & beta_x2_pres < cutoff) / sum(real_main_effects == 0)
    spec_early[cutoff] <- sum(real_main_effects == 0 & beta_early_pres < cutoff) / sum(real_main_effects == 0)
    spec_late[cutoff] <- sum(real_main_effects == 0 & beta_late_pres < cutoff) / sum(real_main_effects == 0)
    spec_coop[cutoff] <- sum(real_main_effects == 0 & beta_coop_pres < cutoff) / sum(real_main_effects == 0)
    spec_coop_adap[cutoff] <- sum(real_main_effects == 0 & beta_coop_adap_pres < cutoff) / sum(real_main_effects == 0)
  }
  
  sens_mat <- data.frame(Cutoff = 1:10, Only_X1 = sens_x1, Only_X2 = sens_x2, Early_fusion = sens_early, Late_fusion = sens_late, Coop = sens_coop, AdapCoop = sens_coop_adap)
  spec_mat <- data.frame(Cutoff = 1:10, Only_X1 = spec_x1, Only_X2 = spec_x2, Early_fusion = spec_early, Late_fusion = spec_late, Coop = spec_coop, AdapCoop = spec_coop_adap)
  
  sens_mat <- melt(sens_mat, id.vars = "Cutoff")
  spec_mat <- melt(spec_mat, id.vars = "Cutoff")
  
  levels(sens_mat$variable) <- c("Only X1", "Only X2", "Early fusion", "Late fusion", "CoopPliable", "Adap CoopPliable")
  levels(spec_mat$variable) <- c("Only X1", "Only X2", "Early fusion", "Late fusion", "CoopPliable", "Adap CoopPliable")
  
  return(list(sens_mat = sens_mat, spec_mat = spec_mat))
}

# Create plots for different simulations
plots <- list()
sim_params <- list(
  list(simseed = 100, sigma = 33, tx1 = 2, tx2 =2, factor_strength = 4, snr_avg = 1),
  list(simseed = 100, sigma = 40, tx1 = 4, tx2 =4, factor_strength = 4, snr_avg = 1),
  list(simseed = 100, sigma = 22, tx1 = 0, tx2 =0, factor_strength = 4, snr_avg = 2),
  list(simseed = 100, sigma = 24, tx1 = 0, tx2 =0, factor_strength = 8, snr_avg = 3)
)


for (i in 1:length(sim_params)) {
  params <- sim_params[[i]]
  data <- generate_data(params$simseed, params$sigma, params$tx1, params$tx2,  params$factor_strength, params$snr_avg)
  
  spec_mat <- data$spec_mat
  sens_mat <- data$sens_mat
  
  # Create the specificity plot
  spec_plot <- ggplot(spec_mat, aes(x = Cutoff, y = value, color = variable)) +
    geom_line(size = 1) + # Increase line thickness
    scale_x_continuous(breaks = seq(floor(min(spec_mat$Cutoff)), ceiling(max(spec_mat$Cutoff)), by = 1)) +
    labs(x = "Cutoff", y = "Specificity", color = "Method") + # Set axis labels and legend title
    theme_bw(base_size = 12)+ theme(plot.margin=unit(c(6.5,5.5,11,5.5), 'pt'), panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank(),  # Axis text
                                         legend.title = element_text(face='bold'), # Increase legend title size
                                        # legend.text = element_text(size = 35), # Increase legend text size
                                        # legend.key.size = unit(2, 'cm'),
                                        # legend.key.height = unit(0.5, 'cm'),
                                        legend.position = "bottom" # Position legend at the bottom
    )

  
  # Create the sensitivity plot without the legend
  sens_plot <- ggplot(sens_mat, aes(x = Cutoff, y = value, color = variable)) +
    geom_line(size = 1) + # Increase line thickness
    scale_x_continuous(breaks = seq(floor(min(sens_mat$Cutoff)), ceiling(max(sens_mat$Cutoff)), by = 1)) +
    labs(x = "Cutoff", y = "Sensitivity") +
    theme_bw(base_size = 12)+ theme(plot.margin=unit(c(6.5,5.5,11,5.5), 'pt'), panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank())# Set axis labels

  
  # Add plots to list
  plots[[paste0("spec_plot_", i)]] <- spec_plot
  plots[[paste0("sens_plot_", i)]] <- sens_plot
}

# Extract the legend from the first specificity plot
legend <- get_legend(plots$spec_plot_1)

# Remove the legend from the first specificity plot
plots$spec_plot_1 <- plots$spec_plot_1 + theme(legend.position = "none")
plots$spec_plot_2 <- plots$spec_plot_2 + theme(legend.position = "none")
plots$spec_plot_3 <- plots$spec_plot_3 + theme(legend.position = "none")
plots$spec_plot_4 <- plots$spec_plot_4 + theme(legend.position = "none")

# Remove the legend from the first specificity plot
plots$sens_plot_1 <- plots$sens_plot_1 + theme(legend.position = "none")
plots$sens_plot_2 <- plots$sens_plot_2 + theme(legend.position = "none")
plots$sens_plot_3 <- plots$sens_plot_3 + theme(legend.position = "none")
plots$sens_plot_4 <- plots$sens_plot_4 + theme(legend.position = "none")

# Combine the plots into a 4x2 grid with numbers 1, 2, 3, 4 at the top left of each row
combined_plot <- plot_grid(
  plot_grid(plots$spec_plot_1, plots$sens_plot_1, ncol = 2) + draw_label("1", x = 0, y = 0.99, hjust = 0, vjust = 1, size = 15, fontface="bold"),
  plot_grid(plots$spec_plot_2, plots$sens_plot_2, ncol = 2) + draw_label("2", x = 0, y = 0.99, hjust = 0, vjust = 1, size = 15, fontface="bold"),
  plot_grid(plots$spec_plot_3, plots$sens_plot_3, ncol = 2) + draw_label("3", x = 0, y = 0.99, hjust = 0, vjust = 1, size = 15, fontface="bold"),
  plot_grid(plots$spec_plot_4, plots$sens_plot_4, ncol = 2) + draw_label("4", x = 0, y = 0.99, hjust = 0, vjust = 1, size = 15, fontface="bold"),
  ncol = 1, rel_heights = c(1, 1, 1, 1)
)

# Add the legend below the combined plot
combined_plot_with_legend <- plot_grid(combined_plot, legend, ncol = 1, rel_heights = c(1, 0.1))

# ggsave("combined_sens_spec_plots.pdf" ,combined_plot_with_legend, , width = 7, height = 10)
# ggsave("combined_sens_spec_plots.eps" ,combined_plot_with_legend, , width = 7, height = 10, dpi=800)
# # Save the combined plot to a PDF file
# pdf("combined_sens_spec_plots.pdf", width = 7, height = 10) # Adjust width and height as needed
# print(combined_plot_with_legend)
# dev.off()







