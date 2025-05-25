################################################
################# PACKAGES #####################
################################################

library(readr)
library(tidyr)
library(plyr)
library(stringr)
library(forcats)

library(ggplot2)
library(gridExtra)
library(patchwork)

################################################
################# FUNCTIONS ####################
################################################


# Generate graphs of posterior distributions of differences between conditions
plot_diff <- function(input_fit, input_fit_name, parameter, conditions, density_color) {
  
  # Get the name of the experiment
  experiment <- paste0("(Experiment ", str_extract(input_fit_name, "[0-9]"), ")")
  
  
  # Get the posterior distribution of condition differences
  post_diff <- data.frame(as.matrix(input_fit, pars=parameter))
  names(post_diff) <- c("post_diff")
  
  # Get posterior density and credible interval
  density_diff <- with(density(post_diff$post_diff), data.frame(x, y))
  ci_diff <- bayestestR::ci(post_diff$post_diff, method="HDI")
  
  # Plot density with credible interval shaded
  ggplot(post_diff, aes(x = post_diff)) +
    geom_density(color = rgb(0,0,0,0), fill = density_color, size = 1, alpha = .2) +
    scale_x_continuous(name = paste("RT Advantage for Incidental Exposure")) +
    scale_y_continuous(name = paste(conditions, experiment, sep = " ")) +
    geom_area(data =  density_diff[density_diff$x>ci_diff$CI_low & density_diff$x<ci_diff$CI_high,], aes(x = x, y = y), fill = density_color, alpha = .5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = .9) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.title = element_text(size=12, face="bold"),
          axis.text = element_text(size=8),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
}



# Generate graphs of the RT data
plot_data <- function(input_data, input_data_name, outcome, rt_limits, incidental_color) {
  
  # Get the name of the experiment
  experiment <- paste0("Experiment ", str_extract(input_data_name, "[0-9]"))
  
  
  ggplot(input_data[input_data$Accuracy == 1,], aes_string(x = "Block", y = outcome)) +
    stat_summary(aes(color = Condition), geom = "line", fun = "mean", position = position_dodge(width = .1), show.legend=F) +
    stat_summary(aes(group = Condition), geom = "linerange", fun.data = "mean_se", position = position_dodge(width = .1), alpha = .8, color = "black", show.legend=F) +
    scale_x_continuous(name = "Exposure Phase Block", breaks = 1:8) +
    scale_y_continuous(name = paste(experiment, "\nRT on Accurate Trials")) +
    scale_color_manual(values = c(incidental_color, "grey")) +
    scale_fill_manual(values = c(incidental_color, "grey")) +
    coord_cartesian(ylim = rt_limits) +
    theme_bw() +
    theme(panel.grid=element_blank(),
          axis.title=element_text(size=12, face="bold"),
          axis.text=element_text(size=8))
}


################################################
############ LOAD AND FORMAT DATA ##############
################################################

# This section reads in the data and the model fits


# Define file paths and object names
exp_names <- c("1A", "1B", "2", "3")
file_names <- paste0("data/exposure/exposure_", exp_names, ".rds")
object_names <- c("exp1A", "exp1B", "exp2", "exp3")

# Read all datasets into a named list
exp_list <- file_names %>%
  set_names(object_names) %>%
  map(readRDS)

# Apply the same RT filter used in analysis
exp_list$exp3 <- exp_list$exp3 %>%
  group_by(ID) %>%
  dplyr::mutate(RT_filter = ifelse(RT > 2000 | Accuracy == 0, NA, RT))


# Identify the winning models from the analysis and 
# read the fits into a list
fit_file_names <- list.files(path = "stan_fits/exposure/", full.names = T)
fit_object_names <- paste0("fit_", 
                           str_extract(fit_file_names, pattern = "acc|rt"), "_exp",
                           str_extract(fit_file_names, pattern = "[0-9]+[AB]*"))

fit_list <- fit_file_names %>%
  set_names(fit_object_names) %>%
  map(readRDS)


# Reorder by dependent variable and experiment
fit_ID <- gsub("fit_", "", fit_object_names)
desired_order <- c("rt_exp1A", "rt_exp1B", "rt_exp2", "rt_exp3",
                   "acc_exp1A", "acc_exp1B", "acc_exp2")  # Expand this if needed

# Step 5: Reorder list by custom experiment order
fit_list <- fit_list[order(factor(fit_ID, levels = desired_order))]

######################################################
############### CONDITION DIFFERENCES ################
######################################################

# Compute HDIs of condition differences and return a tibble with model name
hdi_summary <- fit_list %>%
  imap_dfr(~ {
    samples <- as.matrix(.x) %>%
      `[`(, grep("slope_diff", colnames(.))) %>%  # keep only columns matching "slope_diff"
      as.data.frame()
    
    bayestestR::ci(samples, method = "HDI", ci = 0.89) %>%
      mutate(model = .y)
  })


hdi_summary

######################################################
######## GRAPHS: MODEL CONDITION DIFFERENCES #########
######################################################

# Plot distributions of condition differences with credible intervals highlighted

# Use plot_diff function defined in functions block above.
# First define the arguments to be used in the function for each experiment
plot_diff_meta <- tribble(
  ~fit_index, ~plot_name,         ~parameter,         ~conditions,   ~density_color,  
  1,          "plot_diff_exp1A",  "post_slope_diff",   "Base vs Dense",  rgb(0, 0.1, 1), 
  2,          "plot_diff_exp1B",  "post_slope_diff",   "Base vs Sparse", rgb(0.8, 0, 0.1), 
  3,          "plot_diff_exp2",   "post_slope_diff",   "Base vs Dense",  rgb(0, 0.1, 1),     
  4,          "plot_diff_exp3",   "slope_diff",        "Base vs Dense",  rgb(0, 0.1, 1),     
)


# Create named list of plots
plot_diff_list <- pmap(
  plot_diff_meta,
  function(fit_index, plot_name, parameter, conditions, density_color) {
    fit <- fit_list[[fit_index]]
    fit_name <- names(fit_list)[fit_index]
    plot_diff(fit, fit_name, parameter, conditions, density_color)
  }
) %>%
  set_names(plot_diff_meta$plot_name)


grid.arrange(
  grobs = plot_diff_list,
  layout_matrix = rbind(c(1,2),
                        c(3,4))
)

######################################################
################### GRAPHS: DATA #####################
######################################################


# Just the data - RT on accurate trials across trials

plot_data_meta <- tribble(
  ~exp_index, ~plot_name,         ~outcome,         ~rt_limits,   ~incidental_color,  
  1,          "plot_data_exp1A",  "RT",   c(230, 380),  rgb(0, 0.1, 1), 
  2,          "plot_data_exp1B",  "RT",   c(230, 380), rgb(0.8, 0, 0.1), 
  3,          "plot_data_exp2",   "RT",   c(230, 380),  rgb(0, 0.1, 1),     
  4,          "plot_data_exp3",   "RT_filter",         c(400, 900),  rgb(0, 0.1, 1),     
)


# Create named list of plots
plot_data_list <- pmap(
  plot_data_meta,
  function(exp_index, plot_name, outcome, rt_limits, incidental_color) {
    exp <- exp_list[[exp_index]]
    exp_name <- names(exp_list)[exp_index]
    plot_data(exp, exp_name, outcome, rt_limits, incidental_color)
  }
) %>%
  set_names(plot_data_meta$plot_name)

grid.arrange(
  grobs = plot_data_list,
  layout_matrix = rbind(c(1,2),
                        c(3,4))
)

