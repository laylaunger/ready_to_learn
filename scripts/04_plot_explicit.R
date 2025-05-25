# Title: Explicit Phase Analysis
# Author: Layla Unger
# R version: 4.3.2
# Packages: readr, stringr, dplyr, tidyr, purrr, ggplot2, patchwork, gridExtra, ggpubr


################################################
################ DESCRIPTION ###################
################################################

# This script is generates graphs to visualize the results of the
# Explicit phase. To run this script, you need to have fit the 
# models and stored their results from the script 02_analysis_explicit.R

################################################
################# PACKAGES #####################
################################################

library(readr)
library(tidyr)
library(plyr)
library(stringr)

library(ggplot2)
library(gridExtra)
library(patchwork)

################################################
################# FUNCTIONS ####################
################################################

# Generate graphs of posterior distributions of differences between conditions
# from a fit object. Only the winning models are loaded below, so the input_fit
# object will always be the winning model for an experiment. This function
# checks whether the winning model was one with condition differences in
# the slope or intercept, and plots those differences if so.
plot_diff <- function(input_fit, input_fit_name, conditions, density_color, exp_contrast = NULL) {
  
  # Nothing to graph if the input_fit is the baseline model (no condition differences)
  if(grepl("base", input_fit_name)) {
    print("No diff to plot")
    return()
  }
  
  # Get the name of the experiment
  experiment <- paste("Experiment", str_extract(input_fit_name, "[0-9]$"))
  
  # If the winning model had condition differences in slope, plot those;
  # otherwise, plot intercept
  type_diff <- ifelse(grepl("slope", input_fit_name), "(slope)", "(intercept)")
  
  # Identify the parameter that differed across conditions
  parameter <- ifelse(type_diff == "(slope)", "post_slope_diff", "post_alpha_diff")
  
  # For experiment 1, there was a contrast between both Base & Dense and Base & Sparse
  # The exp_contrast argument focuses on Base vs Dense by default, but can be set to 2
  # for experiment 1 to focus on Base vs Sparse
  if(!is.null(exp_contrast) && !is.na(exp_contrast)) {
    parameter <- paste0(parameter, "[", exp_contrast, "]")
  }
  
  post_diff <- data.frame(as.matrix(input_fit, pars = parameter))
  
  names(post_diff) <- "post_diff"
  
  density_diff <- with(density(post_diff$post_diff), data.frame(x, y))
  
  ci_diff <- bayestestR::ci(post_diff$post_diff, method = "HDI")
  
  ggplot(post_diff, aes(x = post_diff)) +
    geom_area(
      data = density_diff[density_diff$x > ci_diff$CI_low & density_diff$x < ci_diff$CI_high,],
      aes(x = x, y = y),
      fill = density_color, alpha = .8
    ) +
    geom_density(color = "black", fill = density_color, size = 1, alpha = .1) +
    scale_x_continuous(name = paste("Category Learning Advantage for\nCategory Exposure", type_diff)) +
    scale_y_continuous(name = experiment) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = .9) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 8),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
}


# Generate graphs of data
plot_data <- function(input_data, input_data_name, incidental_color, condition_contrast = "Dense") {
  
  # Get the name of the experiment
  experiment <- paste("Experiment", str_extract(input_data_name, "[0-9]"))
  
  # Focus on the condition contrast of interest
  input_data <- input_data %>%
    dplyr::filter(grepl(condition_contrast, Condition))
  
  ggplot(input_data, aes(x = factor(Block), y = Accuracy)) +
    stat_summary(
      aes(group = Condition, color = Condition),
      geom = "line", fun = "mean", position = position_dodge(width = .1), show.legend = FALSE
    ) +
    stat_summary(
      aes(group = Condition),
      geom = "linerange", fun.data = "mean_se", position = position_dodge(width = .1)
    ) +
    stat_summary(
      aes(group = Condition, color = Condition),
      geom = "point", fun = "mean", position = position_dodge(width = .1), show.legend = FALSE
    ) +
    scale_x_discrete(name = "Block") +
    scale_y_continuous(
      name = paste(experiment, "Accuracy", sep = "\n"),
      breaks = seq(.4, 1, by = .2)
    ) +
    scale_color_manual(values = c(incidental_color, "grey")) +
    coord_cartesian(ylim = c(.4, 1)) +
    geom_hline(yintercept = .5, linetype = "dashed") +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 8)
    )
}


################################################
############ LOAD AND FORMAT DATA ##############
################################################

# This section reads in the data and the model fits
# for just the winning models.

# Define file paths and object names for data
exp_file_names <- paste0("data/explicit/explicit_", 1:5, ".rds")
exp_object_names <- paste0("exp", 1:5)

# Read all datasets into a named list
exp_list <- exp_file_names %>%
  set_names(exp_object_names) %>%
  map(readRDS)


# Identify the winning models from the analysis and 
# read the fits into a list
exp_winners <- read_csv(file = "stan_fits/explicit/explicit_winners.csv", show_col_types = FALSE)

fit_file_names <- paste0("stan_fits/explicit/fit_", exp_winners$winner, "_exp", exp_winners$exp, ".rds")
fit_object_names <- paste0("fit_", exp_winners$winner, "_exp", exp_winners$exp)

fit_list <- fit_file_names %>%
  set_names(fit_object_names) %>%
  map(readRDS)


######################################################
############# GRAPHS: LEGEND FOR PLOTS ###############
######################################################

# Create a legend to accompany the plots for all experiments. Use exp1 data because this
# includes all the conditions that were used across all experiments so that the legend
# will include a key for all conditions
legend_data <- exp_list$exp1
legend_data$Condition <- factor(legend_data$Condition, levels = c("Dense", "Sparse", "Base Dense", "Base Sparse"))
legend <- ggplot(legend_data, aes(x = Condition, y = Accuracy)) +
  stat_summary(aes(fill = Condition), geom = "bar") +
  scale_fill_manual(breaks = c("Dense", "Sparse"), 
                    values = rgb(c(0, .8, 0, 0), c(.1, 0, 0, 0), c(1, .1, 0, 0)),
                    name = "", labels = c("Dense Structure", "Sparse Structure")) +
  theme(legend.position="bottom",
        legend.direction = "horizontal") +
  guides(fill = guide_legend(title.position="top", title.hjust = 0.5))
legend <- as_ggplot(get_legend(legend))

######################################################
######## GRAPHS: MODEL CONDITION DIFFERENCES #########
######################################################

# Plot distributions of condition differences with credible intervals highlighted

# Use plot_diff function defined in functions block above.
# First define the arguments to be used in the function for each experiment
plot_diff_meta <- tribble(
  ~fit_index, ~plot_name,         ~conditions,            ~density_color,         ~exp_contrast,
  1,          "plot_diff_exp1A",  "Dense vs Baseline",    rgb(0, 0.1, 1),         1,
  1,          "plot_diff_exp1B",  "Sparse vs Baseline",   rgb(0.8, 0, 0.1),       2,
  2,          "plot_diff_exp2",   "Dense vs Baseline",    rgb(0, 0.1, 1),         NA,
  3,          "plot_diff_exp3",   "Dense vs Baseline",    rgb(0, 0.1, 1),         NA,
  4,          "plot_diff_exp4",   "Dense vs Baseline",    rgb(0, 0.1, 1),         NA,
  5,          "plot_diff_exp5",   "Dense vs Baseline",    rgb(0, 0.1, 1),         NA
)

# Create named list of plots
plot_diff_list <- pmap(
  plot_diff_meta,
  function(fit_index, plot_name, conditions, density_color, exp_contrast) {
    fit <- fit_list[[fit_index]]
    fit_name <- names(fit_list)[fit_index]
    plot_diff(fit, fit_name, conditions, density_color, exp_contrast)
  }
) %>%
  set_names(plot_diff_meta$plot_name)


# Now plot! 

# Combine all plots using wrap_plots (from patchwork)
plot_diff_grid <- wrap_plots(plot_diff_list, ncol = 2)

# Combine legend and plot grid, with legend on top
plot_diff_all <- legend / plot_diff_grid + plot_layout(heights = c(0.15, 1))

# Display it
print(plot_diff_all)

######################################################
################### GRAPHS: DATA #####################
######################################################

# Just the data, binned into blocks

# Use plot_data function defined in functions block above.
# First define the arguments to be used in the function for each experiment
plot_data_meta <- tribble(
  ~data_index, ~plot_name,                   ~incidental_color,  ~condition_contrast,      
  1,          "plot_data_exp1A",     rgb(0, 0.1, 1),         "Dense",
  1,          "plot_data_exp1B",     rgb(0.8, 0, 0.1),       "Sparse",
  2,          "plot_data_exp2",      rgb(0, 0.1, 1),         "Dense",
  3,          "plot_data_exp3",      rgb(0, 0.1, 1),         "Dense",
  4,          "plot_data_exp4",      rgb(0, 0.1, 1),         "Dense",
  5,          "plot_data_exp5",      rgb(0, 0.1, 1),         "Dense",
)

# Create named list of plots
plot_data_list <- pmap(
  plot_data_meta,
  function(data_index, plot_name, incidental_color, condition_contrast) {
    exp <- exp_list[[data_index]]
    exp_name <- names(exp_list)[data_index]
    plot_data(exp, exp_name, incidental_color, condition_contrast)
  }
) %>%
  set_names(plot_data_meta$plot_name)


# Now plot! 

# Combine all plots using wrap_plots (from patchwork)
plot_data_grid <- wrap_plots(plot_data_list, ncol = 2)

# Combine legend and plot grid, with legend on top
plot_data_all <- legend / plot_data_grid + plot_layout(heights = c(0.15, 1))

# Display it
print(plot_data_all)
