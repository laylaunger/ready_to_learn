# Title: Explicit Phase Analysis
# Author: Layla Unger
# R version: 4.3.2
# Packages: readr, stringr, dplyr, tidyr, purrr, ggplot2, cowplot, gridExtra, ggpubr, rstan, bayestestR, ggmcmc


################################################
################ DESCRIPTION ###################
################################################

# This script is part of a project that investigated whether mere passive exposure
# to new categories fostered subsequent explicit learning.

# This project included a set of 5 experiments that each involved two phases: an
# Exposure phase involving passive exposure, and an Explicit phase that tested
# subsequent learning. This script analyzes data from the Explicit phase.

# In the explicit phase, participants learned two categories, "Flurps" and "Jalets", 
# with conventional supervision (i.e., categorized an exemplar, then received 
# corrective feedback after each trial)
# In each experiment, the Explicit Phase occurred after an "Exposure Phase" in 
# one of two conditions: (1) Baseline, or  (2) Incidental exposure to categories.
# In Experiment 1, categories were either Dense or Sparse. In all subsequent experiments, 
# categorieswere dense

# Goal of analyses: Test how condition affects category learning across trials
# Specifically, test hypothesis that Incidental exposure to dense categories facilitates 
# subsequent category learning in comparison to Baseline, whereas exposure to sparse 
# categories does not

# Data for each experiment are stored in a long-
# format dataframe. Important columns are:
# ID
# Condition
# Accuracy (0 or 1)
# Running_Trial

# Analysis steps currently in script:
# Steps currently in script:
# (1) Define a set of new functions that will be used in the script
# (2) Load a set of hierarchical Bayesian models that will be used to analyze the data.
#     In these models,the outcome variable (categorization Accuracy) of a logistic regression 
#     is predicted by trial number. 
#     Across the models, the intercept and slope for each participant are each taken 
#     from either a single distribution, or one of two distributions based on condition. When 
#     a parameter is taken from one of two distributions based on condition, for each iteration, 
#     the model compares the Baseline and Incidental slope distributions. This comparison is used 
#     to compare whether changes Accuracy across trials vary across the two conditions.
# (3) Load and format Accuracy data from each experiment.
# (4) Fit the models to the data
# (5) For each experiment, compare the fits of the different models using leave-one-out 
#     cross-validation, and identify the best-fitting model
# (6) From the best-fitting models, extract differences between the distributions for slopes in 
#     the two conditions for each sample from the posterior. 
# (7) Evaluate differences between distributions between the two conditions for each
#     experiment.


################################################
################# PACKAGES #####################
################################################

library(tidyr)
library(plyr)
library(stringr)

library(ggplot2)
library(cowplot)
library(gridExtra)
library(ggpubr)

library(rstan)
library(bayestestR)
library(loo)

rstan_options(auto_write = TRUE)
options(mc.cores = 3)

################################################
################# FUNCTIONS ####################
################################################

# This section defines functions that will be used in the script.

############## FORMATTING DATA #################

# Stan models need data in a list format. 

# Generate object in which acc data used to fit models is in list
# format used as input for stan models
gen_explicit_data <- function(experiment, num_comparisons) {
  
  #Identify each participant's condition
  sub_condition <- experiment %>%
    group_by(ID) %>%
    dplyr::summarise(Cond = unique(Condition_num)) %>%
    pull(Cond)
  
  dat <- list(num_sub = length(unique(experiment$ID_num)), #number of participants
              num_trials = length(unique(experiment$Running_Trial)), #number of trials for each participant
              num_conditions =length(unique(experiment$Condition_num)), #number of conditions
              num_comparisons = num_comparisons, #number of comparisons between conditions (base vs dense only, or also base vs sparse)
              sub_condition = sub_condition,
              trial = matrix(experiment$Running_Trial, nrow=length(unique(experiment$ID)), ncol=length(unique(experiment$Running_Trial)),byrow=T), #matrix of trial numbers for each participant
              acc = matrix(experiment$Accuracy, nrow=length(unique(experiment$ID)), ncol=length(unique(experiment$Running_Trial)), byrow=T) #matrix of accuracies for each participant
  )
  return(dat)
}

################################################
############ LOAD AND FORMAT DATA ##############
################################################

# This section reads in and formats data so that it can be fit using rstan.

# Define file paths and object names
file_names <- paste0("data/explicit/explicit_", 1:5, ".rds")
object_names <- paste0("exp", 1:5)

# Read all datasets into a named list
exp_list <- file_names %>%
  set_names(object_names) %>%
  map(readRDS)

# When fitting stan models, we'll need to be able to index participants
# and conditions. So note that the variables ID and Condition
# have corresponding ID_num and Condition_num variables in which they
# are formatted as consecutive integers starting at 1
head(exp_list[[1]])


# For fitting the stan models, data need to be formatted in a list where each 
# element is a variable that is input to the model. Here, we use the function
# defined above, gen_explicit_data(), to turn data for each experiment into a list.

# Define number of comparisons between conditions for each experiment
num_comparisons_lookup <- c(exp1 = 2, exp2 = 1, exp3 = 1, exp4 = 1, exp5 = 1)

# Apply gen_explicit_data() to each experiment in the list
dat_exp_list <- imap(exp_list, function(df, name) {
    gen_explicit_data(df, num_comparisons_lookup[[name]])})


################################################
################ FIT MODELS ####################
################################################

# The models take a long time to fit, so it is recommended
# to save the fits and load them in future. Comment-out
# lines that fit models after fitting as needed.

# Fits are done one at a time. It is not recommended
# that you combine this in a tidy pipeline because fitting
# multiple stan models in quick succession can sometimes crash R.

niter <- 8000
burnin <- 1000
nchains <- 1

######################### COMPILE MODELS ############################

# The models are saved as stan files. Here we load and compile them
# before fitting.

path_models <- "stan_models/explicit/explicit_model_"

int_slope_acc_compiled <- stan_model(file = paste(path_models, "intercept_slope.stan", sep = ""))
int_acc_compiled <- stan_model(file = paste(path_models, "intercept.stan", sep = ""))
slope_acc_compiled <- stan_model(file = paste(path_models, "slope.stan", sep = ""))
base_acc_compiled <- stan_model(file = paste(path_models, "base.stan", sep = ""))


############################# EXP 1 ################################
# Intercept and slope
fit_int_slope_exp1 <- sampling(object = int_slope_acc_compiled,
                           data=dat_exp_list$exp1,iter=niter+burnin,
                           warmup=burnin,chains=nchains,
                           control = list(adapt_delta = 0.99, max_treedepth=15))
saveRDS(fit_int_slope_exp1, file="stan_fits/explicit/fit_int_slope_exp1.rds")

# fit_int_slope_exp1 <- readRDS(file="stan_fits/explicit/fit_int_slope_exp1.rds")

# Intercept only
fit_int_exp1 <- sampling(object= int_acc_compiled,
                         data=dat_exp_list$exp1,iter=niter+burnin,
                         warmup=burnin,chains=nchains,
                         control = list(adapt_delta = 0.99, max_treedepth=15))
saveRDS(fit_int_exp1, file="stan_fits/explicit/fit_int_exp1.rds")

# fit_int_exp1 <- readRDS(file="stan_fits/explicit/fit_int_exp1.rds")

# Slope only
fit_slope_exp1 <- sampling(object = slope_acc_compiled,
                           data=dat_exp_list$exp1,iter=niter+burnin,
                           warmup=burnin,chains=nchains,
                           control = list(adapt_delta = 0.99, max_treedepth = 15))
saveRDS(fit_slope_exp1, file="stan_fits/explicit/fit_slope_exp1.rds")

# fit_slope_exp1 <- readRDS(file="stan_fits/explicit/fit_slope_exp1.rds")

# Baseline
fit_base_exp1 <- sampling(object = base_acc_compiled,
                          data=dat_exp_list$exp1,iter=niter+burnin,
                          warmup=burnin,chains=nchains,
                          control = list(adapt_delta = 0.99, max_treedepth = 15))
saveRDS(fit_base_exp1, file="stan_fits/explicit/fit_base_exp1.rds")

# fit_base_exp1 <- readRDS(file="stan_fits/explicit/fit_base_exp1.rds")

############################# EXP 2 ################################

# Intercept and slope
fit_int_slope_exp2 <- sampling(object = int_slope_acc_compiled,
                               data=dat_exp_list$exp2,iter=niter+burnin,
                               warmup=burnin,chains=nchains,
                               control = list(adapt_delta = 0.99, max_treedepth=15))
saveRDS(fit_int_slope_exp2, file="stan_fits/explicit/fit_int_slope_exp2.rds")

# fit_int_slope_exp2 <- readRDS(file="stan_fits/explicit/fit_int_slope_exp2.rds")

# Intercept only
fit_int_exp2 <- sampling(object= int_acc_compiled,
                         data=dat_exp_list$exp2,iter=niter+burnin,
                         warmup=burnin,chains=nchains,
                         control = list(adapt_delta = 0.99, max_treedepth=15))
saveRDS(fit_int_exp2, file="stan_fits/explicit/fit_int_exp2.rds")

# fit_int_exp2 <- readRDS(file="stan_fits/explicit/fit_int_exp2.rds")

# Slope only
fit_slope_exp2 <- sampling(object = slope_acc_compiled,
                           data=dat_exp_list$exp2,iter=niter+burnin,
                           warmup=burnin,chains=nchains,
                           control = list(adapt_delta = 0.99, max_treedepth = 15))
saveRDS(fit_slope_exp2, file="stan_fits/explicit/fit_slope_exp2.rds")

# fit_slope_exp2 <- readRDS(file="stan_fits/explicit/fit_slope_exp2.rds")

# Baseline
fit_base_exp2 <- sampling(object = base_acc_compiled,
                          data=dat_exp_list$exp2,iter=niter+burnin,
                          warmup=burnin,chains=nchains,
                          control = list(adapt_delta = 0.99, max_treedepth = 15))
saveRDS(fit_base_exp2, file="stan_fits/explicit/fit_base_exp2.rds")

# fit_base_exp2 <- readRDS(file="stan_fits/explicit/fit_base_exp2.rds")

############################# EXP 3 ################################

# Intercept and slope
fit_int_slope_exp3 <- sampling(object = int_slope_acc_compiled,
                               data=dat_exp_list$exp3,iter=niter+burnin,
                               warmup=burnin,chains=nchains,
                               control = list(adapt_delta = 0.99, max_treedepth=15))
saveRDS(fit_int_slope_exp3, file="stan_fits/explicit/fit_int_slope_exp3.rds")

# fit_int_slope_exp3 <- readRDS(file="stan_fits/explicit/fit_int_slope_exp3.rds")

# Intercept only
fit_int_exp3 <- sampling(object= int_acc_compiled,
                         data=dat_exp_list$exp3,iter=niter+burnin,
                         warmup=burnin,chains=nchains,
                         control = list(adapt_delta = 0.99, max_treedepth=15))
saveRDS(fit_int_exp3, file="stan_fits/explicit/fit_int_exp3.rds")

# fit_int_exp3 <- readRDS(file="stan_fits/explicit/fit_int_exp3.rds")

# Slope only
fit_slope_exp3 <- sampling(object = slope_acc_compiled,
                           data=dat_exp_list$exp3,iter=niter+burnin,
                           warmup=burnin,chains=nchains,
                           control = list(adapt_delta = 0.99, max_treedepth = 15))
saveRDS(fit_slope_exp3, file="stan_fits/explicit/fit_slope_exp3.rds")

# fit_slope_exp3 <- readRDS(file="stan_fits/explicit/fit_slope_exp3.rds")

# Baseline
fit_base_exp3 <- sampling(object = base_acc_compiled,
                          data=dat_exp_list$exp3,iter=niter+burnin,
                          warmup=burnin,chains=nchains,
                          control = list(adapt_delta = 0.99, max_treedepth = 15))
saveRDS(fit_base_exp3, file="stan_fits/explicit/fit_base_exp3.rds")

# fit_base_exp3 <- readRDS(file="stan_fits/explicit/fit_base_exp3.rds")

############################# EXP 4 ################################

# Intercept and slope
fit_int_slope_exp4 <- sampling(object = int_slope_acc_compiled,
                               data=dat_exp_list$exp4,iter=niter+burnin,
                               warmup=burnin,chains=nchains,
                               control = list(adapt_delta = 0.99, max_treedepth=15))
saveRDS(fit_int_slope_exp4, file="stan_fits/explicit/fit_int_slope_exp4.rds")

# fit_int_slope_exp4 <- readRDS(file="stan_fits/explicit/fit_int_slope_exp4.rds")

# Intercept only
fit_int_exp4 <- sampling(object= int_acc_compiled,
                         data=dat_exp_list$exp4,iter=niter+burnin,
                         warmup=burnin,chains=nchains,
                         control = list(adapt_delta = 0.99, max_treedepth=15))
saveRDS(fit_int_exp4, file="stan_fits/explicit/fit_int_exp4.rds")

# fit_int_exp4 <- readRDS(file="stan_fits/explicit/fit_int_exp4.rds")

# Slope only
fit_slope_exp4 <- sampling(object = slope_acc_compiled,
                           data=dat_exp_list$exp4,iter=niter+burnin,
                           warmup=burnin,chains=nchains,
                           control = list(adapt_delta = 0.99, max_treedepth = 15))
saveRDS(fit_slope_exp4, file="stan_fits/explicit/fit_slope_exp4.rds")

# fit_slope_exp4 <- readRDS(file="stan_fits/explicit/fit_slope_exp4.rds")

# Baseline
fit_base_exp4 <- sampling(object = base_acc_compiled,
                          data=dat_exp_list$exp4,iter=niter+burnin,
                          warmup=burnin,chains=nchains,
                          control = list(adapt_delta = 0.99, max_treedepth = 15))
saveRDS(fit_base_exp4, file="stan_fits/explicit/fit_base_exp4.rds")

# fit_base_exp4 <- readRDS(file="stan_fits/explicit/fit_base_exp4.rds")

############################# EXP 5 ################################

# Intercept and slope
fit_int_slope_exp5 <- sampling(object = int_slope_acc_compiled,
                               data=dat_exp_list$exp5,iter=niter+burnin,
                               warmup=burnin,chains=nchains,
                               control = list(adapt_delta = 0.99, max_treedepth=15))
saveRDS(fit_int_slope_exp5, file="stan_fits/explicit/fit_int_slope_exp5.rds")

# fit_int_slope_exp5 <- readRDS(file="stan_fits/explicit/fit_int_slope_exp5.rds")

# Intercept only
fit_int_exp5 <- sampling(object= int_acc_compiled,
                         data=dat_exp_list$exp5,iter=niter+burnin,
                         warmup=burnin,chains=nchains,
                         control = list(adapt_delta = 0.99, max_treedepth=15))
saveRDS(fit_int_exp5, file="stan_fits/explicit/fit_int_exp5.rds")

# fit_int_exp5 <- readRDS(file="stan_fits/explicit/fit_int_exp5.rds")

# Slope only
fit_slope_exp5 <- sampling(object = slope_acc_compiled,
                           data=dat_exp_list$exp5,iter=niter+burnin,
                           warmup=burnin,chains=nchains,
                           control = list(adapt_delta = 0.99, max_treedepth = 15))
saveRDS(fit_slope_exp5, file="stan_fits/explicit/fit_slope_exp5.rds")

# fit_slope_exp5 <- readRDS(file="stan_fits/explicit/fit_slope_exp5.rds")

# Baseline
fit_base_exp5 <- sampling(object = base_acc_compiled,
                          data=dat_exp_list$exp5,iter=niter+burnin,
                          warmup=burnin,chains=nchains,
                          control = list(adapt_delta = 0.99, max_treedepth = 15))
saveRDS(fit_base_exp5, file="stan_fits/explicit/fit_base_exp5.rds")

# fit_base_exp5 <- readRDS(file="stan_fits/explicit/fit_base_exp5.rds")


################################################
############### ASSESS MODELS ##################
################################################

############## COMPARE MODELS ##################

# Define the model types and experiment numbers
model_types <- c("base", "int", "slope", "int_slope")
experiments <- 1:5

# Compute LOO for each model and store in long format
loo_long <- expand.grid(exp = experiments, model = model_types) %>%
  mutate(
    fit_name = paste0("fit_", model, "_exp", exp),
    loo_obj = map(fit_name, ~ loo(get(.x)))
  )

# Reshape to wide format with one row per experiment
loo_wide <- loo_long %>%
  select(exp, model, loo_obj) %>%
  pivot_wider(names_from = model, values_from = loo_obj)

# Perform loo_compare across models for each experiment
loo_comparisons <- loo_wide %>%
  mutate(comparison = pmap(list(int_slope, int, slope, base), loo_compare)) %>%
  select(exp, comparison)

# Combine the comparison results into a dataframe that includes columns
# indicating the name of the model and its ranking in the results of
# loo_compare. The model with rank 1 is the winner.
loo_results <- map2_dfr(
  loo_comparisons$exp,
  loo_comparisons$comparison,
  
  function(exp, comp) {
    model_key <- data.frame(model = c("int_slope", "int", "slope", "base"))
    rownames(model_key) <- c("model1", "model2", "model3", "model4")
    
    comp_df <- as_tibble(comp, rownames = "model_id") %>%
      full_join(as_tibble(model_key, rownames = "model_id"), by = "model_id") %>%
      mutate(exp = exp,
             ranking = row_number())
    
    return(comp_df)
  }
)

# Look at the winning model for each experiment
loo_winners <- loo_results %>%
  group_by(exp) %>%
  dplyr::summarise(winner = model[ranking == 1])

loo_winners

# Write the winner results to a file
write.csv(loo_winners, file = "stan_fits/explicit/explicit_winners.csv", row.names = F)


################################################
############ CONDITION DIFFERENCES #############
################################################

# Combine fits and the parameter containing condition differences
# in each posterior sample into a tibble
posterior_info <- tibble::tibble(
  exp = 1:5,
  model_obj = list(fit_slope_exp1, fit_slope_exp2, fit_int_exp3, fit_slope_exp4, fit_slope_exp5),
  parameter = c("post_slope_diff", "post_slope_diff", "post_alpha_diff", "post_slope_diff", "post_slope_diff"),
  rename_cols = list(
    c("post_diff_dense", "post_diff_sparse"),
    c("post_diff_dense"),
    c("post_diff_dense"),
    c("post_diff_dense"),
    c("post_diff_dense")
  )
)


########### INCIDENTAL VS BASELINE #############

# Extract samples from posterior for differences between incidental 
# and baseline conditions from best fitting model for each experiment 
# (slope for experiments 1, 2, and 4, intercept for experiment 3)

# Extract posterior samples
posterior_samples <- posterior_info %>%
  mutate(
    samples = map2(model_obj, parameter, ~ as.data.frame(as.matrix(.x, pars = .y))),
    samples = map2(samples, rename_cols, ~ setNames(.x, .y))
  )


# Calculate 89% credible intervals for the posterior distributions
# of condition differences
hdi_results <- posterior_samples %>%
  mutate(
    hdi = map(samples, ~ map_dfr(.x, ~ bayestestR::ci(.x, method = "HDI", ci = .89), .id = "term"))
  ) %>%
  select(exp, hdi) %>%
  unnest(hdi)

# Look:
hdi_results

# From each posterior, get the probability that the difference
# was greater than 0
prob_diff <- posterior_samples %>%
  mutate(
    prob_diff = map(samples, ~ map_dbl(.x, ~ mean(.x > 0)))
  ) %>%
  transmute(
    exp,
    term = map(rename_cols, ~ .x),
    prob_diff
  ) %>%
  unnest(c(term, prob_diff))



# Combine HDI and probability results
posterior_summary <- left_join(hdi_results, prob_diff, by = c("exp", "term"))

# Look at results:
posterior_summary


############## DENSE VS SPARSE ################

# Learnability of dense vs sparse categories in baseline condition in Experiment 1
post_slope_exp1 <- data.frame(as.matrix(fit_slope_exp1, pars=c("condition_slope_mean[1]", "condition_slope_mean[3]")))

names(post_slope_exp1) <- c("Base_Dense", "Base_Sparse")

post_slope_exp1$Diff <- post_slope_exp1$Base_Dense - post_slope_exp1$Base_Sparse

bayestestR::ci(post_slope_exp1$Diff, method="HDI")



