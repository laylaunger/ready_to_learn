# Title: Exposure Phase Analysis
# Author: Layla Unger
# R version: 4.3.2
# Packages: readr, stringr, psycho, dplyr, tidyr, purrr, ggplot2, rstan, bayestestR, ggmcmc

################################################
################ DESCRIPTION ###################
################################################

# This script is part of a project that investigated whether mere passive exposure
# to new categories fostered subsequent explicit learning.

# This project included a set of 5 experiments that each involved two phases: an
# Exposure phase involving passive exposure, and an Explicit phase that tested
# subsequent learning. This script analyzes data from the Exposure phase.

# Note that only Experiments 1 - 3 included an "Exposure Phase" in which behavioral
# data were recorded. Thus, this script focuses on Experiments 1-3.

# In the Exposure phase, participants completed a cover task. During the cover task, 
# participants saw a stimulus "jump" to either the left or right. After the stimulus 
# had jumped, participants responded indicating whether it jumped to the left or right.
# Across experiments, participants completed this cover task under one of two exposure 
# conditions:
# -- Baseline, in which the stimuli had no category structure, and 
# -- Incidental, in which, unbeknownst to participants, the stimuli belonged to one of 
#    two categories.

# In Experiment 1, categories were either Dense or Sparse. Data for these two 
# "category structure" conditions are stored as Experiment1A and Experiment1B. Only 
# Dense categories were used in Experiments 2 and 3.

# Goal of analyses: Test whether exposure condition affected performance on the cover task.
# Specifically, test whether incidental exposure to dense or sparse categories caused 
# participants to speed up more rapidly than baseline exposure.

# Type of analyses: Bayesian hierarchical modeling

# Data for each experiment are stored in a long- format dataframe. Important columns are:
# ID
# Condition
# Accuracy (0 or 1)
# Running_Trial (Trial 1-30)
# RT

# Steps currently in script:

# (1) Define a set of new functions that will be used in the script to format data for
# Bayesian hierarchical modeling, extract model fits and plot results

# (2) Load a set of hierarchical Bayesian models that will be used to analyze the data. 
# In these models, the outcome variable (RT or Accuracy) is predicted by trial number. 
# Each participant has a slope for change in the outcome variable across trials. Participant 
# slopes are drawn from one of two distributions, depending on whether they belonged to the 
# Baseline or Incidental condition. For each iteration, the model compares the Baseline and 
# Incidental slope distributions. This comparison is used to compare whether changes in RT or 
# Accuracy across trials vary across the two conditions.

# (3) Load and format RT and Accuracy data from each experiment.

# (4) Fit the models to the data

# (5) From the models, extract differences between the distributions for slopes in the two
# conditions for each sample from the posterior.

# (6) Evaluate differences between slope distributions between the two conditions for each
# experiment.

################################################
################# PACKAGES #####################
################################################

library(readr)
library(stringr)
library(psycho)

library(dplyr)
library(tidyr)
library(purrr)

library(ggplot2)

library(rstan)
library(bayestestR)
library(ggmcmc)

rstan_options(auto_write = TRUE)
options(mc.cores = 3)

################################################
################# FUNCTIONS ####################
################################################

# This section defines functions that will be used in the script.

############## FORMATTING DATA #################

# Stan models need data in a list format. 

# This function generates the list format needed to analyse RT data.
# A different list format is needed for experiments 1A, 1B and 2 versus
# experiment 3. This function includes an argument for specifying the
# format.
gen_rt_data <- function(experiment, standard_rt = TRUE) {
  
  #Identify each participant's condition
  sub_condition <- experiment %>%
    group_by(ID) %>%
    dplyr::summarise(Cond = unique(Condition_num)) %>%
    pull(Cond)
  
  if(standard_rt) {
    dat_experiment <- list(num_sub = length(unique(experiment$ID_num)), #number of participants
                           num_trials = length(unique(experiment$Running_Trial)), #number of trials for each participant
                           num_conditions = length(unique(experiment$Condition_num)), #number of conditions
                           sub_condition = sub_condition,
                           trial = matrix(experiment$Running_Trial, nrow = length(unique(experiment$ID)), ncol = length(unique(experiment$Running_Trial)) , byrow = T), #matrix of trial numbers for each participant
                           trial_acc  =  matrix(experiment$Acc_num, nrow = length(unique(experiment$ID)), ncol = length(unique(experiment$Running_Trial)) , byrow = T), #matrix of trial accuracies for each participant
                           rt = matrix(experiment$RT, nrow = length(unique(experiment$ID)), ncol = length(unique(experiment$Running_Trial)), byrow = T) #matrix of rts for each participant
    )
  } else {
    data_experiment = list(num_sub = length(unique(exp3$ID_num)), #number of participants
                           num_obs = nrow(exp3[!(is.na(exp3$RT_filter)),]), #total number of non-filtered RT observations across particpants
                           num_trials = length(unique(exp3$Running_Trial)), # number of trials in exposure phase
                           num_conditions = length(unique(exp3$Condition_num)), #number of conditions
                           sub_condition = sub_condition, #condition for each participant
                           id = exp3$ID_num[!is.na(exp3$RT_filter)], #vector of ids across observations for non-filtered RTs
                           trial = exp3$Running_Trial[!is.na(exp3$RT_filter)], #vector of trial numbers across observations for non-filtered RTs
                           rt = exp3$RT_filter[!is.na(exp3$RT_filter)] #vector of RTs across observations for non-filtered RTs
    )
  }
}

# This function generates the list format needed to analyse RT data.
gen_acc_data <- function(experiment) {
  
  #Identify each participant's condition
  sub_condition <- experiment %>%
    group_by(ID) %>%
    dplyr::summarise(Cond = unique(Condition_num)) %>%
    pull(Cond)
  
  dat_experiment <- list(num_sub = length(unique(experiment$ID_num)), #number of participants
                         num_trials = length(unique(experiment$Running_Trial)), #number of trials for each participant
                         num_conditions = length(unique(experiment$Condition_num)), #number of conditions
                         sub_condition = sub_condition, #condition for each participant
                         trial = matrix(experiment$Running_Trial, nrow = length(unique(experiment$ID)), ncol = length(unique(experiment$Running_Trial)),byrow = T), #matrix of trial numbers for each participant
                         acc = matrix(experiment$Accuracy, nrow = length(unique(experiment$ID)), ncol = length(unique(experiment$Running_Trial)), byrow = T) #matrix of accuracies for each participant
  )
}

######### SIGNAL DETECTION MEASURES ############

# For experiments that used a 1-back task, for each participant,
# this function calculates numbers of hits (detecting repetitions),
# correct rejections (correct rejections of non-repetitions), 
# false alarms (identifying non-repetitions as repetitions), and
# misses (failing to detect repetitions)
# Use to calculate hit rates and false alarm rates, hit - false alarms,
# and dprime

extract_signal <- function(input) {
  n_hit <- length(input$Accuracy[input$Rep == 1 & input$Accuracy == 1])
  n_cr <- length(input$Accuracy[input$Rep == 0 & input$Accuracy == 1])
  n_fa <- length(input$Accuracy[input$Rep == 0 & input$Accuracy == 0])
  n_miss <- length(input$Accuracy[input$Rep ==1 & input$Accuracy == 0])
  hit <- n_hit / (n_hit + n_miss)
  fa <- n_fa / (n_fa + n_cr)
  indeces <- psycho::dprime(n_hit, n_fa, n_miss, n_cr)
  d.prime <- indeces$dprime
  hit_fa <- hit - fa
  output <- data.frame(hit, fa, d.prime, hit_fa)
  return(output)
}



################################################
############ LOAD AND FORMAT DATA ##############
################################################

# This section reads in and formats data so that it can be fit using rstan. 
# Note that the data for Experiment 1 is divided into two objects: exp1A, 
# which is the Dense Structure condition, and exp1B, which is the Sparse 
# Structure condition.

# Define file paths and object names
exp_names <- c("1A", "1B", "2", "3")
file_names <- paste0("data/exposure/exposure_", exp_names, ".rds")
object_names <- c("exp1A", "exp1B", "exp2", "exp3")

# Read all datasets into a named list
exp_list <- file_names %>%
  set_names(object_names) %>%
  map(readRDS)

# Note that the data are in the same format for experiments 1A, 1B, and 2, 
# but slightly different for experiment 3. 
# This is because experiments 1A, 1B, and 2 had progressively strict limits 
# on reaction time, and therefore: 
# (1) Have a sizable minority of inaccurate trials (~15%-20%), and 
# (2) No outlier reaction times.
# In contrast, experiment 3 had no reaction time limit, and therefore has the 
# opposite pattern: 
# (1) Almost no inaccurate trials (0.6%), and 
# (2) Some outlier reaction times.
# 
# Therefore, RT data used to fit models from experiments 1A, 1B, and 2 include all 
# RTs, such that inaccurate trials are simply fit to a different distribution from 
# accurate trials. In addition, accuracy data from these experiments are also used 
# to fit a model for accuracy.
# In contrast, RT data used to fit models from experiment 3 include only accurate 
# trial RTs (with very long RTs removed), and inaccurate trials are not modeled. 

exp_list$exp3 <- exp_list$exp3 %>%
  group_by(ID) %>%
  dplyr::mutate(RT_filter = ifelse(RT > 2000 | Accuracy == 0, NA, RT))


# When fitting stan models, we'll need to be able to index participants,
# conditions and accuracies. So note that these variables have
# corresponding ID_num and Condition_num variables in which they
# are formatted as consecutive integers starting at 1
# There are always 2 conditions, corresponding to the Baseline
# and Category conditions
head(exp_list[[1]])


# Data used to fit the models

# Note: the progressively strict RT limits in experiments 1A, 1B, and 2 meant 
# that changes in RT over trials in these experiments were approximately linear. 
# Therefore, these changes were modeled using linear regression.
# In contrast, without RT limits in experiment 3, changes in RT showed a typical 
# pattern of speeding up rapidly at first, and then asymptoting. Therefore, these 
# changes were modeled by predicting RT as the outcome of an asymptotic function with 
# intercept, rate, and asymptote parameters.

# Reaction time data

# Create a named logical vector for that identifies, for each experiment, whether
# to use standard RT (all except experiment 3)
rt_flags <- c(exp1A = TRUE, exp1B = TRUE, exp2 = TRUE, exp3 = FALSE)

# Apply the function using map2
dat_rt_list <- map2(exp_list, rt_flags, gen_rt_data)

# Accuracy data

dat_acc_list <- map(exp_list[names(exp_list) != "exp3"], ~ gen_acc_data(.x))


######################################################
#################### FIT MODELS ######################
######################################################

niter <- 100
burnin <- 100
nchains <- 1

# The models take a long time to fit, so it is recommended
# to save the fits and load them in future. Comment-out
# lines that fit models after fitting as needed.

# Fits are done one at a time. It is not recommended
# that you combine this in a tidy pipeline because fitting
# multiple stan models in quick succession can sometimes crash R.

################### COMPILE MODELS ###################

# The models are saved as stan files. Here we load and compile 
# them before fitting.

path_models <- "stan_models/exposure/exposure_model_"

model_slope_acc_compiled <- stan_model(file = paste(path_models, "slope_acc.stan", sep = ""))
model_slope_rt_compiled <- stan_model(file = paste(path_models, "slope_rt.stan", sep = ""))
model_exponent_rt_compiled <- stan_model(file = paste(path_models, "exponent_rt.stan", sep = ""))

################### EXPERIMENT 1A ####################


fit_slope_rt_exp1A <- sampling(object = model_slope_rt_compiled,
                               data = dat_rt_list$exp1A, iter = niter + burnin,
                               warmup = burnin, chains = nchains,
                               verbose = TRUE, control = list(adapt_delta = 0.99, max_treedepth=15))
saveRDS(fit_slope_rt_exp1A, file = "stan_fits/exposure/fit_slope_rt_exp1A.rds")

# fit_slope_rt_exp1A <- readRDS(file = "stan_fits/exposure/fit_slope_rt_exp1A.rds")
traceplot(fit_slope_rt_exp1A, pars="condition_slope_mean")


fit_slope_acc_exp1A <- sampling(object = model_slope_acc_compiled,
                               data = dat_acc_list$exp1A, iter = niter + burnin,
                               warmup = burnin, chains = nchains,
                               verbose = TRUE, control = list(adapt_delta = 0.99, max_treedepth=15))
saveRDS(fit_slope_acc_exp1A, file = "stan_fits/exposure/fit_slope_acc_exp1A.rds")

# fit_slope_acc_exp1A <- readRDS(file = "stan_fits/exposure/fit_slope_acc_exp1A.rds")
traceplot(fit_slope_acc_exp1A, pars="condition_slope_mean")

################### EXPERIMENT 1B ####################

fit_slope_rt_exp1B <- sampling(object = model_slope_rt_compiled,
                               data = dat_rt_list$exp1B, iter = niter + burnin,
                               warmup = burnin, chains = nchains,
                               verbose = TRUE, control = list(adapt_delta = 0.99, max_treedepth=15))
saveRDS(fit_slope_rt_exp1B, file = "stan_fits/exposure/fit_slope_rt_exp1B.rds")

# fit_slope_rt_exp1B <- readRDS(file = "stan_fits/exposure/fit_slope_rt_exp1B.rds")
traceplot(fit_slope_rt_exp1B, pars="condition_slope_mean")


fit_slope_acc_exp1B <- sampling(object = model_slope_acc_compiled,
                                data = dat_acc_list$exp1B, iter = niter + burnin,
                                warmup = burnin, chains = nchains,
                                verbose = TRUE, control = list(adapt_delta = 0.99, max_treedepth=15))
saveRDS(fit_slope_acc_exp1B, file = "stan_fits/exposure/fit_slope_acc_exp1B.rds")

# fit_slope_acc_exp1B <- readRDS(file = "stan_fits/exposure/fit_slope_acc_exp1B.rds")
traceplot(fit_slope_acc_exp1B, pars="condition_slope_mean")

################### EXPERIMENT 2 ####################

fit_slope_rt_exp2 <- sampling(object = model_slope_rt_compiled,
                              data = dat_rt_list$exp2, iter = niter + burnin,
                              warmup = burnin, chains = nchains,
                              verbose = TRUE, control = list(adapt_delta = 0.99, max_treedepth=15))
saveRDS(fit_slope_rt_exp2, file = "stan_fits/exposure/fit_slope_rt_exp2.rds")

# fit_slope_rt_exp2 <- readRDS(file = "stan_fits/exposure/fit_slope_rt_exp2.rds")
traceplot(fit_slope_rt_exp2, pars="condition_slope_mean")


fit_slope_acc_exp2 <- sampling(object = model_slope_acc_compiled,
                               data = dat_acc_list$exp2, iter = niter + burnin,
                               warmup = burnin, chains = nchains,
                               verbose = TRUE, control = list(adapt_delta = 0.99, max_treedepth=15))
saveRDS(fit_slope_acc_exp2, file = "stan_fits/exposure/fit_slope_acc_exp2.rds")

# fit_slope_acc_exp2 <- readRDS(file = "stan_fits/exposure/fit_slope_acc_exp2.rds")
traceplot(fit_slope_acc_exp2, pars="condition_slope_mean")

############################# EXP 3 ################################

fit_exponent_rt_exp3 <- sampling(object = model_exponent_rt_compiled,
                              data = dat_rt_list$exp3, iter = niter + burnin,
                              warmup = burnin, chains = nchains,
                              control = list(adapt_delta = 0.99, max_treedepth=15))
saveRDS(fit_exponent_rt_exp3, file = "stan_fits/exposure/fit_exponent_rt_exp3.rds")

# fit_exponent_acc_exp3 <- readRDS(file = "stan_fits/exposure/fit_exponent_rt_exp3.rds")
traceplot(fit_exponent_rt_exp3, pars=c("condition_rate_mean", "condition_asymptote_mean"))


################################################
############# SIGNAL DETECTION #################
################################################

# Experiments 1-3 used tasks that afforded a test of implicit learning 
# based on RT over trials. In contrast, Experiments 4 and 5 used 1-back 
# tasks that did not afford such tests
# This section contains code to check whether participants in Experiments 
# 4 and 5 were on-task (i.e., performed well at detecting 1-back repetitions, 
# assessed using signal detection indeces)

exp4 <- readRDS(file="data/exposure/exposure_4.rds")
exp5 <- readRDS(file="data/exposure/exposure_5.rds")

exp4_signal <- exp4 %>%
  group_by(ID, Condition) %>%
  group_modify(~ extract_signal(.))

exp5_signal <- exp5 %>%
  group_by(ID, Condition) %>%
  group_modify(~ extract_signal(.))

# Identify number of participants with inattentive performance (fewer hits than 
# false alarms). These participants are removed from analyses.
nrow(exp4_signal[exp4_signal$hit_fa <= 0,])
nrow(exp5_signal[exp5_signal$hit_fa <= 0,])

exp4_signal <- exp4_signal[exp4_signal$hit_fa > 0,]
exp5_signal <- exp5_signal[exp5_signal$hit_fa > 0,]

# Check that remaining participants were on-task (above chance
# hit - fa)
t.test(exp4_signal$hit_fa, mu = 0)
t.test(exp5_signal$hit_fa, mu = 0)

# Descriptive stats
mean(exp4_signal$hit_fa)
sd(exp4_signal$hit_fa)

mean(exp5_signal$hit_fa)
sd(exp5_signal$hit_fa)
