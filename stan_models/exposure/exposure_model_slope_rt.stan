// Name: model_slope_rt
// Model for change in RT across trials in experiments with progressively decreasing limits on RT, 
// in which changes in RT across trials were approximately linear. 
// Due to linear changes in RT, this model predicts RT for each trial for each participant as the 
// outcome of a linear regression with an intercept and slope for trial number. To test condition
// differences, slopes for change in RT across trials are sampled from different distributions based
// on between-subjects condition.

data {
  int<lower=1> num_sub; // number of subjects
  int<lower=1> num_trials; // number of trials
  int<lower=1> num_conditions; // number of conditons
  int<lower=1> sub_condition[num_sub]; // condition for each subject to be used as *indeces* to identify params for different conditions
  matrix<lower=0>[num_sub,num_trials] trial; // trial number predictor - matrix of trial numbers for each subject
  int<lower=0> trial_acc[num_sub,num_trials]; // accuracies for each trial to be used as *indeces* to identify separate params for acc & inacc trials. So acc = 1, inacc = 2
  real<lower=0> rt[num_sub,num_trials]; // reaction times
  }
parameters {
  real sub_alpha[num_sub];  // Intercept for each participant
  vector[num_sub] sub_condition_slope; // Participants' trial number slopes for accurate trials only. To be drawn from one of two distributions, based on sub_condition
  vector[num_sub] sub_inacc_slope; // Participants' trial number slopes for inaccurate trials
  
  real hyper_alpha_mean;  // mean of the distribution for subject intercepts
  real<lower=0> hyper_alpha_sigma2; // variance of the distribution for subject intercepts
  
  real condition_slope_mean[num_conditions]; // means of the 2 distributions for particpant trial number slopes***
  real<lower=0> condition_slope_sigma2[num_conditions]; // variances of the 2 distributions for particpant trial number slopes
  
  real inacc_slope_mean; // Mean of distribution for participant trial number slopes on *inaccurate* trials
  real<lower=0> inacc_slope_sigma2; // Variance of distribution for participant trial number slopes on *inaccurate* trials
  
  real<lower=0> rt_sigma2; // Error term for predicting trial-level RTs
  }

transformed parameters {
  real<lower=0> hyper_alpha_sigma;
  real<lower=0> condition_slope_sigma[num_conditions];
  real<lower=0> inacc_slope_sigma;
  
  real<lower=0> rt_sigma;
  
  matrix[num_sub, 2] sub_slope; // Trial number slopes for each participant. [,1] = slope for accurate trials based on condition; [,2] = slope for inaccurate trials
 
  hyper_alpha_sigma = sqrt(hyper_alpha_sigma2);
  condition_slope_sigma = sqrt(condition_slope_sigma2);
  inacc_slope_sigma = sqrt(inacc_slope_sigma2);
  
  rt_sigma = sqrt(rt_sigma2);
  
  sub_slope = append_col(sub_condition_slope, sub_inacc_slope); // For each subject, combine their slope for accurate trials based on condition, and their slope for inaccurate trials
  }

model {
// Declare predicted rt variable
  real rt_mu[num_sub, num_trials];

// Priors for hyperparameters
  hyper_alpha_mean ~ normal(0,1/.01);
  condition_slope_mean ~ normal(0,1/.01);
  inacc_slope_mean ~ normal(0,1/.01);
  
  hyper_alpha_sigma2 ~ inv_gamma(.01,.01);
  condition_slope_sigma2 ~ inv_gamma(.01,.01);
  inacc_slope_sigma2 ~ inv_gamma(.01,.01);
  rt_sigma2 ~ inv_gamma(.01,.01);
  
// Priors for subject intercept and slope means and variances
  for (i in 1:num_sub){
    sub_alpha[i] ~ normal(hyper_alpha_mean, hyper_alpha_sigma);
    sub_condition_slope[i] ~ normal(condition_slope_mean[sub_condition[i]], condition_slope_sigma[sub_condition[i]]);
    sub_inacc_slope[i] ~ normal(inacc_slope_mean, inacc_slope_sigma);
  }

// Likelihood for reaction time: Predicted by participant's intercept, their slope for accurate trials (based on condition), and their slope for inacc trials
  for (i in 1:num_sub){                           
    for (j in 1:num_trials){
      rt_mu[i,j] = sub_alpha[i] + sub_slope[i,trial_acc[i,j]] * trial[i,j];
      rt[i,j] ~ normal(rt_mu[i,j], rt_sigma);
    }
  }
}

generated quantities {
  real post_slope_diff; // difference between means for slopes for baseline vs dense or sparse conditions
  matrix[num_sub,num_trials] pred_rt; // matrix of predicted RTs for each subject/trial
  
  post_slope_diff = condition_slope_mean[1] - condition_slope_mean[2];
    
  for (i in 1:num_sub) {
    for (j in 1:num_trials) {
      pred_rt[i,j] = normal_rng(sub_alpha[i] + sub_slope[i,trial_acc[i,j]] * trial[i,j], rt_sigma); // predicted RT on jth trial for the ith subject
    }
  }
}
