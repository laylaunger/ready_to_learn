// Name: model_exponent_rt
// Model for change in RT across trials in experiments with no limits on RT,
// in which changes in RT across trials were approximately exponential. 
// Due to nonlinear changes in RT, this model predicts RT for each trial based on
// a function with parameters: asymptote, start, and rate. 
// To test condition differences, the parameters are sampled from different 
// distributions based on between-subjects condition.

data {
  // Observed data
  int<lower=1> num_sub; // number of subjects
  int<lower=1> num_obs; // total number of observations across subjects
  int<lower=1> num_trials; // number of trials in exposure phase
  int<lower=1> num_conditions; // number of conditons
  int<lower=1> sub_condition[num_sub]; // condition for each subject to be used as *indeces* to identify params for different conditions
  int<lower=1> id[num_obs]; // subject id for each observation
  int<lower=1> trial[num_obs]; // trial number predictor across all observations
  real<lower=0> rt[num_obs]; // reaction times across all observations
}

parameters{
  // Participant-level asymptote, gain, rate, and noise  
  real<upper=2000>              sub_asymptote[num_sub];
  real<lower=200, upper=20000>  sub_start[num_sub];
  real<lower=0, upper=50>       sub_rate[num_sub];
    
  real<lower=0>   rt_sigma2; // RT noise
  
  // Group-level asymptote, gain, and rate  
  real<upper=2000>             condition_asymptote_mean[num_conditions];
  real<lower=200, upper=2000>  condition_start_mean[num_conditions];
  real<lower=0, upper=50>      condition_rate_mean[num_conditions];
  
  real<lower=0>   condition_asymptote_sigma2[num_conditions];
  real<lower=0>   condition_start_sigma2[num_conditions];
  real<lower=0>   condition_rate_sigma2[num_conditions];

}

transformed parameters {
  real<lower=0> rt_sigma;

  real<lower=0>   condition_asymptote_sigma[num_conditions];
  real<lower=0>   condition_start_sigma[num_conditions];
  real<lower=0>   condition_rate_sigma[num_conditions];
   
  rt_sigma = sqrt(rt_sigma2);
   
  condition_asymptote_sigma = sqrt(condition_asymptote_sigma2);
  condition_start_sigma = sqrt(condition_start_sigma2);
  condition_rate_sigma = sqrt(condition_rate_sigma2);
  
}

model {
  // Latent RT
  real rt_mu[num_obs];
  
  rt_sigma2 ~ inv_gamma(.01,.01);
  
  // Group-level parameters
  for(c in 1:num_conditions) {
    condition_asymptote_mean[c] ~ normal( 400 , 1/.001 ) T[0 , 2000]; 
    condition_start_mean[c] ~     normal( 1000 , 1/.001 ) T[200, 2000];
    condition_rate_mean[c] ~      normal( 1, 1/.03 ) T[0, 50];
  
    condition_asymptote_sigma2[c] ~ inv_gamma(.01,.01);
    condition_start_sigma2[c] ~ inv_gamma(.01,.01);
    condition_rate_sigma2[c] ~ inv_gamma(.01,.01);
  }
  
  // Participant-level parameters
  
  for(i in 1:num_sub) {
    sub_asymptote[i] ~ normal( condition_asymptote_mean[sub_condition[i]], condition_asymptote_sigma[sub_condition[i]]);
    sub_start[i] ~ normal( condition_start_mean[sub_condition[i]], condition_start_sigma[sub_condition[i]]);
    sub_rate[i] ~ normal( condition_rate_mean[sub_condition[i]], condition_rate_sigma[sub_condition[i]]);
    
  }
  
  //Likelihood: Loop over observations
  for (i in 1:num_obs) {

    rt_mu[i] = sub_asymptote[id[i]] - (sub_asymptote[id[i]] - sub_start[id[i]]) * exp( (1 - trial[i])/sub_rate[id[i]] );

    rt[i] ~   normal( rt_mu[i] , rt_sigma );
    }
}

generated quantities {
  real condition_gain_mean[num_conditions];
  real sub_gain[num_sub];
  real pred_mu[num_obs];
  real pred_rt[num_obs]; 
  real log_lik[num_obs];
  real sub_slope_trial[num_obs];
  real slope_trial_base[num_trials];
  real slope_trial_dense[num_trials];
  real slope_mean_base;
  real slope_mean_dense;
  real slope_head_base;
  real slope_head_dense;
  real slope_tail_base;
  real slope_tail_dense;
  real slope_diff;
  real slope_head_diff;
  real slope_tail_diff;
  real asymptote_diff; 
  real gain_diff;
  real rate_diff;
  
  
  // Calculate gains
  for (c in 1:num_conditions) {
    condition_gain_mean[c] = condition_asymptote_mean[c] - condition_start_mean[c];
  }
  
  for (i in 1:num_sub) {
    sub_gain[i] = sub_asymptote[i] - sub_start[i];
  }
  
  // Calculate trial slopes
  for(t in 1:num_trials) {
    slope_trial_base[t] = (condition_gain_mean[1] * exp( (1 - trial[t])/condition_rate_mean[1] )) / condition_rate_mean[1];
    slope_trial_dense[t] = (condition_gain_mean[2] * exp( (1 - trial[t])/condition_rate_mean[2] )) / condition_rate_mean[2];
  }
  
  
  slope_mean_base = mean(slope_trial_base);
  slope_mean_dense = mean(slope_trial_dense);
  slope_head_base = mean( head(slope_trial_base, 40) );
  slope_head_dense = mean( head(slope_trial_dense, 40) );
  slope_tail_base = mean( tail(slope_trial_base, 40) );
  slope_tail_dense = mean( tail(slope_trial_dense, 40) );
  
  slope_diff = slope_mean_base - slope_mean_dense;
  slope_head_diff = slope_head_base - slope_head_dense;
  slope_tail_diff = slope_tail_base - slope_tail_dense;
  
  asymptote_diff = condition_asymptote_mean[1] - condition_asymptote_mean[2];
  gain_diff =      condition_gain_mean[1] - condition_gain_mean[2];
  rate_diff =      condition_rate_mean[1] - condition_rate_mean[2];
  
  for (i in 1:num_obs){ 
    pred_mu[i] = sub_asymptote[id[i]] - ( sub_asymptote[id[i]] - sub_start[id[i]] ) * exp( (1 - trial[i])/sub_rate[id[i]] );
    pred_rt[i] = normal_rng(pred_mu[i], rt_sigma );
    log_lik[i] = normal_lpdf(rt[i] | pred_mu[i], rt_sigma );
    sub_slope_trial[i] = (sub_gain[id[i]] * exp( (1 - trial[i])/sub_rate[id[i]] )) / sub_rate[id[i]];
  }
}
