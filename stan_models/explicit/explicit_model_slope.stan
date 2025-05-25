// Model name: model_slope
// Model in which slopes but not intercepts are drawn from different distributions based on condition

data {
  int<lower=1> num_sub; // number of subjects
  int<lower=1> num_trials; // number of trials
  int<lower=1> num_conditions; // number of conditons
  int<lower=1> num_comparisons; // number of comparisons between conditions; 2 in exp1 (baseline vs dense and base vs sparse), 1 in other experiments
  int<lower=1> sub_condition[num_sub]; // condition for each subject
  matrix<lower=0>[num_sub,num_trials] trial; // trial number predictor
  int<lower=0,upper=1> acc[num_sub,num_trials]; // accuracies (binary integers - 0 or 1)
  }
parameters {
  real alpha[num_sub];        // Intercept for each participant
  real sub_slope[num_sub]; // Slope for trial number for each particpant
  
  real hyper_alpha;  // mean of the single distribution for subject intercepts
  real<lower=0> hyper_alpha_sigma2; // variance of the single distribution for subject intercepts
  real condition_slope_mean[num_conditions];   // means of the 2 distributions for for subject trial slopes
  real<lower=0> condition_slope_sigma2[num_conditions];  // variances of the 2 distributions for subject trial slopes
  }

transformed parameters {
  real<lower=0> hyper_alpha_sigma;
  real<lower=0> condition_slope_sigma[num_conditions];
  
  hyper_alpha_sigma = sqrt(hyper_alpha_sigma2);
  condition_slope_sigma = sqrt(condition_slope_sigma2);
  }

model {
// Priors for hyperparameters
  hyper_alpha ~ normal(0,1/.01);
  condition_slope_mean ~ normal(0,1/.01);
  
  hyper_alpha_sigma2 ~ inv_gamma(.01,.01);
  condition_slope_sigma2 ~ inv_gamma(.01,.01);


// Priors for intercepts and slopes for each participant
// Use this to draw intercept from one distribution, and slope from one of two distributions based on participant condition
  for (i in 1:num_sub){                           
    alpha[i] ~ normal(hyper_alpha,hyper_alpha_sigma);
    sub_slope[i] ~ normal(condition_slope_mean[sub_condition[i]],condition_slope_sigma[sub_condition[i]]);

// Likelihood for accuracy: Predicted by participant's intercept, their slope for trial type, and their slope for trial num
    for (j in 1:num_trials){                         
    acc[i,j] ~ bernoulli_logit(alpha[i] + sub_slope[i] * trial[i,j]);
    }
  }
}

generated quantities {
  real post_slope_diff[num_comparisons]; // difference between means of conditon distributions for slopes  
  matrix[num_sub,num_trials] log_lik;  // matrix of log-likelihoods for each acc observation
  real ppd_mean_acc[num_sub]; // vector of predicted mean accuracies
  matrix[num_sub,num_trials] pred_prob_acc; // matrix of predicted probabilities of acc response for each subject/trial
  
  //Calculate difference between means of condition distributions for slopes
  if(num_comparisons == 1){
    post_slope_diff[1] = condition_slope_mean[2] - condition_slope_mean[1];
  }
  else {
    post_slope_diff[1] = condition_slope_mean[2] - condition_slope_mean[1];
    post_slope_diff[2] = condition_slope_mean[4] - condition_slope_mean[3];
  }
  
  //Calculate log likelihoods, predicted mean accuracies, and predicted probabilities of a correct response on each trial 
  for (i in 1:num_sub) {
    vector[num_trials] ppd_accuracies; //vector of predicted accuracies for each trial
    
    for (j in 1:num_trials) {
      log_lik[i,j] =  bernoulli_logit_lpmf(acc[i,j] | alpha[i] + sub_slope[i] * trial[i,j] );
      pred_prob_acc[i,j] = inv_logit(alpha[i] + sub_slope[i] * trial[i,j]); // probability of acc response on the jth trial for the ith subject
      ppd_accuracies[j] = bernoulli_logit_rng(alpha[i] + sub_slope[i] * trial[i,j] ); //predicted acc on the jth trial for the ith subject
    }
    ppd_mean_acc[i] = mean(ppd_accuracies); //predicted mean acc for the ith subject
  }
}
