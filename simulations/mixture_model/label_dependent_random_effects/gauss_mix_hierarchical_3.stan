data {
  int<lower=2> ncells;
  int<lower=2> N;
  real y[N];               // observations
  int cell[N];
  simplex[2] theta[ncells];
}
parameters {
  ordered[2] mu;
  ordered[2] mu_obs[ncells];
  vector<lower=0>[2] sigma;  // scales of mixture components
  //vector<lower=0>[2] re_sigma;
  vector<lower=0.01>[2] re_sigma;
}
model {
  vector[2] log_theta[ncells] = log(theta);  // cache log calculation
  for (k in 1:ncells){
    mu_obs[k,1] ~ normal(mu[1], re_sigma[1]);
    mu_obs[k,2] ~ normal(mu[2], re_sigma[2]);
  }
  
  for (n in 1:N) {
    int cur_cell = cell[n];
    vector[2] lps = log_theta[cur_cell];
    vector[2] cur_mu = mu_obs[cur_cell];
    for (k in 1:2)
      lps[k] += normal_lpdf(y[n] | cur_mu[k], sigma[k]);
      target += log_sum_exp(lps);
  }
}
