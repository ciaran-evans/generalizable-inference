data {
  int<lower=2> ncells;
  int<lower=2> N;
  real y[N];               // observations
  int cell[N];
  simplex[2] theta[ncells];
}
parameters {
  ordered[2] mu;
  real diff_obs[ncells];
  vector<lower=0>[2] sigma;  // scales of mixture components
  //vector<lower=0>[2] re_sigma;
  real<lower=0.01> re_sigma;
}
model {
  vector[2] log_theta[ncells] = log(theta);  // cache log calculation
  for (k in 1:ncells){
    diff_obs[k] ~ normal(0, re_sigma);
  }
  
  for (n in 1:N) {
    int cur_cell = cell[n];
    vector[2] lps = log_theta[cur_cell];
    for (k in 1:2)
      lps[k] += normal_lpdf(y[n] | mu[k] + diff_obs[cur_cell], sigma[k]);
    target += log_sum_exp(lps);
  }
}
