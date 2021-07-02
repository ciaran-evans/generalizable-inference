library(dplyr)
library(readr)
library(mgcv)
library(tidyr)
library(lme4)
library(mvtnorm)
library(sn)

expit <- function(p){
  return(exp(p)/(1 + exp(p)))
}

label_shift_estimate <- function(train_y, train_f, test_f, thresh){
  train_props <- c(1 - mean(train_y), mean(train_y))
  train_conf <- table(train_f > thresh, train_y)/length(train_f)
  test_pred_props <- c(1-mean(test_f > thresh), 
                       mean(test_f > thresh))
  w = solve(train_conf) %*% test_pred_props
  muy = train_props * w
  return(muy[2])
}


simulate_training_data <- function(train_pi, n_obs){
  obs_y <- rbinom(n_obs, 1, train_pi)
  train_data <- rsn(n_obs, -0.5, 2, 3)*(1-obs_y) + 
    (8 - rsn(n_obs, 3, 2, 3))*(obs_y)
  
  return(data.frame(z = train_data, y = obs_y))
}

simulate_test_data <- function(n_obs_cell, num_cells = 15){
  test_data <- data.frame(x = c(), z = c(), y = c(),
                          cell = c())
  for(k in 1:num_cells){
    pi_k <- 0.4
    obs_y <- rbinom(n_obs_cell, 1, pi_k)
    re_k_puff <- rnorm(1, 0, 0.5)
    re_k_nonpuff <- rnorm(1, 0, 0.2)
    
    z_k <- rsn(n_obs_cell, 0, 2, 3)*(1 - obs_y) + 
      (8 - rsn(n_obs_cell, 3, 2, 3))*(obs_y)
    x_k <- z_k + re_k_puff*obs_y +
      re_k_nonpuff*(1 - obs_y) + rnorm(n_obs_cell, 0, 0.2)
    test_data <- test_data %>%
      rbind(data.frame(x = x_k, z = z_k, y = obs_y, cell = k))
  }
  return(test_data)
}

train_pi <- 0.2
obs_per_cell = 500
num_ci = 300
num_boot = 200
prevalence_est <- c()
bootstrap_prevalence_est <- matrix(nrow = num_ci, ncol=num_boot)
num_cells = 15

set.seed(1)


for(ci_rep in 1:num_ci){
  tryCatch({
    train_data <- simulate_training_data(train_pi, 1000)
    test_data <- simulate_test_data(obs_per_cell, num_cells)
    
    train_gam <- train_data %>%
      gam(y ~ s(z), data = ., family=binomial())
    
    test_data <- test_data %>%
      mutate(pred = predict.gam(train_gam, test_data,
                                type = "response"))
    
    # now label shift estimation
    test_puff_est <- label_shift_estimate(train_data$y, 
                                          train_gam$fitted.values, 
                                          test_data$pred, 
                                          0.5)
    
    w1 = test_puff_est/mean(train_data$y)
    w2 = (1-test_puff_est)/(1 - mean(train_data$y))
    
    prevalence_est[ci_rep] <- test_puff_est
    
    
    # next, bootstrap the classifier
    Zbar <- predict.gam(train_gam, train_data,
                        type="lpmatrix")
    lambda = train_gam$sp[1]
    S = matrix(0, nrow=ncol(Zbar), ncol=ncol(Zbar))
    S[2:ncol(S), 2:ncol(S)] = train_gam$smooth[[1]]$S[[1]]
    W_Zbar = sweep(Zbar, 1, train_gam$weights, "*")
    beta_var = solve(t(Zbar) %*% W_Zbar + lambda*S)
    beta_var = 0.5*(beta_var + t(beta_var))
    betahat <- train_gam$coefficients
    
    # for each bootstrap sample:
    for(samp in 1:num_boot){
      bootstrap_test_data <- test_data %>%
        sample_n(nrow(test_data), replace=T) %>%
        mutate(pred = pred*w1/(pred*w1 + (1 - pred)*w2),
               puff = rbinom(nrow(test_data), 1, pred))
      
      betasamp <- suppressWarnings(rmvnorm(1, mean=betahat, 
                                           sigma=beta_var))
      Zbar_new <- predict.gam(train_gam,
                              newdata = bootstrap_test_data,
                              type="lpmatrix")
      
      bootstrap_preds <- expit(Zbar_new %*% t(betasamp))
      
      bootstrap_train_f <- expit(Zbar %*% t(betasamp))
      bootstrap_train_y  <- rbinom(nrow(train_data), 1, 
                                   bootstrap_train_f)
      bootstrap_prevalence_est[ci_rep, samp] <- label_shift_estimate(bootstrap_train_y, 
                                                 bootstrap_train_f, 
                                                 bootstrap_preds, 
                                                 0.5)
      
      print(samp)
    }
  }, error = function(e){})
  
  print(ci_rep)
}


write.table(prevalence_est, file="prevalence_scenario_2c_estimates.txt")
write.table(bootstrap_prevalence_est, file="prevalence_scenario_2c_ci_reps.txt")
