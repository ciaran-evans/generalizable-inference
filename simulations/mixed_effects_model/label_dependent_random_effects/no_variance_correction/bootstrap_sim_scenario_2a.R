library(dplyr)
library(readr)
library(mgcv)
library(tidyr)
library(lme4)
library(mvtnorm)

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
  train_data <- rnorm(n_obs, -0.5, 1)*obs_y + 
    rnorm(n_obs, 3, 1)*(1-obs_y)
  
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
    
    z_k <- rnorm(n_obs_cell, 0, 1)*obs_y + 
      rnorm(n_obs_cell, 3, 1)*(1-obs_y)
    x_k <- 3 - z_k + re_k_puff*obs_y +
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
test_mean_est <- c()
bootstrap_mean_est <- matrix(nrow = num_ci, ncol=num_boot)
num_cells = 15

set.seed(2)

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
    
    me_mod_puff <- test_data %>%
      mutate(pred = pred*w1/(pred*w1 + (1 - pred)*w2)) %>%
      lmer(x ~ (1|cell), weights = pred,
           data=.)
    
    puff_intercept <- summary(me_mod_puff)$coefficients[1,1]
    puff_re_sd <- sqrt(summary(me_mod_puff)$varcor[[1]][1])
    
    test_mean_est[ci_rep] <- puff_intercept
    
    me_mod_nonpuff <- test_data %>%
      mutate(pred = pred*w1/(pred*w1 + (1 - pred)*w2)) %>%
      lmer(x ~ (1|cell), weights = 1-pred,
           data=.)
    
    nonpuff_intercept <- summary(me_mod_nonpuff)$coefficients[1,1]
    nonpuff_re_sd <- sqrt(summary(me_mod_nonpuff)$varcor[[1]][1])
    
    
    test_data <- test_data %>%
      mutate(puff_resid = puff_intercept + x - 
               fitted(me_mod_puff),
             nonpuff_resid = nonpuff_intercept + x - 
               fitted(me_mod_nonpuff))
    
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
      bootstrap_test_data <- data.frame(x = c(),
                                        z = c(),
                                        y = c(),
                                        cell = c(),
                                        pred = c())
      
      for(k in 1:num_cells){
        re_puff <- rnorm(1, 0, puff_re_sd)
        re_nonpuff <- rnorm(1, 0, nonpuff_re_sd)
        
        bootstrap_test_data <- bootstrap_test_data %>%
          rbind(test_data %>%
                  sample_n(obs_per_cell, replace=T) %>%
                  mutate(pred = pred*w1/(pred*w1 + (1 - pred)*w2),
                         puff = rbinom(obs_per_cell, 1, pred)) %>%
                  mutate(puff_resid = puff_resid + re_puff,
                         nonpuff_resid = nonpuff_resid + re_nonpuff,
                         resid_dist = puff*puff_resid + (1-puff)*nonpuff_resid,
                         cell = k))
      }
      
      betasamp <- suppressWarnings(rmvnorm(1, mean=betahat, 
                                           sigma=beta_var))
      Zbar_new <- predict.gam(train_gam,
                              newdata = bootstrap_test_data,
                              type="lpmatrix")
      
      bootstrap_preds <- expit(Zbar_new %*% t(betasamp))
      
      bootstrap_train_f <- expit(Zbar %*% t(betasamp))
      bootstrap_train_y  <- rbinom(nrow(train_data), 1, 
                                   bootstrap_train_f)
      bootstrap_puff_est <- label_shift_estimate(bootstrap_train_y, 
                                                 bootstrap_train_f, 
                                                 bootstrap_preds, 
                                                 0.5)
      
      if(!is.na(bootstrap_puff_est) & bootstrap_puff_est > 0 & bootstrap_puff_est < 1){
        bootstrap_w1 = bootstrap_puff_est/mean(bootstrap_train_y)
        bootstrap_w2 = (1-bootstrap_puff_est)/(1 - mean(bootstrap_train_y))
        
        bootstrap_me_mod <- bootstrap_test_data %>%
          mutate(pred = bootstrap_preds, 
                 pred = pred*bootstrap_w1/(pred*bootstrap_w1 + 
                                             (1 - pred)*bootstrap_w2)) %>%
          lmer(resid_dist ~ (1|cell), weights = pred,
               data=.)
        
        bootstrap_mean_est[ci_rep, samp] <- summary(bootstrap_me_mod)$coefficients[1,1]
        
        #print(samp)
      }
    }
  }, error = function(e){})
  
  print(ci_rep)
}

write.table(test_mean_est, file="ls_bootstrap_scenario_2a_estimates.txt")
write.table(bootstrap_mean_est, file="ls_bootstrap_scenario_2a_ci_reps.txt")
