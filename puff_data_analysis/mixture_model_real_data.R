library(dplyr)
library(readr)
library(mgcv)
library(tidyr)
library(rstan)
library(mvtnorm)

# helpful functions
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

### Import data

full_data <- read_csv("puff_experiment_data.csv") %>%
  mutate(intens_ratio = log(intens_ratio),
         snr = log(snr),
         smoothness = log(smoothness),
         cargo = as.character(cargo)) %>%
  group_by(cargo) %>%
  mutate(intens_ratio_scaled = scale(intens_ratio),
         conv_perim_scaled = scale(conv_perim),
         conv_area_scaled = scale(conv_area),
         noise_scaled = scale(noise),
         snr_scaled = scale(snr)) %>%
  ungroup()


# construct training and test data
train_data <- full_data %>%
  filter(cargo != "3",
         !(cell %in% c("1-1", "1-5",
                       "2-1", "2-3")))

validation_data <- full_data %>%
  filter(cell %in% c("1-1", "1-5",
                     "2-1", "2-3"))

test_data <- full_data %>%
  filter(cargo == "3")

# now fit the model
puff_gam <- train_data %>%
  gam(puff ~ s(conv_area_scaled) + s(conv_perim_scaled) +
        s(noise_scaled) + s(snr_scaled) +
        s(intens_ratio_scaled),
      data =., family=binomial())

test_data <- test_data %>%
  mutate(pred = predict.gam(puff_gam, test_data,
                            type = "response"))


# now label shift estimation
test_puff_est <- label_shift_estimate(train_data$puff, 
                                     puff_gam$fitted.values, 
                                     test_data$pred, 
                                     0.5)

w1 = test_puff_est/mean(train_data$puff)
w2 = (1-test_puff_est)/(1 - mean(train_data$puff))

# fit the model

full_mix <- stan_model('gauss_mix_hierarchical_assisted.stan')

input_data <- list(ncells = length(unique(test_data$cell)),
                   N = length(test_data$smoothness),
                   y = test_data$smoothness,
                   cell = as.numeric(as.factor(test_data$cell)),
                   theta = cbind(rep(test_puff_est, 7),
                                 rep(1 - test_puff_est, 7)))

best_value <- -Inf
for(attempt in 1:15){
  tryCatch({
    optimize <- optimizing(full_mix, data=input_data,
                           algorithm="LBFGS")
    
    if(optimize$value > best_value){
      best_value <- optimize$value
      best_pars <- optimize$par
    }
  }, error = function(e){})
}

puff_intercept = best_pars["mu[1]"]
nonpuff_intercept = best_pars["mu[2]"]
puff_re_sd = best_pars["re_sigma[1]"]
nonpuff_re_sd = best_pars["re_sigma[2]"]

test_data <- test_data %>%
  mutate(cell = as.numeric(as.factor(cell)),
         puff_est = best_pars[paste("mu_obs[", 
                                    cell,
                                    ",1]",
                                    sep="")],
         nonpuff_est = best_pars[paste("mu_obs[", 
                                       cell,
                                       ",2]",
                                       sep="")],
         puff_resid = puff_intercept + smoothness - puff_est,
         nonpuff_resid = nonpuff_intercept + smoothness - nonpuff_est)

# next, bootstrap the classifier
Zbar <- predict.gam(puff_gam, train_data,
                    type="lpmatrix")
lambda = puff_gam$sp[1]
S = matrix(0, nrow=ncol(Zbar), ncol=ncol(Zbar))
S[2:ncol(S), 2:ncol(S)] = puff_gam$smooth[[1]]$S[[1]]
W_Zbar = sweep(Zbar, 1, puff_gam$weights, "*")
beta_var = solve(t(Zbar) %*% W_Zbar + lambda*S)
beta_var = 0.5*(beta_var + t(beta_var))
betahat <- puff_gam$coefficients

# now bootstrap data
bootstrap_mean_est <- c()
num_boot <- 200
set.seed(212)

for(samp in 1:num_boot){
  tryCatch({
    bootstrap_test_data <- data.frame(conv_area_scaled = c(), 
                                      conv_perim_scaled = c(),
                                      noise_scaled = c(), 
                                      snr_scaled = c(),
                                      intens_ratio_scaled = c(),
                                      puff_resid = c(),
                                      nonpuff_resid = c(),
                                      cell = c(),
                                      puff = c())
    
    for(k in 1:7){
      re_puff <- rnorm(1, 0, puff_re_sd)
      re_nonpuff <- rnorm(1, 0, nonpuff_re_sd)
      
      bootstrap_test_data <- bootstrap_test_data %>%
        rbind(test_data %>%
                sample_n(7500, replace=T) %>%
                mutate(pred = pred*w1/(pred*w1 + (1 - pred)*w2),
                       puff = rbinom(7500, 1, pred)) %>%
                mutate(puff_resid = puff_resid + re_puff,
                       nonpuff_resid = nonpuff_resid + re_nonpuff,
                       smoothness = puff*puff_resid + (1-puff)*nonpuff_resid,
                       cell = k))
    }
    
    
    betasamp <- suppressWarnings(rmvnorm(1, mean=betahat, 
                                         sigma=beta_var))
    Zbar_new <- predict.gam(puff_gam,
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
    
    bootstrap_w1 = bootstrap_puff_est/mean(bootstrap_train_y)
    bootstrap_w2 = (1-bootstrap_puff_est)/(1 - mean(bootstrap_train_y))
    
    input_data <- list(ncells = length(unique(bootstrap_test_data$cell)),
                       N = length(bootstrap_test_data$smoothness),
                       y = bootstrap_test_data$smoothness,
                       cell = bootstrap_test_data$cell,
                       theta = cbind(rep(bootstrap_puff_est, 7),
                                     rep(1 - bootstrap_puff_est, 7)))
    
    best_bootstrap_value <- -Inf
    for(attempt in 1:15){
      tryCatch({
        optimize <- optimizing(full_mix, data=input_data)
        if(optimize$value > best_bootstrap_value){
          best_bootstrap_value <- optimize$value
          best_bootstrap_pars <- optimize$par
        }},
        error = function(e){})
      #print(attempt)
    }
    
    bootstrap_mean_est[samp] <- best_bootstrap_pars["mu[1]"]
    
    print(samp)
  }, error = function(e){})
  
}

output <- c(puff_intercept, bootstrap_mean_est)
write(output, 
      file="mixture_model_real_data_results.txt")