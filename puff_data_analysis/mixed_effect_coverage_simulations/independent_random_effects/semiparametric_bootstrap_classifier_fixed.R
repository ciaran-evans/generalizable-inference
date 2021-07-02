library(dplyr)
library(readr)
library(mgcv)
library(tidyr)
library(lme4)
library(mvtnorm)

args <- commandArgs(trailingOnly = TRUE)
ci_rep = as.numeric(args[1])

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

simulate_training_data <- function(train_pi, 
                                   original_train_data){
  obs_y <- rbinom(nrow(original_train_data), 1, train_pi)
  train_data <- ((original_train_data %>%
                    filter(puff == 1) %>%
                    select(conv_area_scaled, conv_perim_scaled,
                           noise_scaled, snr_scaled,
                           intens_ratio_scaled) %>%
                    sample_n(nrow(original_train_data), replace=T))*matrix(rep(obs_y, 5), 
                                                                           ncol=5) +
                   (original_train_data %>%
                      filter(puff == 0) %>%
                      select(conv_area_scaled, conv_perim_scaled,
                             noise_scaled, snr_scaled,
                             intens_ratio_scaled) %>%
                      sample_n(nrow(original_train_data), replace=T))*matrix(rep(1-obs_y, 5), 
                                                                             ncol=5)) %>%
    as.data.frame() %>%
    mutate(puff = obs_y)
  
  return(train_data)
}

simulate_test_data <- function(original_test_data){
  test_data <- data.frame(conv_area_scaled = c(), 
                          conv_perim_scaled = c(),
                          noise_scaled = c(), 
                          snr_scaled = c(),
                          intens_ratio_scaled = c(),
                          smoothness = c(),
                          puff = c(),
                          cell = c())
  
  for(k in 1:7){
    pi_k <- runif(1, 0.02, 0.08)
    obs_y <- rbinom(7500, 1, pi_k)
    re_k <- rnorm(1, 0, 0.1)
    
    test_data <- test_data %>%
      rbind(((original_test_data %>%
                filter(puff == 1) %>%
                select(conv_area_scaled, conv_perim_scaled,
                       noise_scaled, snr_scaled,
                       intens_ratio_scaled, smoothness) %>%
                sample_n(7500, replace=T) %>%
                mutate(smoothness = smoothness + re_k))*matrix(rep(obs_y, 6), ncol=6) +
               (original_test_data %>%
                  filter(puff == 0) %>%
                  select(conv_area_scaled, conv_perim_scaled,
                         noise_scaled, snr_scaled,
                         intens_ratio_scaled, smoothness) %>%
                  sample_n(7500, replace=T) %>%
                  mutate(smoothness = smoothness + re_k))*matrix(rep(1-obs_y, 6), ncol=6)) %>%
              as.data.frame() %>%
              mutate(puff = obs_y,
                     cell = k))
  }
  return(test_data)
}


### Import data

full_data <- read_csv("../../puff_experiment_data.csv") %>%
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


# construct original training and test data
original_train_data <- full_data %>%
  filter(cargo != "3",
         !(cell %in% c("1-1", "1-5",
                       "2-1", "2-3")))

original_validation_data <- full_data %>%
  filter(cell %in% c("1-1", "1-5",
                     "2-1", "2-3"))

original_test_data <- full_data %>%
  filter(cargo == "3")


# bootstrap intervals
set.seed(ci_rep)
num_boot <- 200

# first create new training data
train_pi <- runif(1, 0.008, 0.015)
train_data <- simulate_training_data(train_pi, 
                                     original_train_data)


# now fit the model
puff_gam <- train_data %>%
  gam(puff ~ s(conv_area_scaled) + s(conv_perim_scaled) +
        s(noise_scaled) + s(snr_scaled) +
        s(intens_ratio_scaled),
      data =., family=binomial())

# and we'll create new test data
test_data <- simulate_test_data(original_test_data)

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


# now fit mixed effects models for this set of simulated data
me_mod_puff <- test_data %>%
  mutate(pred = pred*w1/(pred*w1 + (1 - pred)*w2)) %>%
  lmer(smoothness ~ (1|cell), weights = pred,
       data=.)

puff_intercept <- summary(me_mod_puff)$coefficients[1,1]
puff_re_sd <- sqrt(summary(me_mod_puff)$varcor[[1]][1])

me_mod_nonpuff <- test_data %>%
  mutate(pred = pred*w1/(pred*w1 + (1 - pred)*w2)) %>%
  lmer(smoothness ~ (1|cell), weights = 1-pred,
       data=.)

nonpuff_intercept <- summary(me_mod_nonpuff)$coefficients[1,1]
nonpuff_re_sd <- sqrt(summary(me_mod_nonpuff)$varcor[[1]][1])

### now we want to bootstrap the data 

# first, set up the residuals for bootstrapping
test_data <- test_data %>%
  mutate(puff_resid = puff_intercept + smoothness - 
           fitted(me_mod_puff),
         nonpuff_resid = nonpuff_intercept + smoothness - 
           fitted(me_mod_nonpuff))

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

# for each bootstrap sample:
bootstrap_mean_est <- c()
for(samp in 1:num_boot){
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
    re_k <- rnorm(1, 0, nonpuff_re_sd)
    
    bootstrap_test_data <- bootstrap_test_data %>%
      rbind(test_data %>%
              sample_n(7500, replace=T) %>%
              mutate(pred = pred*w1/(pred*w1 + (1 - pred)*w2),
                     puff = rbinom(7500, 1, pred)) %>%
              mutate(puff_resid = puff_resid + re_k,
                     nonpuff_resid = nonpuff_resid + re_k,
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
  
  bootstrap_me_mod <- bootstrap_test_data %>%
    mutate(pred = bootstrap_preds, 
           pred = pred*bootstrap_w1/(pred*bootstrap_w1 + 
                                       (1 - pred)*bootstrap_w2)) %>%
    lmer(smoothness ~ (1|cell), weights = pred,
         data=.)
  
  bootstrap_mean_est[samp] <- summary(bootstrap_me_mod)$coefficients[1,1]
}

output <- c(puff_intercept, bootstrap_mean_est)

write(output, 
      file=paste("semiparametric_bootstrap_fixed_",
                 ci_rep, ".txt", sep = ""))
