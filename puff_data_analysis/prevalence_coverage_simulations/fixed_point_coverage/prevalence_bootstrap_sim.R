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

label_shift_estimate <- function(train_y, test_f){
  possible_vals <- seq(0.001, 0.15, by = 0.0001)
  possible_means <- c()
  for(val in possible_vals){
    temp_w1 <- val/mean(train_y)
    temp_w2 <- (1 - val)/(1 - mean(train_y))
    temp_corr <- test_f*temp_w1/(test_f*temp_w1 + 
                                   (1-test_f)*temp_w2)
    possible_means <- c(possible_means, mean(temp_corr))
  }
  prev_est <- possible_vals[which.min(abs(possible_means - possible_vals))]
  prev_est
}

simulate_data <- function(puff_prop, 
                          original_data){
  obs_y <- rbinom(nrow(original_data), 1, puff_prop)
  new_data <- ((original_data %>%
                  filter(puff == 1) %>%
                  select(conv_area_scaled, conv_perim_scaled,
                         noise_scaled, snr_scaled,
                         intens_ratio_scaled) %>%
                  sample_n(nrow(original_data), replace=T))*matrix(rep(obs_y, 5), 
                                                                   ncol=5) +
                 (original_data %>%
                    filter(puff == 0) %>%
                    select(conv_area_scaled, conv_perim_scaled,
                           noise_scaled, snr_scaled,
                           intens_ratio_scaled) %>%
                    sample_n(nrow(original_data), replace=T))*matrix(rep(1-obs_y, 5), 
                                                                     ncol=5)) %>%
    as.data.frame() %>%
    mutate(puff = obs_y)
  
  return(new_data)
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
train_pi <- 0.01
train_data <- simulate_data(train_pi, original_train_data)


# now fit the model
puff_gam <- train_data %>%
  gam(puff ~ s(conv_area_scaled) + s(conv_perim_scaled) +
        s(noise_scaled) + s(snr_scaled) +
        s(intens_ratio_scaled),
      data =., family=binomial())

# and we'll create new test data
test_pi <- 0.05
test_data <- simulate_data(test_pi, original_test_data)

test_data <- test_data %>%
  mutate(pred = predict.gam(puff_gam, test_data,
                            type = "response"))

# now label shift estimation
test_puff_est <- label_shift_estimate(train_data$puff,
                                     test_data$pred)
#test_puff_est
bad_est <- mean(test_data$pred)

w1 = test_puff_est/mean(train_data$puff)
w2 = (1-test_puff_est)/(1 - mean(train_data$puff))


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
  bootstrap_test_data <- test_data %>%
    sample_n(nrow(test_data), replace=T) %>%
    mutate(pred = pred*w1/(pred*w1 + (1 - pred)*w2),
           puff = rbinom(nrow(test_data), 1, pred))
  
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
                                             bootstrap_preds)
  
  bootstrap_mean_est[samp] <- bootstrap_puff_est
  
  #print(samp)
}

output <- c(test_puff_est, bad_est, bootstrap_mean_est)

write(output, 
      file=paste("prevalence_bootstrap_",
                 ci_rep, ".txt", sep = ""))
