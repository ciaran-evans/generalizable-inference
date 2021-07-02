library(dplyr)
library(readr)
library(mgcv)
library(tidyr)
library(lme4)
library(mvtnorm)

# helpful functions
expit <- function(p){
  return(exp(p)/(1 + exp(p)))
}

label_shift_estimate_fixed_point <- function(a, b, 
                                             test_f, train_prop){
  possible_vals <- seq(a, b, by = 0.0001)
  diffs <- c()
  for(v in possible_vals){
    w1 <- v/train_prop
    w2 <- (1-v)/(1 - train_prop)
    corrected_preds <- w1*test_f/(w1*test_f + w2*(1 - test_f))
    diffs <- c(diffs, abs(mean(corrected_preds) - v))
  }
  return(possible_vals[which.min(diffs)])
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


### separate data into training, validation, and test

train_data <- full_data %>%
  filter(cargo != "3",
         !(cell %in% c("1-1", "1-5",
                       "2-1", "2-3")))

validation_data <- full_data %>%
  filter(cell %in% c("1-1", "1-5",
                     "2-1", "2-3"))

test_data <- full_data %>%
  filter(cargo == "3")

set.seed(2)

### fit logistic gam to training data, then 
### apply to validation and test data

puff_gam <- train_data %>%
  gam(puff ~ s(conv_area_scaled) + s(conv_perim_scaled) +
        s(noise_scaled) + s(snr_scaled) +
        s(intens_ratio_scaled),
      data =., family=binomial())

validation_data <- validation_data %>%
  mutate(pred = predict.gam(puff_gam, validation_data,
                            type = "response"))

test_data <- test_data %>%
  mutate(pred = predict.gam(puff_gam, test_data,
                            type = "response"))


### Label shift estimation

## Fixed point method
train_y  <- train_data$puff
test_f <- test_data$pred
test_puff_est_fp <- label_shift_estimate_fixed_point(0.001, 0.15, 
                                                     test_f, mean(train_y))

w1 = test_puff_est_fp/mean(train_data$puff)
w2 = (1-test_puff_est_fp)/(1 - mean(train_data$puff))



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
num_boot <- 200
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
  bootstrap_puff_est <- label_shift_estimate_fixed_point(0.001, 
                                                    0.15,
                                                    bootstrap_preds,
                                                    mean(bootstrap_train_y))
  
  bootstrap_mean_est[samp] <- bootstrap_puff_est
  
  #print(samp)
}


output <- c(test_puff_est_fp, bootstrap_mean_est)

write(output, 
      file="prevalence_real_data_results_fp.txt")
