library(dplyr)
library(readr)
library(ggplot2)
library(mgcv)
library(tidyr)
library(gridExtra)
library(xtable)
library(rstan)

### Helper functions

logit <- function(p){
  return(log(p/(1 - p)))
}

expit <- function(p){
  return(exp(p)/(1 + exp(p)))
}

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# label shift estimation from Lipton et al. (2018)
label_shift_estimate_lipton <- function(train_y, train_f, test_f, thresh){
  train_props <- c(1 - mean(train_y), mean(train_y))
  train_conf <- table(train_f > thresh, train_y)/length(train_f)
  test_pred_props <- c(1-mean(test_f > thresh), 
                       mean(test_f > thresh))
  w = solve(train_conf) %*% test_pred_props
  muy = train_props * w
  return(muy[2])
}

# label shift estimation with fixed-point method
label_shift_estimate_fixed_point <- function(a, b, test_f, train_prop){
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


################################################
###### Data, classifier, and assumptions #######
################################################

### Table 3 (puffs by cell and cargo)
full_data %>%
  group_by(cargo, cell) %>%
  summarize(puff_num = sum(puff),
            puff_prop = mean(puff)) %>%
  xtable::xtable(digits=3)


### Figure 3 (feature distributions)
p1 <- full_data %>%
  mutate(puff = ifelse(puff == 1, "Puff", "Nonpuff")) %>%
  ggplot(aes(x = conv_perim_scaled, color=cargo, 
             lty = as.factor(puff))) +
  geom_density(lwd = 1) +
  scale_color_manual(values=c("#fcbe4bff", "#748cc2ff", "#d64f96ff")) +
  theme_bw() +
  labs(x = "ConvexPerimeter",
       y = "Density", color = "Condition", lty="") +
  theme(legend.position="none")

p2 <- full_data %>%
  mutate(puff = ifelse(puff == 1, "Puff", "Nonpuff")) %>%
  ggplot(aes(x = noise_scaled, color=cargo, 
             lty = as.factor(puff))) +
  geom_density(lwd = 1) +
  scale_color_manual(values=c("#fcbe4bff", "#748cc2ff", "#d64f96ff")) +
  theme_bw() +
  labs(x = "Noise",
       y = "Density", color = "Condition", lty="")  +
  theme(legend.position="none")

p3 <- full_data %>%
  mutate(puff = ifelse(puff == 1, "Puff", "Nonpuff")) %>%
  ggplot(aes(x = smoothness, color=cargo, 
             lty = as.factor(puff))) +
  geom_density(lwd = 1) +
  scale_color_manual(values=c("#fcbe4bff", "#748cc2ff", "#d64f96ff")) +
  theme_bw() +
  labs(x = "Smoothness",
       y = "Density", color = "Condition", lty="") 

legend <- get_legend(p3)

p3 <- p3 + theme(legend.position="none")

pdf(file="feature_distributions.pdf",
    width=16, height=4)
grid.arrange(p1, p2, p3, legend, ncol=4, widths = c(3, 3, 3, 1))
dev.off()



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

## Method from Lipton et al. (2018)
train_y  <- train_data$puff
train_f <- puff_gam$fitted.values
thresh <- 0.5

test_f <- test_data$pred
test_puff_est_lipton <- label_shift_estimate_lipton(train_y, train_f, 
                                             test_f, thresh)

## Fixed point method
test_puff_est_fp <- label_shift_estimate_fixed_point(0.001, 0.15, 
                                                     test_f, mean(train_y))


## True test prevalence
mean(test_data$puff)


### Figure 4 (calibration plots)
w1 = test_puff_est_fp/mean(train_data$puff)
w2 = (1-test_puff_est_fp)/(1 - mean(train_data$puff))

p1 <- validation_data %>%
  mutate(bin = .bincode(pred, 
                        breaks=c(-1, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.1))) %>%
  group_by(bin) %>%
  summarize(mean_prob = mean(pred),
            true_prop = mean(puff),
            puff_count = sum(puff),
            event_count = n(),
            lower = binom.test(puff_count, event_count, 
                               mean_prob, conf.level=0.99)$conf.int[1],
            upper = binom.test(puff_count, event_count, 
                               mean_prob, conf.level=0.99)$conf.int[2]) %>%
  ggplot(aes(x = mean_prob, y = true_prop)) +
  geom_point() +
  geom_errorbar(aes(ymin=lower, ymax=upper), color="orange") +
  geom_abline(slope=1, intercept=0, lty=2, color="red") +
  scale_y_continuous(limits=c(0,1)) +
  scale_x_continuous(limits=c(0,1)) +
  theme_bw() +
  labs(x = "Bin mean of classifier predictions", 
       y = "Proportion of events that are puffs",
       title = "Calibration plot for validation data (Conditions 1 and 2))")

p2 <- test_data %>%
  mutate(bin = .bincode(pred, 
                        breaks=c(-1, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.1))) %>%
  group_by(bin) %>%
  summarize(mean_prob = mean(pred),
            true_prop = mean(puff),
            puff_count = sum(puff),
            event_count = n(),
            lower = binom.test(puff_count, event_count, 
                               mean_prob, conf.level=0.99)$conf.int[1],
            upper = binom.test(puff_count, event_count, 
                               mean_prob, conf.level=0.99)$conf.int[2]) %>%
  ggplot(aes(x = mean_prob, y = true_prop)) +
  geom_point() +
  geom_errorbar(aes(ymin=lower, ymax=upper), color="orange") +
  geom_abline(slope=1, intercept=0, lty=2, color="red") +
  scale_y_continuous(limits=c(0,1)) +
  scale_x_continuous(limits=c(0,1)) +
  theme_bw() +
  labs(x = "Bin mean of classifier predictions", 
       y = "Proportion of events that are puffs",
       title = "Calibration plot for test data (Condition 3)")

p3 <- test_data %>%
  mutate(pred = pred*w1/(pred*w1 + (1 - pred)*w2)) %>%
  mutate(bin = .bincode(pred, 
                        breaks=c(-1, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.1))) %>%
  group_by(bin) %>%
  summarize(mean_prob = mean(pred),
            true_prop = mean(puff),
            puff_count = sum(puff),
            event_count = n(),
            lower = binom.test(puff_count, event_count, 
                               mean_prob, conf.level=0.99)$conf.int[1],
            upper = binom.test(puff_count, event_count, 
                               mean_prob, conf.level=0.99)$conf.int[2]) %>%
  ggplot(aes(x = mean_prob, y = true_prop)) +
  geom_point() +
  geom_errorbar(aes(ymin=lower, ymax=upper), color="orange") +
  geom_abline(slope=1, intercept=0, lty=2, color="red") +
  scale_y_continuous(limits=c(0,1)) +
  scale_x_continuous(limits=c(0,1)) +
  theme_bw() +
  labs(x = "Bin mean of classifier predictions", 
       y = "Proportion of events that are puffs",
       title = "Calibration plot for test data (Condition 3), label shift adjustment")

pdf(file = "classifier_calibration_plots.pdf",
    width = 18, height=5)
grid.arrange(p1, p2, p3, ncol=3)
dev.off()



### Figure 5 (checking A3)
full_test_gam <- test_data %>%
  gam(puff ~ s(conv_area_scaled) + s(conv_perim_scaled) +
        s(noise_scaled) + s(snr_scaled) +
        s(intens_ratio_scaled) + s(smoothness),
      data =., family=binomial())

p1 <- test_data %>%
  mutate(pred = pred*w1/(pred*w1 + (1 - pred)*w2),
         true_est = full_test_gam$fitted.values) %>%
  ggplot(aes(x = pred, y = true_est)) +
  geom_point(alpha = 0.3) +
  geom_abline(slope=1, intercept=0, color="red", lwd=1) +
  theme_bw() +
  labs(x = "Classifier predicted probability, with label shift correction",
       y = "Estimate of true probability from test data",
       title = "Assessment of classifier predictions on test data (Condition 3)")

p2 <- full_data %>%
  group_by(cargo, cell) %>%
  summarize(puff_mean = mean(smoothness[puff == 1]),
            nonpuff_mean = mean(smoothness[puff == 0])) %>%
  ggplot(aes(x = puff_mean, y = nonpuff_mean, color = cargo)) +
  geom_point(size=3) + 
  scale_color_manual(values=c("#fcbe4bff", "#748cc2ff", "#d64f96ff")) +
  theme_bw() +
  labs(x = "Mean Smoothness for puffs in each cell",
       y = "Mean Smoothness for nonpuffs in each cell",
       color = "Condition",
       title = "Relationship between class-conditional means")

pdf(file = "classifier_prediction_assessment.pdf",
    width=16, height=6)
grid.arrange(p1, p2, ncol=2)
dev.off()



##################################################
###### Inference with classifier predictions #####
##################################################

### Table 4
test_data %>%
  group_by(cell) %>%
  rename(feature = smoothness) %>%
  mutate(pred2 = pred*w1/(pred*w1 + (1 - pred)*w2)) %>%
  summarize(puff_mean = mean(feature[puff == 1]),
            weighted_mean_raw = sum(feature * pred)/sum(pred),
            weighted_mean_label_shift = sum(feature * pred2)/sum(pred2),
            thresh_mean_raw = mean(feature[pred > 0.5]),
            thresh_mean_label_shift = mean(feature[pred2 > 0.5]),
            nonpuff_mean = mean(feature[puff == 0]))

test_data %>%
  group_by(cell) %>%
  rename(feature = smoothness) %>%
  mutate(pred2 = pred*w1/(pred*w1 + (1 - pred)*w2)) %>%
  summarize(puff_mean = mean(feature[puff == 1]),
            weighted_mean_raw = sum(feature * pred)/sum(pred),
            weighted_mean_label_shift = sum(feature * pred2)/sum(pred2),
            thresh_mean_raw = mean(feature[pred > 0.5]),
            thresh_mean_label_shift = mean(feature[pred2 > 0.5]),
            nonpuff_mean = mean(feature[puff == 0])) %>%
  ungroup() %>%
  summarize(puff_mean = mean(puff_mean),
            weighted_mean_raw = mean(weighted_mean_raw),
            weighted_mean_label_shift = mean(weighted_mean_label_shift),
            thresh_mean_raw = mean(thresh_mean_raw),
            thresh_mean_label_shift = mean(thresh_mean_label_shift),
            nonpuff_mean = mean(nonpuff_mean))


#################################################
### Inference for prevalence                  ###
### (95% confidence intervals for prevalence) ###
#################################################

## Using Lipton's label shift correction method
## (prevalence_real_data_results_lipton.txt is generated by 
##  prevalence_ci_real_data_lipton.R)
lipton_results <- read_delim("prevalence_real_data_results_lipton.txt", 
                         " ", escape_double = FALSE, col_names = FALSE, 
                         trim_ws = TRUE) %>%
  unlist()
lipton_results <- lipton_results[!is.na(lipton_results)]
lower <- 2*lipton_results[1] - quantile(lipton_results[2:201],
                                    0.975, na.rm=T)
upper <- 2*lipton_results[1] - quantile(lipton_results[2:201],
                                    0.025, na.rm=T)

## Using fixed point method
## (prevalence_real_data_results_fp.txt is generated by 
##  prevalence_ci_real_data_fixed_point.R)
fp_results <- read_delim("prevalence_real_data_results_fp.txt", 
           " ", escape_double = FALSE, col_names = FALSE, 
           trim_ws = TRUE) %>%
  unlist()
fp_results <- fp_results[!is.na(fp_results)]
lower <- 2*fp_results[1] - quantile(fp_results[2:201],
                                         0.975, na.rm=T)
upper <- 2*fp_results[1] - quantile(fp_results[2:201],
                                         0.025, na.rm=T)

### Simulations to assess coverage

## Using Lipton method
prevalence_ests <- c()
bootstrap_sim_results <- matrix(nrow=300, ncol=200)
for(r in 1:300){
  tryCatch({
    fname = paste("prevalence_coverage_simulations/lipton_coverage/prevalence_bootstrap_",
                  r, ".txt", sep="")
    temp <- read_delim(fname, 
                       " ", escape_double = FALSE, col_names = FALSE, 
                       trim_ws = TRUE) %>%
      unlist()
    temp <- temp[!is.na(temp)]
    bootstrap_sim_results[r,] <- temp[3:202]
    prevalence_ests[r] <- temp[1]
    print(r)
  }, error = function(e){})
  
}

covers <- c()
for(r in 1:300){
  lower <- 2*prevalence_ests[r] - quantile(bootstrap_sim_results[r,],
                                           0.975, na.rm=T)
  upper <- 2*prevalence_ests[r] - quantile(bootstrap_sim_results[r,],
                                           0.025, na.rm=T)
  covers <- c(covers,
              lower <= 0.05 & upper >= 0.05)
}
mean(covers, na.rm=T)


## Using fixed point method
prevalence_ests <- c()
bootstrap_sim_results <- matrix(nrow=300, ncol=200)
for(r in 1:300){
  tryCatch({
    fname = paste("prevalence_coverage_simulations/fixed_point_coverage/prevalence_bootstrap_",
                  r, ".txt", sep="")
    temp <- read_delim(fname, 
                       " ", escape_double = FALSE, col_names = FALSE, 
                       trim_ws = TRUE) %>%
      unlist()
    temp <- temp[!is.na(temp)]
    bootstrap_sim_results[r,] <- temp[3:202]
    prevalence_ests[r] <- temp[1]
    print(r)
  }, error = function(e){})
  
}

covers <- c()
for(r in 1:300){
  lower <- 2*prevalence_ests[r] - quantile(bootstrap_sim_results[r,],
                                           0.975, na.rm=T)
  upper <- 2*prevalence_ests[r] - quantile(bootstrap_sim_results[r,],
                                           0.025, na.rm=T)
  covers <- c(covers,
              lower <= 0.05 & upper >= 0.05)
}
mean(covers, na.rm=T)



#################################################
###### Inference for feature distributions   ####
###### (mixed effects models for Smoothness  ####
######  with classifier probability weights) ####
#################################################

### Label dependent random effects
## (the confidence interval is created by the file 
##  me_model_real_data_label_dependent.R)
label_dependent_results <- read_delim("me_model_real_data_results_label_dependent.txt", 
                         " ", escape_double = FALSE, col_names = FALSE, 
                         trim_ws = TRUE) %>%
  unlist()
label_dependent_results <- label_dependent_results[!is.na(label_dependent_results)]
lower <- 2*label_dependent_results[1] - quantile(label_dependent_results[2:201],
                                    0.975, na.rm=T)
upper <- 2*label_dependent_results[1] - quantile(label_dependent_results[2:201],
                                    0.025, na.rm=T)

## Simulation to assess coverage for label dependent approach


### Label independent random effects
## (the confidence interval is created by the file 
##  me_model_real_data_label_independent.R)
label_independent_results <- read_delim("me_model_real_data_results_label_independent.txt", 
                                      " ", escape_double = FALSE, col_names = FALSE, 
                                      trim_ws = TRUE) %>%
  unlist()
label_independent_results <- label_independent_results[!is.na(label_independent_results)]

label_independent_results[1]

lower <- 2*label_independent_results[1] - quantile(label_independent_results[2:201],
                                                 0.975, na.rm=T)
upper <- 2*label_independent_results[1] - quantile(label_independent_results[2:201],
                                                 0.025, na.rm=T)

## Simulation to assess coverage for label independent approach

bootstrap_ests <- c()
bootstrap_sim_results <- matrix(nrow=300, ncol=200)
for(r in 1:300){
  tryCatch({
    fname = paste("mixed_effect_coverage_simulations/independent_random_effects/semiparametric_bootstrap_fixed_",
                  r, ".txt", sep="")
    temp <- read_delim(fname, 
                       " ", escape_double = FALSE, col_names = FALSE, 
                       trim_ws = TRUE) %>%
      unlist()
    temp <- temp[!is.na(temp)]
    bootstrap_sim_results[r,] <- temp[2:201]
    bootstrap_ests[r] <- temp[1]
    print(r)
  }, error = function(e){})
  
}

covers <- c()
for(r in 1:300){
  lower <- 2*bootstrap_ests[r] - quantile(bootstrap_sim_results[r,],
                                          0.975, na.rm=T)
  upper <- 2*bootstrap_ests[r] - quantile(bootstrap_sim_results[r,],
                                          0.025, na.rm=T)
  covers <- c(covers,
              lower <= -2.86 & upper >= -2.86)
}
mean(covers, na.rm=T)


##########################################################
###### Mixture models                           ##########
###### (inference with gaussian mixture models) ##########
##########################################################

### Figure 8 (appendix)
# (the mixture model struggles if we 
# don't use good mixing proportions)
full_mix <- stan_model('gauss_mix_hierarchical_unassisted.stan')

input_data <- list(ncells = length(unique(test_data$cell)),
                   N = length(test_data$smoothness),
                   y = test_data$smoothness,
                   cell = as.numeric(as.factor(test_data$cell)))

set.seed(4)
optimize <- optimizing(full_mix, data=input_data)

eval_pts <- seq(-4, -1, 0.1)
fitted_vals_puff <- data.frame(x = c(), cell = c(), y = c())
fitted_vals_nonpuff <- data.frame(x = c(), cell = c(), y = c())
for(c in unique(test_data$cell)){
  ind = which(levels(as.factor(test_data$cell)) == c)
  fitted_vals_puff <- rbind(fitted_vals_puff,
                            data.frame(x = eval_pts,
                                       cell = as.character(c),
                                       y = dnorm(eval_pts,
                                                 optimize$par[paste("mu_obs[", 
                                                                    ind,
                                                                    ",1]", sep="")],
                                                 optimize$par["sigma[1]"])))
  fitted_vals_nonpuff <- rbind(fitted_vals_nonpuff,
                               data.frame(x = eval_pts,
                                          cell = as.character(c),
                                          y = dnorm(eval_pts,
                                                    optimize$par[paste("mu_obs[", 
                                                                       ind,
                                                                       ",2]", sep="")],
                                                    optimize$par["sigma[2]"])))
}

fitted_vals <- fitted_vals_puff %>%
  mutate(puff = "Puff") %>%
  rbind(fitted_vals_nonpuff %>%
          mutate(puff = "Nonpuff"))


pdf(file="bad_mixture_estimates.pdf",
    width=16, height=4)
test_data %>%
  mutate(puff = ifelse(puff == 1, "Puff", "Nonpuff")) %>%
  ggplot() +
  geom_density(aes(x = smoothness, 
                   color=as.factor(puff)), lwd = 1) +
  geom_line(aes(x = x, y = y, lty = puff), color="black", data = fitted_vals) +
  facet_wrap(~cell, ncol=4) +
  theme_bw() + 
  labs(x = "Smoothness", y = "Density",
       color = "True labels",
       lty = "Estimated distributions",
       title = "Mixture model estimation with test data (Condition 3)")
dev.off()


### Figure 9 (appendix)
# (we can noticeably improve the fit 
#  if we use the correct mixing proportions)
full_mix <- stan_model('gauss_mix_hierarchical_assisted.stan')

input_data <- list(ncells = length(unique(test_data$cell)),
                   N = length(test_data$smoothness),
                   y = test_data$smoothness,
                   cell = as.numeric(as.factor(test_data$cell)),
                   theta = cbind(rep(test_puff_est_lipton, 7),
                                 rep(1 - test_puff_est_lipton, 7)))

set.seed(7)
optimize <- optimizing(full_mix, data=input_data)
optimize$par
fitted_vals <- dnorm(test_data$smoothness,
                     optimize$par[paste("mu_obs[", 
                                        as.numeric(as.factor(test_data$cell)),
                                        ",1]", sep="")],
                     optimize$par["sigma[1]"])*test_data$puff + 
  dnorm(test_data$smoothness,
        optimize$par[paste("mu_obs[", 
                           as.numeric(as.factor(test_data$cell)),
                           ",2]", sep="")],
        optimize$par["sigma[2]"])*(1-test_data$puff)


pdf(file="good_mixture_estimates.pdf",
    width=16, height=4)
test_data %>%
  mutate(puff = ifelse(puff == 1, "Puff", "Nonpuff"),
         cell_numeric = as.numeric(as.factor(cell)),
         norm_est = fitted_vals) %>%
  ggplot(aes(x = smoothness)) +
  geom_density(aes(color=as.factor(puff)), lwd = 1) +
  geom_line(aes(y = norm_est, lty=as.factor(puff)), color="black") +
  facet_wrap(~cell, ncol=4) +
  theme_bw() + 
  labs(x = "Smoothness", y = "Density",
       color = "True labels",
       lty = "Estimated distributions",
       title = "Mixture model estimation with test data (Condition 3)")
dev.off()


### Table 6 (appendix)
test_data %>%
  mutate(cell_numeric = as.numeric(as.factor(cell))) %>%
  group_by(cell, cell_numeric) %>%
  summarize(puff_mean = mean(smoothness[puff == 1]),
            nonpuff_mean = mean(smoothness[puff == 0])) %>%
  left_join(data.frame(cell_numeric = 1:7,
                       mixture_puff_mean = c(optimize$par["mu_obs[1,1]"],
                                             optimize$par["mu_obs[2,1]"],
                                             optimize$par["mu_obs[3,1]"],
                                             optimize$par["mu_obs[4,1]"],
                                             optimize$par["mu_obs[5,1]"],
                                             optimize$par["mu_obs[6,1]"],
                                             optimize$par["mu_obs[7,1]"]),
                       mixture_nonpuff_mean = c(optimize$par["mu_obs[1,2]"],
                                                optimize$par["mu_obs[2,2]"],
                                                optimize$par["mu_obs[3,2]"],
                                                optimize$par["mu_obs[4,2]"],
                                                optimize$par["mu_obs[5,2]"],
                                                optimize$par["mu_obs[6,2]"],
                                                optimize$par["mu_obs[7,2]"])),
            by = c("cell_numeric"))


### Confidence interval from mixture model
## (the confidence interval is created by the file 
##  mixture_model_real_data.R)
mixture_results <- read_delim("mixture_model_real_data_results.txt", 
                                        " ", escape_double = FALSE, col_names = FALSE, 
                                        trim_ws = TRUE) %>%
  unlist()
mixture_results <- mixture_results[!is.na(mixture_results)]
lower <- 2*mixture_results[1] - quantile(mixture_results[2:201],
                                                   0.975, na.rm=T)
upper <- 2*mixture_results[1] - quantile(mixture_results[2:201],
                                                   0.025, na.rm=T)

### Simulation to assess mixture model coverage

bootstrap_ests <- c()
bootstrap_sim_results <- matrix(nrow=300, ncol=200)
for(r in 1:300){
  fname = paste("mixture_model_coverage_simulations/mixture_model_bootstrap_",
                r, ".txt", sep="")
  temp <- read_delim(fname, 
                     " ", escape_double = FALSE, col_names = FALSE, 
                     trim_ws = TRUE) %>%
    unlist()
  temp <- temp[!is.na(temp)]
  bootstrap_sim_results[r,] <- temp[2:201]
  bootstrap_ests[r] <- temp[1]
  print(r)
}

covers <- c()
lowers <- c()
uppers <- c()
for(r in 1:300){
  lower <- 2*bootstrap_ests[r] - quantile(bootstrap_sim_results[r,],
                                          0.975, na.rm=T)
  upper <- 2*bootstrap_ests[r] - quantile(bootstrap_sim_results[r,],
                                          0.025, na.rm=T)
  covers <- c(covers,
              lower <= -2.86 & upper >= -2.86)
  lowers[r] <- lower
  uppers[r] <- upper
  #print(c(lower, upper))
}
mean(covers)
