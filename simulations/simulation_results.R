### Scenarios (see Simulations section in the paper)
### There are 3 main scenarios (numbered 1, 2, 3)
###   Scenario 1: (A1), (A2), and (A3) are satisfied
###   Scenario 2: (A2) and (A3) are satisfied
###   Scenario 3: (A1) and (A2) are satisfied
### For each scenario, data either come from a normal or 
### skewed normal distribution. a denotes normal, c denotes
### skewed normal.
### 
### Data can also be generated with label-independent 
### random effects (Table 1 in the paper) or label-dependent 
### random effects (Table 2 in the paper).


########################################
######## Prevalence   ##################
########################################

### Results are in files 
### prevalence_scenario_xx_estimates.txt (point estimates at 
###  each repetition of the simulation) and 
### prevalence_scenario_xx_ci_reps.txt (bootstrap replications 
###  at each simulation repetition, to construct CIs)
### These files are created by 
### prevalence_scenario_xx.R

prevalence_coverage <- vector(length=6)
prevalence_coverage_sd <- vector(length=6)
prevalence_means <- vector(length=6)
prevalence_sd <- vector(length=6)

prevalence_1a_ests <- read_delim("prevalence/prevalence_scenario_1a_estimates.txt", 
                                 " ", escape_double = FALSE, col_names = FALSE, 
                                 trim_ws = TRUE, skip = 1)$X2
prevalence_1a_bootstrap <- read_delim("prevalence/prevalence_scenario_1a_ci_reps.txt", 
                                      " ", escape_double = FALSE, col_names = FALSE, 
                                      trim_ws = TRUE, skip = 1)[,-1]

covers <- c()
for(r in 1:300){
  lower <- 2*prevalence_1a_ests[r] - quantile(prevalence_1a_bootstrap[r,],
                                              0.975, na.rm=T)
  upper <- 2*prevalence_1a_ests[r] - quantile(prevalence_1a_bootstrap[r,],
                                              0.025, na.rm=T)
  covers[r] <- lower <= 0.4 & upper >= 0.4
  #print(r)
}

prevalence_coverage[1] <- mean(covers, na.rm=T)
prevalence_coverage_sd[1] <- sqrt(mean(covers, 
                                       na.rm=T)*(1-mean(covers, 
                                                        na.rm=T))/sum(!is.na(covers)))
prevalence_means[1] <- mean(prevalence_1a_ests, na.rm=T)
prevalence_sd[1] <- sd(prevalence_1a_ests, na.rm=T)/sqrt(sum(!is.na(prevalence_1a_ests)))




prevalence_1c_ests <- read_delim("prevalence/prevalence_scenario_1c_estimates.txt", 
                                 " ", escape_double = FALSE, col_names = FALSE, 
                                 trim_ws = TRUE, skip = 1)$X2
prevalence_1c_bootstrap <- read_delim("prevalence/prevalence_scenario_1c_ci_reps.txt", 
                                      " ", escape_double = FALSE, col_names = FALSE, 
                                      trim_ws = TRUE, skip = 1)[,-1]

covers <- c()
for(r in 1:300){
  lower <- 2*prevalence_1c_ests[r] - quantile(prevalence_1c_bootstrap[r,],
                                              0.975, na.rm=T)
  upper <- 2*prevalence_1c_ests[r] - quantile(prevalence_1c_bootstrap[r,],
                                              0.025, na.rm=T)
  covers[r] <- lower <= 0.4 & upper >= 0.4
  #print(r)
}

prevalence_coverage[2] <- mean(covers, na.rm=T)
prevalence_coverage_sd[2] <- sqrt(mean(covers, 
                                       na.rm=T)*(1-mean(covers, 
                                                        na.rm=T))/sum(!is.na(covers)))
prevalence_means[2] <- mean(prevalence_1c_ests, na.rm=T)
prevalence_sd[2] <- sd(prevalence_1c_ests, na.rm=T)/sqrt(sum(!is.na(prevalence_1c_ests)))



prevalence_2a_ests <- read_delim("prevalence/prevalence_scenario_2a_estimates.txt", 
                                 " ", escape_double = FALSE, col_names = FALSE, 
                                 trim_ws = TRUE, skip = 1)$X2
prevalence_2a_bootstrap <- read_delim("prevalence/prevalence_scenario_2a_ci_reps.txt", 
                                      " ", escape_double = FALSE, col_names = FALSE, 
                                      trim_ws = TRUE, skip = 1)[,-1]

covers <- c()
for(r in 1:300){
  lower <- 2*prevalence_2a_ests[r] - quantile(prevalence_2a_bootstrap[r,],
                                              0.975, na.rm=T)
  upper <- 2*prevalence_2a_ests[r] - quantile(prevalence_2a_bootstrap[r,],
                                              0.025, na.rm=T)
  covers[r] <- lower <= 0.4 & upper >= 0.4
  print(r)
}

prevalence_coverage[3] <- mean(covers, na.rm=T)
prevalence_coverage_sd[3] <- sqrt(mean(covers, 
                                       na.rm=T)*(1-mean(covers, 
                                                        na.rm=T))/sum(!is.na(covers)))
prevalence_means[3] <- mean(prevalence_2a_ests, na.rm=T)
prevalence_sd[3] <- sd(prevalence_2a_ests, na.rm=T)/sqrt(sum(!is.na(prevalence_2a_ests)))





prevalence_2c_ests <- read_delim("prevalence/prevalence_scenario_2c_estimates.txt", 
                                 " ", escape_double = FALSE, col_names = FALSE, 
                                 trim_ws = TRUE, skip = 1)$X2
prevalence_2c_bootstrap <- read_delim("prevalence/prevalence_scenario_2c_ci_reps.txt", 
                                      " ", escape_double = FALSE, col_names = FALSE, 
                                      trim_ws = TRUE, skip = 1)[,-1]

covers <- c()
for(r in 1:300){
  lower <- 2*prevalence_2c_ests[r] - quantile(prevalence_2c_bootstrap[r,],
                                              0.975, na.rm=T)
  upper <- 2*prevalence_2c_ests[r] - quantile(prevalence_2c_bootstrap[r,],
                                              0.025, na.rm=T)
  covers[r] <- lower <= 0.4 & upper >= 0.4
  print(r)
}

prevalence_coverage[4] <- mean(covers, na.rm=T)
prevalence_coverage_sd[4] <- sqrt(mean(covers, 
                                       na.rm=T)*(1-mean(covers, 
                                                        na.rm=T))/sum(!is.na(covers)))
prevalence_means[4] <- mean(prevalence_2c_ests, na.rm=T)
prevalence_sd[4] <- sd(prevalence_2c_ests, na.rm=T)/sqrt(sum(!is.na(prevalence_2c_ests)))



prevalence_3a_ests <- read_delim("prevalence/prevalence_scenario_3a_estimates.txt", 
                                 " ", escape_double = FALSE, col_names = FALSE, 
                                 trim_ws = TRUE, skip = 1)$X2
prevalence_3a_bootstrap <- read_delim("prevalence/prevalence_scenario_3a_ci_reps.txt", 
                                      " ", escape_double = FALSE, col_names = FALSE, 
                                      trim_ws = TRUE, skip = 1)[,-1]

covers <- c()
for(r in 1:300){
  lower <- 2*prevalence_3a_ests[r] - quantile(prevalence_3a_bootstrap[r,],
                                              0.975, na.rm=T)
  upper <- 2*prevalence_3a_ests[r] - quantile(prevalence_3a_bootstrap[r,],
                                              0.025, na.rm=T)
  covers[r] <- lower <= 0.4 & upper >= 0.4
  print(r)
}

prevalence_coverage[5] <- mean(covers, na.rm=T)
prevalence_coverage_sd[5] <- sqrt(mean(covers, 
                                       na.rm=T)*(1-mean(covers, 
                                                        na.rm=T))/sum(!is.na(covers)))
prevalence_means[5] <- mean(prevalence_3a_ests, na.rm=T)
prevalence_sd[5] <- sd(prevalence_3a_ests, na.rm=T)/sqrt(sum(!is.na(prevalence_3a_ests)))


prevalence_3c_ests <- read_delim("prevalence/prevalence_scenario_3c_estimates.txt", 
                                 " ", escape_double = FALSE, col_names = FALSE, 
                                 trim_ws = TRUE, skip = 1)$X2
prevalence_3c_bootstrap <- read_delim("prevalence/prevalence_scenario_3c_ci_reps.txt", 
                                      " ", escape_double = FALSE, col_names = FALSE, 
                                      trim_ws = TRUE, skip = 1)[,-1]

covers <- c()
for(r in 1:300){
  lower <- 2*prevalence_3c_ests[r] - quantile(prevalence_3c_bootstrap[r,],
                                              0.975, na.rm=T)
  upper <- 2*prevalence_3c_ests[r] - quantile(prevalence_3c_bootstrap[r,],
                                              0.025, na.rm=T)
  covers[r] <- lower <= 0.4 & upper >= 0.4
  print(r)
}

prevalence_coverage[6] <- mean(covers, na.rm=T)
prevalence_coverage_sd[6] <- sqrt(mean(covers, 
                                       na.rm=T)*(1-mean(covers, 
                                                        na.rm=T))/sum(!is.na(covers)))
prevalence_means[6] <- mean(prevalence_3c_ests, na.rm=T)
prevalence_sd[6] <- sd(prevalence_3c_ests, na.rm=T)/sqrt(sum(!is.na(prevalence_3c_ests)))



######################################################################
###### Mixed effects model, label-independent random effects   #######
######################################################################

### mixed effects, same puff and nonpuff random effects
### Results are in files 
### ls_bootstrap_scenario_xx_fixed_estimates.txt (point estimates at 
###  each repetition of the simulation) and 
### ls_bootstrap_scenario_xx_fixed_ci_reps.txt (bootstrap replications 
###  at each simulation repetition, to construct CIs)
### These files are created by 
### bootstrap_sim_scenario_xx_fixed.R
me_coverage_fixed <- vector(length=6)
me_coverage_fixed_sd <- vector(length=6)
me_fixed_mean <- vector(length=6)
me_fixed_sd <- vector(length=6)


me_ests_1a_fixed <- read_delim("mixed_effects_model/label_independent_random_effects/ls_bootstrap_scenario_1a_fixed_estimates.txt", 
                               " ", escape_double = FALSE, col_names = FALSE, 
                               trim_ws = TRUE, skip = 1)$X2
bootstrap_ests_1a_fixed <- read_delim("mixed_effects_model/label_independent_random_effects/ls_bootstrap_scenario_1a_fixed_ci_reps.txt", 
                                      " ", escape_double = FALSE, col_names = FALSE, 
                                      trim_ws = TRUE, skip = 1)[,-1]

covers <- c()
for(r in 1:300){
  lower <- 2*me_ests_1a_fixed[r] - quantile(bootstrap_ests_1a_fixed[r,],
                                            0.975, na.rm=T)
  upper <- 2*me_ests_1a_fixed[r] - quantile(bootstrap_ests_1a_fixed[r,],
                                            0.025, na.rm=T)
  covers[r] <- lower <= 3 & upper >= 3
  print(r)
}

me_coverage_fixed[1] <- mean(covers, na.rm=T)
me_coverage_fixed_sd[1] <- sqrt(mean(covers, 
                                     na.rm=T)*(1-mean(covers, 
                                                      na.rm=T))/sum(!is.na(covers)))
me_fixed_mean[1] <- mean(me_ests_1a_fixed, na.rm=T)
me_fixed_sd[1] <- sd(me_ests_1a_fixed, na.rm=T)/sqrt(sum(!is.na(me_ests_1a_fixed)))



me_ests_1c_fixed <- read_delim("mixed_effects_model/label_independent_random_effects/ls_bootstrap_scenario_1c_fixed_estimates.txt", 
                               " ", escape_double = FALSE, col_names = FALSE, 
                               trim_ws = TRUE, skip = 1)$X2
bootstrap_ests_1c_fixed <- read_delim("mixed_effects_model/label_independent_random_effects/ls_bootstrap_scenario_1c_fixed_ci_reps.txt", 
                                      " ", escape_double = FALSE, col_names = FALSE, 
                                      trim_ws = TRUE, skip = 1)[,-1]

covers <- c()
for(r in 1:300){
  lower <- 2*me_ests_1c_fixed[r] - quantile(bootstrap_ests_1c_fixed[r,],
                                            0.975, na.rm=T)
  upper <- 2*me_ests_1c_fixed[r] - quantile(bootstrap_ests_1c_fixed[r,],
                                            0.025, na.rm=T)
  covers[r] <- lower <= 3.485 & upper >= 3.485
  print(r)
}

me_coverage_fixed[2] <- mean(covers, na.rm=T)
me_coverage_fixed_sd[2] <- sqrt(mean(covers, 
                                     na.rm=T)*(1-mean(covers, 
                                                      na.rm=T))/sum(!is.na(covers)))
me_fixed_mean[2] <- mean(me_ests_1c_fixed, na.rm=T)
me_fixed_sd[2] <- sd(me_ests_1c_fixed, na.rm=T)/sqrt(sum(!is.na(me_ests_1c_fixed)))



me_ests_2a_fixed <- read_delim("mixed_effects_model/label_independent_random_effects/ls_bootstrap_scenario_2a_fixed_estimates.txt", 
                               " ", escape_double = FALSE, col_names = FALSE, 
                               trim_ws = TRUE, skip = 1)$X2
bootstrap_ests_2a_fixed <- read_delim("mixed_effects_model/label_independent_random_effects/ls_bootstrap_scenario_2a_fixed_ci_reps.txt", 
                                      " ", escape_double = FALSE, col_names = FALSE, 
                                      trim_ws = TRUE, skip = 1)[,-1]

covers <- c()
for(r in 1:300){
  lower <- 2*me_ests_2a_fixed[r] - quantile(bootstrap_ests_2a_fixed[r,],
                                            0.975, na.rm=T)
  upper <- 2*me_ests_2a_fixed[r] - quantile(bootstrap_ests_2a_fixed[r,],
                                            0.025, na.rm=T)
  covers[r] <- lower <= 3 & upper >= 3
  print(r)
}

me_coverage_fixed[3] <- mean(covers, na.rm=T)
me_coverage_fixed_sd[3] <- sqrt(mean(covers, 
                                     na.rm=T)*(1-mean(covers, 
                                                      na.rm=T))/sum(!is.na(covers)))
me_fixed_mean[3] <- mean(me_ests_2a_fixed, na.rm=T)
me_fixed_sd[3] <- sd(me_ests_2a_fixed, na.rm=T)/sqrt(sum(!is.na(me_ests_2a_fixed)))



me_ests_2c_fixed <- read_delim("mixed_effects_model/label_independent_random_effects/ls_bootstrap_scenario_2c_fixed_estimates.txt", 
                               " ", escape_double = FALSE, col_names = FALSE, 
                               trim_ws = TRUE, skip = 1)$X2
bootstrap_ests_2c_fixed <- read_delim("mixed_effects_model/label_independent_random_effects/ls_bootstrap_scenario_2c_fixed_ci_reps.txt", 
                                      " ", escape_double = FALSE, col_names = FALSE, 
                                      trim_ws = TRUE, skip = 1)[,-1]

covers <- c()
for(r in 1:300){
  lower <- 2*me_ests_2c_fixed[r] - quantile(bootstrap_ests_2c_fixed[r,],
                                            0.975, na.rm=T)
  upper <- 2*me_ests_2c_fixed[r] - quantile(bootstrap_ests_2c_fixed[r,],
                                            0.025, na.rm=T)
  covers[r] <- lower <= 3.485 & upper >= 3.485
  print(r)
}

me_coverage_fixed[4] <- mean(covers, na.rm=T)
me_coverage_fixed_sd[4] <- sqrt(mean(covers, 
                                     na.rm=T)*(1-mean(covers, 
                                                      na.rm=T))/sum(!is.na(covers)))
me_fixed_mean[4] <- mean(me_ests_2c_fixed, na.rm=T)
me_fixed_sd[4] <- sd(me_ests_2c_fixed, na.rm=T)/sqrt(sum(!is.na(me_ests_2c_fixed)))



me_ests_3a_fixed <- read_delim("mixed_effects_model/label_independent_random_effects/ls_bootstrap_scenario_3a_fixed_estimates.txt", 
                               " ", escape_double = FALSE, col_names = FALSE, 
                               trim_ws = TRUE, skip = 1)$X2
bootstrap_ests_3a_fixed <- read_delim("mixed_effects_model/label_independent_random_effects/ls_bootstrap_scenario_3a_fixed_ci_reps.txt", 
                                      " ", escape_double = FALSE, col_names = FALSE, 
                                      trim_ws = TRUE, skip = 1)[,-1]

covers <- c()
for(r in 1:300){
  lower <- 2*me_ests_3a_fixed[r] - quantile(bootstrap_ests_3a_fixed[r,],
                                            0.975, na.rm=T)
  upper <- 2*me_ests_3a_fixed[r] - quantile(bootstrap_ests_3a_fixed[r,],
                                            0.025, na.rm=T)
  covers[r] <- lower <= 4 & upper >= 4
  print(r)
}

me_coverage_fixed[5] <- mean(covers, na.rm=T)
me_coverage_fixed_sd[5] <- sqrt(mean(covers, 
                                     na.rm=T)*(1-mean(covers, 
                                                      na.rm=T))/sum(!is.na(covers)))
me_fixed_mean[5] <- mean(me_ests_3a_fixed, na.rm=T)
me_fixed_sd[5] <- sd(me_ests_3a_fixed, na.rm=T)/sqrt(sum(!is.na(me_ests_3a_fixed)))



me_ests_3c_fixed <- read_delim("mixed_effects_model/label_independent_random_effects/ls_bootstrap_scenario_3c_fixed_estimates.txt", 
                               " ", escape_double = FALSE, col_names = FALSE, 
                               trim_ws = TRUE, skip = 1)$X2
bootstrap_ests_3c_fixed <- read_delim("mixed_effects_model/label_independent_random_effects/ls_bootstrap_scenario_3c_fixed_ci_reps.txt", 
                                      " ", escape_double = FALSE, col_names = FALSE, 
                                      trim_ws = TRUE, skip = 1)[,-1]

covers <- c()
for(r in 1:300){
  lower <- 2*me_ests_3c_fixed[r] - quantile(bootstrap_ests_3c_fixed[r,],
                                            0.975, na.rm=T)
  upper <- 2*me_ests_3c_fixed[r] - quantile(bootstrap_ests_3c_fixed[r,],
                                            0.025, na.rm=T)
  covers[r] <- lower <= 4.485 & upper >= 4.485
  print(r)
}

me_coverage_fixed[6] <- mean(covers, na.rm=T)
me_coverage_fixed_sd[6] <- sqrt(mean(covers, 
                                     na.rm=T)*(1-mean(covers, 
                                                      na.rm=T))/sum(!is.na(covers)))
me_fixed_mean[6] <- mean(me_ests_3c_fixed, na.rm=T)
me_fixed_sd[6] <- sd(me_ests_3c_fixed, na.rm=T)/sqrt(sum(!is.na(me_ests_3c_fixed)))




####################################################################
######## Mixed effects model, label dependent random effects, ######
######## no variance correction                               ######
####################################################################

### mixed effects, different puff and nonpuff random effects
### Results are in files 
### ls_bootstrap_scenario_xx_estimates.txt (point estimates at 
###  each repetition of the simulation) and 
### ls_bootstrap_scenario_xx_ci_reps.txt (bootstrap replications 
###  at each simulation repetition, to construct CIs)
### These files are created by 
### bootstrap_sim_scenario_xx.R

me_coverage <- vector(length=6)
me_coverage_sd <- vector(length=6)
me_mean <- vector(length=6)
me_sd <- vector(length=6)


me_ests_1a <- read_delim("mixed_effects_model/label_dependent_random_effects/no_variance_correction/ls_bootstrap_scenario_1a_estimates.txt", 
                         " ", escape_double = FALSE, col_names = FALSE, 
                         trim_ws = TRUE, skip = 1)$X2
bootstrap_ests_1a <- read_delim("mixed_effects_model/label_dependent_random_effects/no_variance_correction/ls_bootstrap_scenario_1a_ci_reps.txt", 
                                " ", escape_double = FALSE, col_names = FALSE, 
                                trim_ws = TRUE, skip = 1)[,-1]

covers <- c()
for(r in 1:300){
  lower <- 2*me_ests_1a[r] - quantile(bootstrap_ests_1a[r,],
                                      0.975, na.rm=T)
  upper <- 2*me_ests_1a[r] - quantile(bootstrap_ests_1a[r,],
                                      0.025, na.rm=T)
  covers[r] <- lower <= 3 & upper >= 3
  print(r)
}

me_coverage[1] <- mean(covers, na.rm=T)
me_coverage_sd[1] <- sqrt(mean(covers, 
                               na.rm=T)*(1-mean(covers, 
                                                na.rm=T))/sum(!is.na(covers)))
me_mean[1] <- mean(me_ests_1a, na.rm=T)
me_sd[1] <- sd(me_ests_1a, na.rm=T)/sqrt(sum(!is.na(me_ests_1a)))



me_ests_1c <- read_delim("mixed_effects_model/label_dependent_random_effects/no_variance_correction/ls_bootstrap_scenario_1c_estimates.txt", 
                         " ", escape_double = FALSE, col_names = FALSE, 
                         trim_ws = TRUE, skip = 1)$X2
bootstrap_ests_1c <- read_delim("mixed_effects_model/label_dependent_random_effects/no_variance_correction/ls_bootstrap_scenario_1c_ci_reps.txt", 
                                " ", escape_double = FALSE, col_names = FALSE, 
                                trim_ws = TRUE, skip = 1)[,-1]

covers <- c()
for(r in 1:300){
  lower <- 2*me_ests_1c[r] - quantile(bootstrap_ests_1c[r,],
                                      0.975, na.rm=T)
  upper <- 2*me_ests_1c[r] - quantile(bootstrap_ests_1c[r,],
                                      0.025, na.rm=T)
  covers[r] <- lower <= 3.485 & upper >= 3.485
  print(r)
}

me_coverage[2] <- mean(covers, na.rm=T)
me_coverage_sd[2] <- sqrt(mean(covers, 
                               na.rm=T)*(1-mean(covers, 
                                                na.rm=T))/sum(!is.na(covers)))
me_mean[2] <- mean(me_ests_1c, na.rm=T)
me_sd[2] <- sd(me_ests_1c, na.rm=T)/sqrt(sum(!is.na(me_ests_1c)))



me_ests_2a <- read_delim("mixed_effects_model/label_dependent_random_effects/no_variance_correction/ls_bootstrap_scenario_2a_estimates.txt", 
                         " ", escape_double = FALSE, col_names = FALSE, 
                         trim_ws = TRUE, skip = 1)$X2
bootstrap_ests_2a <- read_delim("mixed_effects_model/label_dependent_random_effects/no_variance_correction/ls_bootstrap_scenario_2a_ci_reps.txt", 
                                " ", escape_double = FALSE, col_names = FALSE, 
                                trim_ws = TRUE, skip = 1)[,-1]

covers <- c()
for(r in 1:300){
  lower <- 2*me_ests_2a[r] - quantile(bootstrap_ests_2a[r,],
                                      0.975, na.rm=T)
  upper <- 2*me_ests_2a[r] - quantile(bootstrap_ests_2a[r,],
                                      0.025, na.rm=T)
  covers[r] <- lower <= 3 & upper >= 3
  print(r)
}

me_coverage[3] <- mean(covers, na.rm=T)
me_coverage_sd[3] <- sqrt(mean(covers, 
                               na.rm=T)*(1-mean(covers, 
                                                na.rm=T))/sum(!is.na(covers)))
me_mean[3] <- mean(me_ests_2a, na.rm=T)
me_sd[3] <- sd(me_ests_2a, na.rm=T)/sqrt(sum(!is.na(me_ests_2a)))




me_ests_2c <- read_delim("mixed_effects_model/label_dependent_random_effects/no_variance_correction/ls_bootstrap_scenario_2c_estimates.txt", 
                         " ", escape_double = FALSE, col_names = FALSE, 
                         trim_ws = TRUE, skip = 1)$X2
bootstrap_ests_2c <- read_delim("mixed_effects_model/label_dependent_random_effects/no_variance_correction/ls_bootstrap_scenario_2c_ci_reps.txt", 
                                " ", escape_double = FALSE, col_names = FALSE, 
                                trim_ws = TRUE, skip = 1)[,-1]

covers <- c()
for(r in 1:300){
  lower <- 2*me_ests_2c[r] - quantile(bootstrap_ests_2c[r,],
                                      0.975, na.rm=T)
  upper <- 2*me_ests_2c[r] - quantile(bootstrap_ests_2c[r,],
                                      0.025, na.rm=T)
  covers[r] <- lower <= 3.485 & upper >= 3.485
  print(r)
}

me_coverage[4] <- mean(covers, na.rm=T)
me_coverage_sd[4] <- sqrt(mean(covers, 
                               na.rm=T)*(1-mean(covers, 
                                                na.rm=T))/sum(!is.na(covers)))
me_mean[4] <- mean(me_ests_2c, na.rm=T)
me_sd[4] <- sd(me_ests_2c, na.rm=T)/sqrt(sum(!is.na(me_ests_2c)))




me_ests_3a <- read_delim("mixed_effects_model/label_dependent_random_effects/no_variance_correction/ls_bootstrap_scenario_3a_estimates.txt", 
                         " ", escape_double = FALSE, col_names = FALSE, 
                         trim_ws = TRUE, skip = 1)$X2
bootstrap_ests_3a <- read_delim("mixed_effects_model/label_dependent_random_effects/no_variance_correction/ls_bootstrap_scenario_3a_ci_reps.txt", 
                                " ", escape_double = FALSE, col_names = FALSE, 
                                trim_ws = TRUE, skip = 1)[,-1]

covers <- c()
for(r in 1:300){
  lower <- 2*me_ests_3a[r] - quantile(bootstrap_ests_3a[r,],
                                      0.975, na.rm=T)
  upper <- 2*me_ests_3a[r] - quantile(bootstrap_ests_3a[r,],
                                      0.025, na.rm=T)
  covers[r] <- lower <= 4 & upper >= 4
  print(r)
}

me_coverage[5] <- mean(covers, na.rm=T)
me_coverage_sd[5] <- sqrt(mean(covers, 
                               na.rm=T)*(1-mean(covers, 
                                                na.rm=T))/sum(!is.na(covers)))
me_mean[5] <- mean(me_ests_3a, na.rm=T)
me_sd[5] <- sd(me_ests_3a, na.rm=T)/sqrt(sum(!is.na(me_ests_3a)))




me_ests_3c <- read_delim("mixed_effects_model/label_dependent_random_effects/no_variance_correction/ls_bootstrap_scenario_3c_estimates.txt", 
                         " ", escape_double = FALSE, col_names = FALSE, 
                         trim_ws = TRUE, skip = 1)$X2
bootstrap_ests_3c <- read_delim("mixed_effects_model/label_dependent_random_effects/no_variance_correction/ls_bootstrap_scenario_3c_ci_reps.txt", 
                                " ", escape_double = FALSE, col_names = FALSE, 
                                trim_ws = TRUE, skip = 1)[,-1]

covers <- c()
for(r in 1:300){
  lower <- 2*me_ests_3c[r] - quantile(bootstrap_ests_3c[r,],
                                      0.975, na.rm=T)
  upper <- 2*me_ests_3c[r] - quantile(bootstrap_ests_3c[r,],
                                      0.025, na.rm=T)
  covers[r] <- lower <= 4.485 & upper >= 4.485
  print(r)
}

me_coverage[6] <- mean(covers, na.rm=T)
me_coverage_sd[6] <- sqrt(mean(covers, 
                               na.rm=T)*(1-mean(covers, 
                                                na.rm=T))/sum(!is.na(covers)))
me_mean[6] <- mean(me_ests_3c, na.rm=T)
me_sd[6] <- sd(me_ests_3c, na.rm=T)/sqrt(sum(!is.na(me_ests_3c)))



####################################################################
######## Mixed effects model, label dependent random effects, ######
######## variance correction                               ######
####################################################################

### mixed effects, different puff and nonpuff random effects
### Variance correction to adjust for label-dependence
### Results are in files 
### ls_bootstrap_scenario_xx_adj_estimates.txt (point estimates at 
###  each repetition of the simulation) and 
### ls_bootstrap_scenario_xx_adj_ci_reps.txt (bootstrap replications 
###  at each simulation repetition, to construct CIs)
### These files are created by 
### bootstrap_sim_scenario_xx_adj.R

me_coverage_adj <- vector(length=6)
me_coverage_adj_sd <- vector(length=6)
me_adj_mean <- vector(length=6)
me_adj_sd <- vector(length=6)


me_ests_1a_adj <- read_delim("mixed_effects_model/label_dependent_random_effects/variance_correction/ls_bootstrap_scenario_1a_adj_estimates.txt", 
                             " ", escape_double = FALSE, col_names = FALSE, 
                             trim_ws = TRUE, skip = 1)$X2
bootstrap_ests_1a_adj <- read_delim("mixed_effects_model/label_dependent_random_effects/variance_correction/ls_bootstrap_scenario_1a_adj_ci_reps.txt", 
                                    " ", escape_double = FALSE, col_names = FALSE, 
                                    trim_ws = TRUE, skip = 1)[,-1]

covers <- c()
for(r in 1:300){
  lower <- 2*me_ests_1a_adj[r] - quantile(bootstrap_ests_1a_adj[r,],
                                          0.975, na.rm=T)
  upper <- 2*me_ests_1a_adj[r] - quantile(bootstrap_ests_1a_adj[r,],
                                          0.025, na.rm=T)
  covers[r] <- lower <= 3 & upper >= 3
  print(r)
}

me_coverage_adj[1] <- mean(covers, na.rm=T)
me_coverage_adj_sd[1] <- sqrt(mean(covers, 
                                   na.rm=T)*(1-mean(covers, 
                                                    na.rm=T))/sum(!is.na(covers)))
me_adj_mean[1] <- mean(me_ests_1a_adj, na.rm=T)
me_adj_sd[1] <- sd(me_ests_1a_adj, na.rm=T)/sqrt(sum(!is.na(me_ests_1a_adj)))



me_ests_1c_adj <- read_delim("mixed_effects_model/label_dependent_random_effects/variance_correction/ls_bootstrap_scenario_1c_adj_estimates.txt", 
                             " ", escape_double = FALSE, col_names = FALSE, 
                             trim_ws = TRUE, skip = 1)$X2
bootstrap_ests_1c_adj <- read_delim("mixed_effects_model/label_dependent_random_effects/variance_correction/ls_bootstrap_scenario_1c_adj_ci_reps.txt", 
                                    " ", escape_double = FALSE, col_names = FALSE, 
                                    trim_ws = TRUE, skip = 1)[,-1]

covers <- c()
for(r in 1:300){
  lower <- 2*me_ests_1c_adj[r] - quantile(bootstrap_ests_1c_adj[r,],
                                          0.975, na.rm=T)
  upper <- 2*me_ests_1c_adj[r] - quantile(bootstrap_ests_1c_adj[r,],
                                          0.025, na.rm=T)
  covers[r] <- lower <= 3.485 & upper >= 3.485
  print(r)
}

me_coverage_adj[2] <- mean(covers, na.rm=T)
me_coverage_adj_sd[2] <- sqrt(mean(covers, 
                                   na.rm=T)*(1-mean(covers, 
                                                    na.rm=T))/sum(!is.na(covers)))
me_adj_mean[2] <- mean(me_ests_1c_adj, na.rm=T)
me_adj_sd[2] <- sd(me_ests_1c_adj, na.rm=T)/sqrt(sum(!is.na(me_ests_1c_adj)))



me_ests_2a_adj <- read_delim("mixed_effects_model/label_dependent_random_effects/variance_correction/ls_bootstrap_scenario_2a_adj_estimates.txt", 
                             " ", escape_double = FALSE, col_names = FALSE, 
                             trim_ws = TRUE, skip = 1)$X2
bootstrap_ests_2a_adj <- read_delim("mixed_effects_model/label_dependent_random_effects/variance_correction/ls_bootstrap_scenario_2a_adj_ci_reps.txt", 
                                    " ", escape_double = FALSE, col_names = FALSE, 
                                    trim_ws = TRUE, skip = 1)[,-1]

covers <- c()
for(r in 1:300){
  lower <- 2*me_ests_2a_adj[r] - quantile(bootstrap_ests_2a_adj[r,],
                                          0.975, na.rm=T)
  upper <- 2*me_ests_2a_adj[r] - quantile(bootstrap_ests_2a_adj[r,],
                                          0.025, na.rm=T)
  covers[r] <- lower <= 3 & upper >= 3
  print(r)
}

me_coverage_adj[3] <- mean(covers, na.rm=T)
me_coverage_adj_sd[3] <- sqrt(mean(covers, 
                                   na.rm=T)*(1-mean(covers, 
                                                    na.rm=T))/sum(!is.na(covers)))
me_adj_mean[3] <- mean(me_ests_2a_adj, na.rm=T)
me_adj_sd[3] <- sd(me_ests_2a_adj, na.rm=T)/sqrt(sum(!is.na(me_ests_2a_adj)))



me_ests_2c_adj <- read_delim("mixed_effects_model/label_dependent_random_effects/variance_correction/ls_bootstrap_scenario_2c_adj_estimates.txt", 
                             " ", escape_double = FALSE, col_names = FALSE, 
                             trim_ws = TRUE, skip = 1)$X2
bootstrap_ests_2c_adj <- read_delim("mixed_effects_model/label_dependent_random_effects/variance_correction/ls_bootstrap_scenario_2c_adj_ci_reps.txt", 
                                    " ", escape_double = FALSE, col_names = FALSE, 
                                    trim_ws = TRUE, skip = 1)[,-1]

covers <- c()
for(r in 1:300){
  lower <- 2*me_ests_2c_adj[r] - quantile(bootstrap_ests_2c_adj[r,],
                                          0.975, na.rm=T)
  upper <- 2*me_ests_2c_adj[r] - quantile(bootstrap_ests_2c_adj[r,],
                                          0.025, na.rm=T)
  covers[r] <- lower <= 3.485 & upper >= 3.485
  print(r)
}

me_coverage_adj[4] <- mean(covers, na.rm=T)
me_coverage_adj_sd[4] <- sqrt(mean(covers, 
                                   na.rm=T)*(1-mean(covers, 
                                                    na.rm=T))/sum(!is.na(covers)))
me_adj_mean[4] <- mean(me_ests_2c_adj, na.rm=T)
me_adj_sd[4] <- sd(me_ests_2c_adj, na.rm=T)/sqrt(sum(!is.na(me_ests_2c_adj)))



me_ests_3a_adj <- read_delim("mixed_effects_model/label_dependent_random_effects/variance_correction/ls_bootstrap_scenario_3a_adj_estimates.txt", 
                             " ", escape_double = FALSE, col_names = FALSE, 
                             trim_ws = TRUE, skip = 1)$X2
bootstrap_ests_3a_adj <- read_delim("mixed_effects_model/label_dependent_random_effects/variance_correction/ls_bootstrap_scenario_3a_adj_ci_reps.txt", 
                                    " ", escape_double = FALSE, col_names = FALSE, 
                                    trim_ws = TRUE, skip = 1)[,-1]

covers <- c()
for(r in 1:300){
  lower <- 2*me_ests_3a_adj[r] - quantile(bootstrap_ests_3a_adj[r,],
                                          0.975, na.rm=T)
  upper <- 2*me_ests_3a_adj[r] - quantile(bootstrap_ests_3a_adj[r,],
                                          0.025, na.rm=T)
  covers[r] <- lower <= 4 & upper >= 4
  print(r)
}

me_coverage_adj[5] <- mean(covers, na.rm=T)
me_coverage_adj_sd[5] <- sqrt(mean(covers, 
                                   na.rm=T)*(1-mean(covers, 
                                                    na.rm=T))/sum(!is.na(covers)))
me_adj_mean[5] <- mean(me_ests_3a_adj, na.rm=T)
me_adj_sd[5] <- sd(me_ests_3a_adj, na.rm=T)/sqrt(sum(!is.na(me_ests_3a_adj)))



me_ests_3c_adj <- read_delim("mixed_effects_model/label_dependent_random_effects/variance_correction/ls_bootstrap_scenario_3c_adj_estimates.txt", 
                             " ", escape_double = FALSE, col_names = FALSE, 
                             trim_ws = TRUE, skip = 1)$X2
bootstrap_ests_3c_adj <- read_delim("mixed_effects_model/label_dependent_random_effects/variance_correction/ls_bootstrap_scenario_3c_adj_ci_reps.txt", 
                                    " ", escape_double = FALSE, col_names = FALSE, 
                                    trim_ws = TRUE, skip = 1)[,-1]

covers <- c()
for(r in 1:300){
  lower <- 2*me_ests_3c_adj[r] - quantile(bootstrap_ests_3c_adj[r,],
                                          0.975, na.rm=T)
  upper <- 2*me_ests_3c_adj[r] - quantile(bootstrap_ests_3c_adj[r,],
                                          0.025, na.rm=T)
  covers[r] <- lower <= 4.485 & upper >= 4.485
  print(r)
}

me_coverage_adj[6] <- mean(covers, na.rm=T)
me_coverage_adj_sd[6] <- sqrt(mean(covers, 
                                   na.rm=T)*(1-mean(covers, 
                                                    na.rm=T))/sum(!is.na(covers)))
me_adj_mean[6] <- mean(me_ests_3c_adj, na.rm=T)
me_adj_sd[6] <- sd(me_ests_3c_adj, na.rm=T)/sqrt(sum(!is.na(me_ests_3c_adj)))




#############################################################
###### Mixture model, label-independent random effects ######
#############################################################

### mixture model, same puff and nonpuff random effects
### Results are in files 
### ls_mixture_scenario_xx_fixed_estimates.txt (point estimates at 
###  each repetition of the simulation) and 
### ls_mixture_scenario_xx_fixed_ci_reps.txt (bootstrap replications 
###  at each simulation repetition, to construct CIs)
### These files are created by 
### mixture_sim_scenario_xx_fixed.R

mixture_coverage_fixed <- vector(length=6)
mixture_coverage_fixed_sd <- vector(length=6)
mixture_fixed_mean <- vector(length=6)
mixture_fixed_sd <- vector(length=6)


mixture_ests_1a_fixed <- read_delim("mixture_model/label_independent_random_effects/ls_mixture_scenario_1a_fixed_estimates.txt", 
                                    " ", escape_double = FALSE, col_names = FALSE, 
                                    trim_ws = TRUE, skip = 1)$X2
bootstrap_ests_1a_fixed <- read_delim("mixture_model/label_independent_random_effects/ls_mixture_scenario_1a_fixed_ci_reps.txt", 
                                      " ", escape_double = FALSE, col_names = FALSE, 
                                      trim_ws = TRUE, skip = 1)[,-1]

covers <- c()
for(r in 1:300){
  lower <- 2*mixture_ests_1a_fixed[r] - quantile(bootstrap_ests_1a_fixed[r,],
                                                 0.975, na.rm=T)
  upper <- 2*mixture_ests_1a_fixed[r] - quantile(bootstrap_ests_1a_fixed[r,],
                                                 0.025, na.rm=T)
  covers[r] <- lower <= 3 & upper >= 3
  print(r)
}

mixture_coverage_fixed[1] <- mean(covers, na.rm=T)
mixture_coverage_fixed_sd[1] <- sqrt(mean(covers, 
                                          na.rm=T)*(1-mean(covers, 
                                                           na.rm=T))/sum(!is.na(covers)))
mixture_fixed_mean[1] <- mean(mixture_ests_1a_fixed, na.rm=T)
mixture_fixed_sd[1] <- sd(mixture_ests_1a_fixed, na.rm=T)/sqrt(sum(!is.na(mixture_ests_1a_fixed)))



mixture_ests_1c_fixed <- read_delim("mixture_model/label_independent_random_effects/ls_mixture_scenario_1c_fixed_estimates.txt", 
                                    " ", escape_double = FALSE, col_names = FALSE, 
                                    trim_ws = TRUE, skip = 1)$X2
bootstrap_ests_1c_fixed <- read_delim("mixture_model/label_independent_random_effects/ls_mixture_scenario_1c_fixed_ci_reps.txt", 
                                      " ", escape_double = FALSE, col_names = FALSE, 
                                      trim_ws = TRUE, skip = 1)[,-1]

covers <- c()
for(r in 1:300){
  lower <- 2*mixture_ests_1c_fixed[r] - quantile(bootstrap_ests_1c_fixed[r,],
                                                 0.975, na.rm=T)
  upper <- 2*mixture_ests_1c_fixed[r] - quantile(bootstrap_ests_1c_fixed[r,],
                                                 0.025, na.rm=T)
  covers[r] <- lower <= 3.485 & upper >= 3.485
  print(r)
}

mixture_coverage_fixed[2] <- mean(covers, na.rm=T)
mixture_coverage_fixed_sd[2] <- sqrt(mean(covers, 
                                          na.rm=T)*(1-mean(covers, 
                                                           na.rm=T))/sum(!is.na(covers)))
mixture_fixed_mean[2] <- mean(mixture_ests_1c_fixed, na.rm=T)
mixture_fixed_sd[2] <- sd(mixture_ests_1c_fixed, na.rm=T)/sqrt(sum(!is.na(mixture_ests_1c_fixed)))



mixture_ests_2a_fixed <- read_delim("mixture_model/label_independent_random_effects/ls_mixture_scenario_2a_fixed_estimates.txt", 
                                    " ", escape_double = FALSE, col_names = FALSE, 
                                    trim_ws = TRUE, skip = 1)$X2
bootstrap_ests_2a_fixed <- read_delim("mixture_model/label_independent_random_effects/ls_mixture_scenario_2a_fixed_ci_reps.txt", 
                                      " ", escape_double = FALSE, col_names = FALSE, 
                                      trim_ws = TRUE, skip = 1)[,-1]

covers <- c()
for(r in 1:300){
  lower <- 2*mixture_ests_2a_fixed[r] - quantile(bootstrap_ests_2a_fixed[r,],
                                                 0.975, na.rm=T)
  upper <- 2*mixture_ests_2a_fixed[r] - quantile(bootstrap_ests_2a_fixed[r,],
                                                 0.025, na.rm=T)
  covers[r] <- lower <= 3 & upper >= 3
  print(r)
}

mixture_coverage_fixed[3] <- mean(covers, na.rm=T)
mixture_coverage_fixed_sd[3] <- sqrt(mean(covers, 
                                          na.rm=T)*(1-mean(covers, 
                                                           na.rm=T))/sum(!is.na(covers)))
mixture_fixed_mean[3] <- mean(mixture_ests_2a_fixed, na.rm=T)
mixture_fixed_sd[3] <- sd(mixture_ests_2a_fixed, na.rm=T)/sqrt(sum(!is.na(mixture_ests_2a_fixed)))



mixture_ests_2c_fixed <- read_delim("mixture_model/label_independent_random_effects/ls_mixture_scenario_2c_fixed_estimates.txt", 
                                    " ", escape_double = FALSE, col_names = FALSE, 
                                    trim_ws = TRUE, skip = 1)$X2
bootstrap_ests_2c_fixed <- read_delim("mixture_model/label_independent_random_effects/ls_mixture_scenario_2c_fixed_ci_reps.txt", 
                                      " ", escape_double = FALSE, col_names = FALSE, 
                                      trim_ws = TRUE, skip = 1)[,-1]

covers <- c()
for(r in 1:300){
  lower <- 2*mixture_ests_2c_fixed[r] - quantile(bootstrap_ests_2c_fixed[r,],
                                                 0.975, na.rm=T)
  upper <- 2*mixture_ests_2c_fixed[r] - quantile(bootstrap_ests_2c_fixed[r,],
                                                 0.025, na.rm=T)
  covers[r] <- lower <= 3.485 & upper >= 3.485
  print(r)
}

mixture_coverage_fixed[4] <- mean(covers, na.rm=T)
mixture_coverage_fixed_sd[4] <- sqrt(mean(covers, 
                                          na.rm=T)*(1-mean(covers, 
                                                           na.rm=T))/sum(!is.na(covers)))
mixture_fixed_mean[4] <- mean(mixture_ests_2c_fixed, na.rm=T)
mixture_fixed_sd[4] <- sd(mixture_ests_2c_fixed, na.rm=T)/sqrt(sum(!is.na(mixture_ests_2c_fixed)))



mixture_ests_3a_fixed <- read_delim("mixture_model/label_independent_random_effects/ls_mixture_scenario_3a_fixed_estimates.txt", 
                                    " ", escape_double = FALSE, col_names = FALSE, 
                                    trim_ws = TRUE, skip = 1)$X2
bootstrap_ests_3a_fixed <- read_delim("mixture_model/label_independent_random_effects/ls_mixture_scenario_3a_fixed_ci_reps.txt", 
                                      " ", escape_double = FALSE, col_names = FALSE, 
                                      trim_ws = TRUE, skip = 1)[,-1]

covers <- c()
for(r in 1:300){
  lower <- 2*mixture_ests_3a_fixed[r] - quantile(bootstrap_ests_3a_fixed[r,],
                                                 0.975, na.rm=T)
  upper <- 2*mixture_ests_3a_fixed[r] - quantile(bootstrap_ests_3a_fixed[r,],
                                                 0.025, na.rm=T)
  covers[r] <- lower <= 4 & upper >= 4
  print(r)
}

mixture_coverage_fixed[5] <- mean(covers, na.rm=T)
mixture_coverage_fixed_sd[5] <- sqrt(mean(covers, 
                                          na.rm=T)*(1-mean(covers, 
                                                           na.rm=T))/sum(!is.na(covers)))
mixture_fixed_mean[5] <- mean(mixture_ests_3a_fixed, na.rm=T)
mixture_fixed_sd[5] <- sd(mixture_ests_3a_fixed, na.rm=T)/sqrt(sum(!is.na(mixture_ests_3a_fixed)))



mixture_ests_3c_fixed <- read_delim("mixture_model/label_independent_random_effects/ls_mixture_scenario_3c_fixed_estimates.txt", 
                                    " ", escape_double = FALSE, col_names = FALSE, 
                                    trim_ws = TRUE, skip = 1)$X2
bootstrap_ests_3c_fixed <- read_delim("mixture_model/label_independent_random_effects/ls_mixture_scenario_3c_fixed_ci_reps.txt", 
                                      " ", escape_double = FALSE, col_names = FALSE, 
                                      trim_ws = TRUE, skip = 1)[,-1]

covers <- c()
for(r in 1:300){
  lower <- 2*mixture_ests_3c_fixed[r] - quantile(bootstrap_ests_3c_fixed[r,],
                                                 0.975, na.rm=T)
  upper <- 2*mixture_ests_3c_fixed[r] - quantile(bootstrap_ests_3c_fixed[r,],
                                                 0.025, na.rm=T)
  covers[r] <- lower <= 4.485 & upper >= 4.485
  print(r)
}

mixture_coverage_fixed[6] <- mean(covers, na.rm=T)
mixture_coverage_fixed_sd[6] <- sqrt(mean(covers, 
                                          na.rm=T)*(1-mean(covers, 
                                                           na.rm=T))/sum(!is.na(covers)))
mixture_fixed_mean[6] <- mean(mixture_ests_3c_fixed, na.rm=T)
mixture_fixed_sd[6] <- sd(mixture_ests_3c_fixed, na.rm=T)/sqrt(sum(!is.na(mixture_ests_3c_fixed)))



#############################################################
###### Mixture model, label-dependent random effects ######
#############################################################

### mixture model, different puff and nonpuff random effects
### Results are in files 
### ls_mixture_scenario_xx_estimates.txt (point estimates at 
###  each repetition of the simulation) and 
### ls_mixture_scenario_xx_ci_reps.txt (bootstrap replications 
###  at each simulation repetition, to construct CIs)
### These files are created by 
### mixture_sim_scenario_xx.R

mixture_coverage_dependent <- vector(length=6)
mixture_coverage_dependent_sd <- vector(length=6)
mixture_dependent_mean <- vector(length=6)
mixture_dependent_sd <- vector(length=6)


mixture_ests_1a_dependent <- read_delim("mixture_model/label_dependent_random_effects/ls_mixture_scenario_1a_estimates.txt", 
                                        " ", escape_double = FALSE, col_names = FALSE, 
                                        trim_ws = TRUE, skip = 1)$X2
bootstrap_ests_1a_dependent <- read_delim("mixture_model/label_dependent_random_effects/ls_mixture_scenario_1a_ci_reps.txt", 
                                          " ", escape_double = FALSE, col_names = FALSE, 
                                          trim_ws = TRUE, skip = 1)[,-1]

covers <- c()
for(r in 1:300){
  lower <- 2*mixture_ests_1a_dependent[r] - quantile(bootstrap_ests_1a_dependent[r,],
                                                     0.975, na.rm=T)
  upper <- 2*mixture_ests_1a_dependent[r] - quantile(bootstrap_ests_1a_dependent[r,],
                                                     0.025, na.rm=T)
  covers[r] <- lower <= 3 & upper >= 3
  print(r)
}

mixture_coverage_dependent[1] <- mean(covers, na.rm=T)
mixture_coverage_dependent_sd[1] <- sqrt(mean(covers, 
                                              na.rm=T)*(1-mean(covers, 
                                                               na.rm=T))/sum(!is.na(covers)))
mixture_dependent_mean[1] <- mean(mixture_ests_1a_dependent, na.rm=T)
mixture_dependent_sd[1] <- sd(mixture_ests_1a_dependent, na.rm=T)/sqrt(sum(!is.na(mixture_ests_1a_dependent)))



mixture_ests_1c_dependent <- read_delim("mixture_model/label_dependent_random_effects/ls_mixture_scenario_1c_estimates.txt", 
                                        " ", escape_double = FALSE, col_names = FALSE, 
                                        trim_ws = TRUE, skip = 1)$X2
bootstrap_ests_1c_dependent <- read_delim("mixture_model/label_dependent_random_effects/ls_mixture_scenario_1c_ci_reps.txt", 
                                          " ", escape_double = FALSE, col_names = FALSE, 
                                          trim_ws = TRUE, skip = 1)[,-1]

covers <- c()
for(r in 1:300){
  lower <- 2*mixture_ests_1c_dependent[r] - quantile(bootstrap_ests_1c_dependent[r,],
                                                     0.975, na.rm=T)
  upper <- 2*mixture_ests_1c_dependent[r] - quantile(bootstrap_ests_1c_dependent[r,],
                                                     0.025, na.rm=T)
  covers[r] <- lower <= 3.485 & upper >= 3.485
  print(r)
}

mixture_coverage_dependent[2] <- mean(covers, na.rm=T)
mixture_coverage_dependent_sd[2] <- sqrt(mean(covers, 
                                              na.rm=T)*(1-mean(covers, 
                                                               na.rm=T))/sum(!is.na(covers)))
mixture_dependent_mean[2] <- mean(mixture_ests_1c_dependent, na.rm=T)
mixture_dependent_sd[2] <- sd(mixture_ests_1c_dependent, na.rm=T)/sqrt(sum(!is.na(mixture_ests_1c_dependent)))



mixture_ests_2a_dependent <- read_delim("mixture_model/label_dependent_random_effects/ls_mixture_scenario_2a_estimates.txt", 
                                        " ", escape_double = FALSE, col_names = FALSE, 
                                        trim_ws = TRUE, skip = 1)$X2
bootstrap_ests_2a_dependent <- read_delim("mixture_model/label_dependent_random_effects/ls_mixture_scenario_2a_ci_reps.txt", 
                                          " ", escape_double = FALSE, col_names = FALSE, 
                                          trim_ws = TRUE, skip = 1)[,-1]

covers <- c()
for(r in 1:300){
  lower <- 2*mixture_ests_2a_dependent[r] - quantile(bootstrap_ests_2a_dependent[r,],
                                                     0.975, na.rm=T)
  upper <- 2*mixture_ests_2a_dependent[r] - quantile(bootstrap_ests_2a_dependent[r,],
                                                     0.025, na.rm=T)
  covers[r] <- lower <= 3 & upper >= 3
  print(r)
}

mixture_coverage_dependent[3] <- mean(covers, na.rm=T)
mixture_coverage_dependent_sd[3] <- sqrt(mean(covers, 
                                              na.rm=T)*(1-mean(covers, 
                                                               na.rm=T))/sum(!is.na(covers)))
mixture_dependent_mean[3] <- mean(mixture_ests_2a_dependent, na.rm=T)
mixture_dependent_sd[3] <- sd(mixture_ests_2a_dependent, na.rm=T)/sqrt(sum(!is.na(mixture_ests_2a_dependent)))



mixture_ests_2c_dependent <- read_delim("mixture_model/label_dependent_random_effects/ls_mixture_scenario_2c_estimates.txt", 
                                        " ", escape_double = FALSE, col_names = FALSE, 
                                        trim_ws = TRUE, skip = 1)$X2
bootstrap_ests_2c_dependent <- read_delim("mixture_model/label_dependent_random_effects/ls_mixture_scenario_2c_ci_reps.txt", 
                                          " ", escape_double = FALSE, col_names = FALSE, 
                                          trim_ws = TRUE, skip = 1)[,-1]

covers <- c()
for(r in 1:300){
  lower <- 2*mixture_ests_2c_dependent[r] - quantile(bootstrap_ests_2c_dependent[r,],
                                                     0.975, na.rm=T)
  upper <- 2*mixture_ests_2c_dependent[r] - quantile(bootstrap_ests_2c_dependent[r,],
                                                     0.025, na.rm=T)
  covers[r] <- lower <= 3.485 & upper >= 3.485
  print(r)
}

mixture_coverage_dependent[4] <- mean(covers, na.rm=T)
mixture_coverage_dependent_sd[4] <- sqrt(mean(covers, 
                                              na.rm=T)*(1-mean(covers, 
                                                               na.rm=T))/sum(!is.na(covers)))
mixture_dependent_mean[4] <- mean(mixture_ests_2c_dependent, na.rm=T)
mixture_dependent_sd[4] <- sd(mixture_ests_2c_dependent, na.rm=T)/sqrt(sum(!is.na(mixture_ests_2c_dependent)))



mixture_ests_3a_dependent <- read_delim("mixture_model/label_dependent_random_effects/ls_mixture_scenario_3a_estimates.txt", 
                                        " ", escape_double = FALSE, col_names = FALSE, 
                                        trim_ws = TRUE, skip = 1)$X2
bootstrap_ests_3a_dependent <- read_delim("mixture_model/label_dependent_random_effects/ls_mixture_scenario_3a_ci_reps.txt", 
                                          " ", escape_double = FALSE, col_names = FALSE, 
                                          trim_ws = TRUE, skip = 1)[,-1]

covers <- c()
for(r in 1:300){
  lower <- 2*mixture_ests_3a_dependent[r] - quantile(bootstrap_ests_3a_dependent[r,],
                                                     0.975, na.rm=T)
  upper <- 2*mixture_ests_3a_dependent[r] - quantile(bootstrap_ests_3a_dependent[r,],
                                                     0.025, na.rm=T)
  covers[r] <- lower <= 4 & upper >= 4
  print(r)
}

mixture_coverage_dependent[5] <- mean(covers, na.rm=T)
mixture_coverage_dependent_sd[5] <- sqrt(mean(covers, 
                                              na.rm=T)*(1-mean(covers, 
                                                               na.rm=T))/sum(!is.na(covers)))
mixture_dependent_mean[5] <- mean(mixture_ests_3a_dependent, na.rm=T)
mixture_dependent_sd[5] <- sd(mixture_ests_3a_dependent, na.rm=T)/sqrt(sum(!is.na(mixture_ests_3a_dependent)))



mixture_ests_3c_dependent <- read_delim("mixture_model/label_dependent_random_effects/ls_mixture_scenario_3c_estimates.txt", 
                                        " ", escape_double = FALSE, col_names = FALSE, 
                                        trim_ws = TRUE, skip = 1)$X2
bootstrap_ests_3c_dependent <- read_delim("mixture_model/label_dependent_random_effects/ls_mixture_scenario_3c_ci_reps.txt", 
                                          " ", escape_double = FALSE, col_names = FALSE, 
                                          trim_ws = TRUE, skip = 1)[,-1]

covers <- c()
for(r in 1:300){
  lower <- 2*mixture_ests_3c_dependent[r] - quantile(bootstrap_ests_3c_dependent[r,],
                                                     0.975, na.rm=T)
  upper <- 2*mixture_ests_3c_dependent[r] - quantile(bootstrap_ests_3c_dependent[r,],
                                                     0.025, na.rm=T)
  covers[r] <- lower <= 4.485 & upper >= 4.485
  print(r)
}

mixture_coverage_dependent[6] <- mean(covers, na.rm=T)
mixture_coverage_dependent_sd[6] <- sqrt(mean(covers, 
                                              na.rm=T)*(1-mean(covers, 
                                                               na.rm=T))/sum(!is.na(covers)))
mixture_dependent_mean[6] <- mean(mixture_ests_3c_dependent, na.rm=T)
mixture_dependent_sd[6] <- sd(mixture_ests_3c_dependent, na.rm=T)/sqrt(sum(!is.na(mixture_ests_3c_dependent)))




#########################################################
############   Tables of results   ######################
#########################################################

### Table 1 (label independent random effects)
df <- data.frame(assumptions = rep(c("(A1), (A2), (A3)", "(A2), (A3)", "(A1), (A2)"), each=2),
                 normal = rep(c("yes", "no"), 3),
                 prevalence_bias = round(prevalence_means, 2) - rep(0.4, 6),
                 prevalence_sd = sprintf("(%.3f)", prevalence_sd),
                 prevalence_coverage = round(prevalence_coverage, 2),
                 prevalence_coverage_sd = sprintf("(%.3f)", 
                                                  prevalence_coverage_sd),
                 mixed_effects_bias = round(me_fixed_mean - c(3, 3.485, 3, 3.485, 4, 4.485), 2),
                 mixed_effects_sd = sprintf("(%.3f)",
                                            me_fixed_sd),
                 mixed_effects_coverage = round(me_coverage_fixed, 2),
                 mixed_effects_coverage_sd = sprintf("(%.3f)", 
                                                     me_coverage_fixed_sd),
                 mixture_bias = round(mixture_fixed_mean - c(3, 3.485, 3, 3.485, 4, 4.485), 2),
                 mixture_sd = sprintf("(%.3f)", mixture_fixed_sd),
                 mixture_coverage = round(mixture_coverage_fixed, 2),
                 mixture_coverage_sd = sprintf("(%.3f)", 
                                               mixture_coverage_fixed_sd)) %>%
  unite("prevalence_bias", prevalence_bias:prevalence_sd, sep=" ") %>%
  unite("prevalence_coverage", prevalence_coverage:prevalence_coverage_sd, 
        sep=" ") %>%
  unite("mixed_effects_bias", mixed_effects_bias:mixed_effects_sd, 
        sep=" ") %>%
  unite("mixed_effects_coverage", mixed_effects_coverage:mixed_effects_coverage_sd, 
        sep=" ") %>%
  unite("mixture_bias", mixture_bias:mixture_sd, 
        sep=" ") %>%
  unite("mixture_coverage", mixture_coverage:mixture_coverage_sd, 
        sep=" ")

df$assumptions <- as.character(df$assumptions)

rle.lengths <- rle(df$assumptions)$lengths
first <- !duplicated(df$assumptions)
df$assumptions[!first] <- ""

# define appearance of \multirow
df$assumptions[first] <-
  paste0("\\midrule\\multirow{", rle.lengths, "}{*}{\\textbf{", df[[1]][first], "}}")


df %>%
  xtable(align = "cc|c|c|c|c|c|c|c",
         caption = "Coverage of bootstrap confidence intervals in simulated data.") %>%
  print(include.rownames=F, include.colnames=F,
        floating=T, hline.after=NULL,
        sanitize.text.function = force,
        add.to.row = list(pos = list(-1,
                                     6),
                          command = c(paste("\\toprule \n",  # NEW row
                                            "\\multirow{2}{*}{Assumptions} & \\multirow{2}{*}{Normal?} & \\multicolumn{2}{c|}{Prevalence} & \\multicolumn{2}{c|}{Mixed Effects Model} & \\multicolumn{2}{c|}{Mixture Model} \\\\ \n",
                                            "& & Mean & Coverage & Mean & Coverage & Mean & Coverage \\\\ \n"
                          ),
                          paste("\\bottomrule \n"  # paste is used as it is more flexible regarding adding lines
                          )
                          )
        ))



### Table 2 (label dependent random effects)
df <- data.frame(assumptions = rep(c("(A1), (A2), (A3)", "(A2), (A3)", "(A1), (A2)"), each=2),
                 normal = rep(c("yes", "no"), 3),
                 mixed_effects_bias = round(me_mean - c(3, 3.485, 3, 3.485, 4, 4.485), 2),
                 mixed_effects_sd = sprintf("(%.3f)",
                                            me_sd),
                 mixed_effects_coverage = round(me_coverage, 2),
                 mixed_effects_coverage_sd = sprintf("(%.3f)", 
                                                     me_coverage_sd),
                 mixed_effects_bias_adj = round(me_adj_mean - c(3, 3.485, 3, 3.485, 4, 4.485), 2),
                 mixed_effects_sd_adj = sprintf("(%.3f)",
                                                me_adj_sd),
                 mixed_effects_coverage_adj = round(me_coverage_adj, 2),
                 mixed_effects_coverage_sd_adj = sprintf("(%.3f)", 
                                                         me_coverage_adj_sd),
                 mixture_bias = round(mixture_dependent_mean - c(3, 3.485, 3, 3.485, 4, 4.485), 2),
                 mixture_sd = sprintf("(%.3f)", mixture_dependent_sd),
                 mixture_coverage = round(mixture_coverage_dependent, 2),
                 mixture_coverage_sd = sprintf("(%.3f)", 
                                               mixture_coverage_dependent_sd)) %>%
  unite("mixed_effects_bias", mixed_effects_bias:mixed_effects_sd, 
        sep=" ") %>%
  unite("mixed_effects_coverage", mixed_effects_coverage:mixed_effects_coverage_sd, 
        sep=" ") %>%
  unite("mixed_effects_bias_adj", mixed_effects_bias_adj:mixed_effects_sd_adj, 
        sep=" ") %>%
  unite("mixed_effects_coverage_adj", mixed_effects_coverage_adj:mixed_effects_coverage_sd_adj, 
        sep=" ") %>%
  unite("mixture_bias", mixture_bias:mixture_sd, 
        sep=" ") %>%
  unite("mixture_coverage", mixture_coverage:mixture_coverage_sd, 
        sep=" ")

df$assumptions <- as.character(df$assumptions)

rle.lengths <- rle(df$assumptions)$lengths
first <- !duplicated(df$assumptions)
df$assumptions[!first] <- ""

# define appearance of \multirow
df$assumptions[first] <-
  paste0("\\midrule\\multirow{", rle.lengths, "}{*}{\\textbf{", df[[1]][first], "}}")


df %>%
  xtable(align = "cc|c|c|c|c|c|c|c",
         caption = "Coverage of bootstrap confidence intervals in simulated data.") %>%
  print(include.rownames=F, include.colnames=F,
        floating=T, hline.after=NULL,
        sanitize.text.function = force,
        add.to.row = list(pos = list(-1,
                                     6),
                          command = c(paste("\\toprule \n",  # NEW row
                                            "\\multirow{2}{*}{Assumptions} & \\multirow{2}{*}{Normal?} & \\multicolumn{2}{c|}{Mixed Effects Model} & \\multicolumn{2}{c|}{Mixed Effects, Variance Adjustment} & \\multicolumn{2}{c|}{Mixture Model} \\\\ \n",
                                            "& & Mean & Coverage & Mean & Coverage & Mean & Coverage \\\\ \n"
                          ),
                          paste("\\bottomrule \n"  # paste is used as it is more flexible regarding adding lines
                          )
                          )
        ))
