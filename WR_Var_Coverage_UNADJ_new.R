rm(list=ls())
library(MASS)
library(nnet)
library(ggplot2)
library(parallel)
library(foreach)
library(doParallel)

#############################################################
######### Variance Coverage of unadjusted_new estimator###########
#############################################################
var_ratio = function(uR,uS,sigR2,sigS2,cov_RS){
  WR_appro_var = (uR^2/uS^2)*(sigR2/uR^2-2*cov_RS/(uR*uS)+sigS2/uS^2)
  return(WR_appro_var)
}

args <- commandArgs(trailingOnly = TRUE)
k <- as.integer(args[1])
if (is.na(k)) k <- 1
paste("Scenario:",k)

numCores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK",8))
if (is.na(numCores)) numCores <- 8
#registerDoMC(cores=numCores)

# define scenarios
scenarios <- read.table("C:/Users/user/Dropbox/research/covariate adjustment win ratio/WR_Params.txt", header=TRUE, sep="")
scenarios <- subset(scenarios, scenario == k)

scenario <- k
sds <- c(scenarios$sd_x1, scenarios$sd_x2, scenarios$sd_x3)
bern_param <- c(scenarios$param_x4, scenarios$param_x5, scenarios$param_x6)


trt_eff1 <- scenarios$trt_eff1
bi_trt <- c(1,-1,1,-1,1,-1)*scenarios$bi_trt
bi_ctrl <- c(1,-1,1,-1,1,-1)*scenarios$bi_ctrl

b01 <- 1
b02 <- 0.05
step_size <- scenarios$step_size #(sample size increase step)
sim_num <- 1000 #simulation iteration

start_time <- Sys.time()
##################################Simulation##################################
#vector of outcomes as a factor with ordered levels
outcomes_3lvl <- factor(c("first", "second", "third"), 
                        levels = c("first", "second", "third"), 
                        ordered = TRUE)

#Order: first < second < third, the larger the better
inv_logit <- function(logit) exp(logit)/(1 + exp(logit))

WP_trt_sim_unadj <- numeric(sim_num)
WP_ctrl_sim_unadj <- numeric(sim_num)
WR_sim_unadj <- numeric(sim_num)
WD_sim_unadj <- numeric(sim_num)

theory_var_unadj_trt <- numeric(sim_num)
theory_var_unadj_ctrl <- numeric(sim_num)
theory_var_WR <- numeric(sim_num)
theory_var_WD <- numeric(sim_num)
Cov_trt_ctrl = numeric(sim_num)

#for balanced design
#WP_trt_sim_true <- 0.2953701
#WP_ctrl_sim_true <- 0.2050298
#WR_sim_true <- 1.447066
#WD_sim_true <- 0.0903403
WP_trt_sim_true <-  0.2959808
WP_ctrl_sim_true <- 0.2048627
WR_sim_true <- 1.45124
WD_sim_true <- 0.0911181
n_count <- 400

for (count_temp in 1:sim_num){
  print(count_temp)
  set.seed(count_temp)
  # covariates
  x1 <- rnorm(n_count, mean = 1, sd = sds[1])
  x2 <- rnorm(n_count, mean = 0.9, sd = sds[2])
  x3 <- rnorm(n_count, mean = 0.8, sd = sds[3])
  
  x4 <- rbinom(n_count,1,bern_param[1])
  x5 <- rbinom(n_count,1,bern_param[2])
  x6 <- rbinom(n_count,1,bern_param[3])
  
  df_cov <- data.frame(x1, x2, x3, x4, x5, x6)
  treatment_assignment <- rbinom(n_count, 1, 0.7) #0.5 for balanced; 0.7 for unbalanced
  
  trt_cov <- df_cov[treatment_assignment == 1, ]
  ctrl_cov <- df_cov[treatment_assignment == 0, ]
  
  trt_cov_quad <- trt_cov^2
  ctrl_cov_quad <- ctrl_cov^2
  
  bi_trt_quad <- 2*bi_trt
  bi_ctrl_quad <- 2*bi_ctrl
  
  combined_trt_quad <- cbind(trt_cov, trt_cov_quad)
  combined_ctrl_quad <- cbind(ctrl_cov, ctrl_cov_quad)
  
  logodds1_trt <- b01 + as.matrix(combined_trt_quad) %*% c(bi_trt, bi_trt_quad) + trt_eff1
  logodds2_trt <- b02 + as.matrix(combined_trt_quad) %*% c(bi_trt, bi_trt_quad) + trt_eff1
  
  logodds1_ctrl <- b01 + as.matrix(combined_ctrl_quad) %*% c(bi_ctrl, bi_ctrl_quad)
  logodds2_ctrl <- b02 + as.matrix(combined_ctrl_quad) %*% c(bi_ctrl, bi_ctrl_quad)
  
  ## Probability Trt
  prob_2to3_trt <- inv_logit(logodds1_trt)
  prob_3_trt <- inv_logit(logodds2_trt)
  prob_1_trt <- 1 - prob_2to3_trt
  prob_2_trt <- prob_2to3_trt - prob_3_trt
  
  ## Probability Ctrl
  prob_2to3_ctrl <- inv_logit(logodds1_ctrl)
  prob_3_ctrl <- inv_logit(logodds2_ctrl)
  prob_1_ctrl <- 1 - prob_2to3_ctrl
  prob_2_ctrl <- prob_2to3_ctrl - prob_3_ctrl
  
  #generate random outcomes
  outcomes_trt <- c()
  for (i in 1:nrow(trt_cov)) {
    outcomes_trt[i] <- sample(
      outcomes_3lvl, 
      size = 1,
      prob = c(prob_1_trt[i], prob_2_trt[i], prob_3_trt[i])
    )
  }
  
  outcomes_ctrl <- c()
  for (i in 1:nrow(ctrl_cov)) {
    outcomes_ctrl[i] <- sample(
      outcomes_3lvl, 
      size = 1,
      prob = c(prob_1_ctrl[i], prob_2_ctrl[i], prob_3_ctrl[i])
    )
  }
  
  ########################Unadjusted_new Estimator#####################################
  df_trt <- data.frame(outcomes_trt)
  df_ctrl <- data.frame(outcomes_ctrl)
  
  colnames(df_trt) <- "outcomes_comb"
  colnames(df_ctrl) <- "outcomes_comb"
  
  df_comb <- rbind(df_trt, df_ctrl)
  df_comb$treatment <- c(rep(1,nrow(trt_cov)), rep(0,nrow(ctrl_cov)))
  
  df_comb$x1 <- c(x1[treatment_assignment == 1], x1[treatment_assignment == 0])
  df_comb$x2 <- c(x2[treatment_assignment == 1], x2[treatment_assignment == 0])
  df_comb$x3 <- c(x3[treatment_assignment == 1], x3[treatment_assignment == 0])
  df_comb$x4 <- c(x4[treatment_assignment == 1], x4[treatment_assignment == 0])
  df_comb$x5 <- c(x5[treatment_assignment == 1], x5[treatment_assignment == 0])
  df_comb$x6 <- c(x6[treatment_assignment == 1], x6[treatment_assignment == 0])
  
  # Define the unadjusted estimator comparison function
  compare_rows_unadj <- function(i) {
    g_hfunc_unadj1_vec <- numeric(n_count)
    g_hfunc_unadj2_vec <- numeric(n_count)
    for (j in 1:n_count) {
      if (df_comb[i,1] > df_comb[j,1]) {
        g_hfunc_unadj1_vec[j] = df_comb[i,"treatment"] * (1-df_comb[j,"treatment"])
        g_hfunc_unadj2_vec[j] = 0
      } else if (df_comb[i,1] < df_comb[j,1]) {
        g_hfunc_unadj1_vec[j] = 0
        g_hfunc_unadj2_vec[j] = df_comb[i,"treatment"] * (1-df_comb[j,"treatment"])
      }
    }
    return(list(unadj1 = g_hfunc_unadj1_vec, unadj2 = g_hfunc_unadj2_vec))
  }
  # Parallelize the outer loop
  results_unadj = apply(matrix(1:n_count,nrow=1),2,compare_rows_unadj)
  
  # Aggregate results
  g_hfunc_unadj1 <- do.call(rbind, lapply(results_unadj, function(x) x$unadj1))
  g_hfunc_unadj2 <- do.call(rbind, lapply(results_unadj, function(x) x$unadj2))
  
  # Calculate tau values
  tau1_unadj <- sum(as.vector(g_hfunc_unadj1), na.rm = TRUE)/
    (sum(df_comb[,"treatment"]) * sum(1-df_comb[,"treatment"]))
  tau2_unadj <- sum(as.vector(g_hfunc_unadj2), na.rm = TRUE)/
    (sum(df_comb[,"treatment"]) * sum(1-df_comb[,"treatment"]))
  
  #point estimate
  WP_trt_sim_unadj[count_temp] <- tau1_unadj
  WP_ctrl_sim_unadj[count_temp] <- tau2_unadj
  WR_sim_unadj[count_temp] <- tau1_unadj / tau2_unadj
  WD_sim_unadj[count_temp] <- tau1_unadj - tau2_unadj
  
  ###########unadjusted_new Theoretical Variance###########
  ###g1_unadj_trt###
  # Pre-calculate constants
  n <- nrow(df_comb)
  adjusted_sum_i_neq_j = sum(df_comb[,"treatment"])*sum(1-df_comb[,"treatment"])
  g1_constant_denominator_unadj <- (1/(n*(n-1))) * adjusted_sum_i_neq_j
  
  outcomes = df_comb$outcomes_comb
  treatment = df_comb$treatment
  
  g1_unadj_trt_fun <- function(i) {
    # Vectorized operation to calculate pairwise comparisons
    pairwise_comparisons <- (1/2) * (
      (treatment[i] * (1 - treatment) * as.numeric(outcomes[i] > outcomes)) + 
        ((1-treatment[i]) * treatment * as.numeric(outcomes[i] < outcomes))
    )
    sum_g1_unadj_trt <- sum(pairwise_comparisons) / g1_constant_denominator_unadj
    return((1/(n-1)) * sum_g1_unadj_trt)
  }
  # Parallelize the computation for g1_IPW_trt
  results_g1_unadj_trt = apply(matrix(1:n,nrow=1),2,g1_unadj_trt_fun)
  
  # Aggregate results
  g1_unadj_trt <- unlist(results_g1_unadj_trt)
  
  ###g1_IPW_ctrl###
  g1_unadj_ctrl_fun <- function(i) {
    # Vectorized operation to calculate pairwise comparisons
    pairwise_comparisons <- (1/2) * (
      (treatment[i] * (1 - treatment) * as.numeric(outcomes[i] < outcomes))+ 
        ((1-treatment[i]) * treatment * as.numeric(outcomes[i] > outcomes))
    )
    sum_g1_unadj_ctrl <- sum(pairwise_comparisons) / g1_constant_denominator_unadj
    return((1/(n-1)) * sum_g1_unadj_ctrl)
  }
  
  # Parallelize the computation for g1_IPW_trt
  results_g1_unadj_ctrl = apply(matrix(1:n,nrow=1),2,g1_unadj_ctrl_fun)
  
  # Aggregate results
  g1_unadj_ctrl <- unlist(results_g1_unadj_ctrl)
  
  unadj_Variance_trt <- function() {
    var_sum <- 0
    for (i in 1:n) {
      var_sum <- var_sum + (2 * (g1_unadj_trt[i] - WP_trt_sim_unadj[count_temp]))^2
    }
    unadj_Var_trt <- var_sum / n
    return(unadj_Var_trt)
  }
  
  unadj_Variance_ctrl <- function() {
    var_sum <- 0
    for (i in 1:n) {
      var_sum <- var_sum + (2 * (g1_unadj_ctrl[i] - WP_ctrl_sim_unadj[count_temp]))^2
    }
    unadj_Var_ctrl <- var_sum / n
    return(unadj_Var_ctrl)
  }
  
  #covariance 
  unadj_cov = function(){
    cov_sum <- 0
    for (i in 1:n) {
      temp_trt = (2 * (g1_unadj_trt[i] - WP_trt_sim_unadj[count_temp]))
      temp_ctrl = (2 * (g1_unadj_ctrl[i] - WP_ctrl_sim_unadj[count_temp]))
      cov_sum = cov_sum + temp_trt*temp_ctrl
    }
    unadj_cov_results <- cov_sum / n
    return(unadj_cov_results)
  }
  
  Var_trt_unadj_calculated <- unadj_Variance_trt()/n
  Var_ctrl_unadj_calculated <- unadj_Variance_ctrl()/n
  cov_unadj_calculated <- unadj_cov()/n
  
  #################unadjusted Estimator Output#########################
  ###Theoretical Variance###
  theory_var_unadj_trt[count_temp] <- Var_trt_unadj_calculated
  theory_var_unadj_ctrl[count_temp] <- Var_ctrl_unadj_calculated
  Cov_trt_ctrl[count_temp] <- cov_unadj_calculated
  theory_var_WR[count_temp] <- var_ratio(WP_trt_sim_unadj[count_temp],WP_ctrl_sim_unadj[count_temp],
                                         theory_var_unadj_trt[count_temp],theory_var_unadj_ctrl[count_temp],
                                         Cov_trt_ctrl[count_temp])
  theory_var_WD[count_temp] <- Var_trt_unadj_calculated+Var_ctrl_unadj_calculated-2*cov_unadj_calculated
}#open on line
end_time <- Sys.time()
run_time <- end_time - start_time
run_time

######Coverage Output######

theory_sd_unadj_trt <- sqrt(theory_var_unadj_trt)
theory_sd_unadj_ctrl <- sqrt(theory_var_unadj_ctrl)
theory_sd_WR = sqrt(theory_var_WR)
theory_sd_WD = sqrt(theory_var_WD)

WP_trt_unadj_coverage <- mean((WP_trt_sim_true > (WP_trt_sim_unadj - qnorm(0.975)*theory_sd_unadj_trt)) & 
                              (WP_trt_sim_true < (WP_trt_sim_unadj + qnorm(0.975)*theory_sd_unadj_trt)))
WP_ctrl_unadj_coverage <- mean((WP_ctrl_sim_true > (WP_ctrl_sim_unadj - qnorm(0.975)*theory_sd_unadj_ctrl)) & 
                               (WP_ctrl_sim_true < (WP_ctrl_sim_unadj + qnorm(0.975)*theory_sd_unadj_ctrl)))
WR_coverage <- mean((WR_sim_true > (WR_sim_unadj - qnorm(0.975)*theory_sd_WR)) & 
                      (WR_sim_true < (WR_sim_unadj + qnorm(0.975)*theory_sd_WR)))
WD_coverage <- mean((WD_sim_true > (WD_sim_unadj - qnorm(0.975)*theory_sd_WD)) & 
                      (WD_sim_true < (WD_sim_unadj + qnorm(0.975)*theory_sd_WD)))

all_est = cbind(WP_trt_sim_unadj,WP_ctrl_sim_unadj,WR_sim_unadj,WD_sim_unadj)
wp_true = c(WP_trt_sim_true,WP_ctrl_sim_true,WR_sim_true,WD_sim_true)
wp_est = apply(all_est,2,mean)
wp_sd = apply(all_est,2,sd)
wp_se = sqrt(apply(cbind(theory_var_unadj_trt,theory_var_unadj_ctrl,theory_var_WR,theory_var_WD),2,mean))
wp_cr = c(WP_trt_unadj_coverage,WP_ctrl_unadj_coverage,WR_coverage,WD_coverage)
res = data.frame(wp_true,wp_est,wp_sd,wp_se,wp_cr)
res
write.csv(res,file="C:/Users/user/Dropbox/research/covariate adjustment win ratio/unadjusted_unb_400.csv")  
