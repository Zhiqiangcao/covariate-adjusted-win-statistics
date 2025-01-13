rm(list=ls())
library(MASS)
library(nnet)
#library(VGAM)
library(ggplot2)
library(parallel)
library(foreach)
library(doParallel)
#numCores <- detectCores()

#############################################################
################ Variance Coverage of AIPW estimator##################
#############################################################
var_ratio = function(uR,uS,sigR2,sigS2,cov_RS){
  WR_appro_var = (uR^2/uS^2)*(sigR2/uR^2-2*cov_RS/(uR*uS)+sigS2/uS^2)
  return(WR_appro_var)
}

args<-commandArgs(trailingOnly = TRUE)
k<-as.integer(args[1])
if (is.na(k)) k <- 1
paste("Scenario:",k)

numCores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK",8))
if (is.na(numCores)) numCores <- 8
#registerDoMC(cores=numCores)

# define scenarios
scenarios <- read.table("C:/Users/82655/Dropbox/research/covariate adjustment win ratio/WR_Params.txt", header=TRUE, sep="")
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
re_count <- 1 #sample size count
sim_num <- 1000 #simulation iteration
# re_count <- scenarios$re_count #sample size count
# sim_num <- scenarios$sim_num #simulation iteration

start_time <- Sys.time()

##################################Simulation##################################
#vector of outcomes as a factor with ordered levels
outcomes_3lvl <- factor(c("first", "second", "third"), 
                        levels = c("first", "second", "third"), 
                        ordered = TRUE)

#Order: first < second < third, the larger the better
inv_logit <- function(logit) exp(logit)/(1 + exp(logit))


WP_trt_sim_AIPW <- numeric(sim_num)
WP_ctrl_sim_AIPW <- numeric(sim_num)
WR_sim_AIPW <- numeric(sim_num)
WD_sim_AIPW <- numeric(sim_num)

theory_var_AIPW_trt <- numeric(sim_num)
theory_var_AIPW_ctrl <- numeric(sim_num)
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

n_count <- 200

for (count_temp in 1:sim_num){
  print(count_temp)
  set.seed(count_temp+10000)
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
  
  ########### AIPW ###########
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
  
  # Calculate propensity scores
  PropScore <- glm(treatment~x1 + x2 + x3 + x4 + x5 + x6, 
                   data = df_comb, family=binomial)
  pi_func <- fitted(PropScore)
  
  ### mu-function
  if((!(1 %in% df_comb[df_comb$treatment==1,]$outcomes_comb) |
      !(2 %in% df_comb[df_comb$treatment==1,]$outcomes_comb) |
      !(3 %in% df_comb[df_comb$treatment==1,]$outcomes_comb))|
     (!(1 %in% df_comb[df_comb$treatment==0,]$outcomes_comb) |
      !(2 %in% df_comb[df_comb$treatment==0,]$outcomes_comb) |
      !(3 %in% df_comb[df_comb$treatment==0,]$outcomes_comb))){
    break
  }
  
  else {
    lp_start <- rep(0, 12)
    th_start <- c(-1, 1)
    start_values <- c(lp_start, th_start)
    dr_trt <- polr(factor(outcomes_comb) ~ x1 + x2 + x3 + x4 + x5 + x6 + 
                    I(x1^2) + I(x2^2) + I(x3^2) + I(x4^2) + I(x5^2) + I(x6^2), 
                   data = df_comb[df_comb$treatment == 1,],
                   start = start_values)
    cond_prob_trt <- predict(dr_trt, newdata = df_comb, type = "probs")
    
    dr_ctrl <- polr(factor(outcomes_comb) ~ x1 + x2 + x3 + x4 + x5 + x6 + 
                    I(x1^2) + I(x2^2) + I(x3^2) + I(x4^2) + I(x5^2) + I(x6^2), 
                    data = df_comb[df_comb$treatment == 0,],
                    start = start_values)
    cond_prob_ctrl <- predict(dr_ctrl, newdata = df_comb, type = "probs")
  }
  
  cl <- makeCluster(numCores - 1)  # Use all cores except one
  registerDoParallel(cl)
  
  # Define the parallelized function
  compare_rows_AIPW <- function(i, df_comb, pi_func, cond_prob_trt, cond_prob_ctrl) {
    n_count <- nrow(df_comb)
    AIPW_num1_vec <- numeric(n_count)
    AIPW_denom1_vec <- numeric(n_count)
    AIPW_num2_vec <- numeric(n_count)
    AIPW_denom2_vec <- numeric(n_count)
    AIPW_mu1_vec <- numeric(n_count)
    AIPW_mu2_vec <- numeric(n_count)
    AIPW_mu_denom1_vec <- numeric(n_count)
    AIPW_mu_denom2_vec <- numeric(n_count)
    
    for(j in 1:n_count){
      if(i != j){
        AIPW_num1_vec[j] = (df_comb[i,'treatment'] * (1 - df_comb[j,'treatment'])) * (1/(pi_func[i] * (1 - pi_func[j]))) *
          ((df_comb[i,1] > df_comb[j,1]) - 
             (cond_prob_trt[i, 2] * cond_prob_ctrl[j, 1] + 
                cond_prob_trt[i, 3] * (cond_prob_ctrl[j, 1] + cond_prob_ctrl[j, 2])))
        
        AIPW_denom1_vec[j] = (df_comb[i,'treatment'] * (1 - df_comb[j,'treatment'])) * (1/(pi_func[i] * (1 - pi_func[j])))
        
        AIPW_num2_vec[j] = (df_comb[i,'treatment'] * (1 - df_comb[j,'treatment'])) * (1/(pi_func[i] * (1 - pi_func[j]))) *
          ((df_comb[i,1] < df_comb[j,1]) - 
             (cond_prob_trt[i, 1] * cond_prob_ctrl[j, 2] + 
                cond_prob_trt[i, 1] * cond_prob_ctrl[j, 3] +
                cond_prob_trt[i, 2] * cond_prob_ctrl[j, 3]))
        
        AIPW_denom2_vec[j] = (df_comb[i,'treatment'] * (1 - df_comb[j,'treatment'])) * (1/(pi_func[i] * (1 - pi_func[j])))
        
        AIPW_mu1_vec[j] = 1 * (cond_prob_trt[i, 2] * cond_prob_ctrl[j, 1] + 
                                 cond_prob_trt[i, 3] * (cond_prob_ctrl[j, 1] + cond_prob_ctrl[j, 2]))
        
        AIPW_mu2_vec[j] = 1 * (cond_prob_trt[i, 1] * cond_prob_ctrl[j, 2] + 
                                 cond_prob_trt[i, 1] * cond_prob_ctrl[j, 3] + 
                                 cond_prob_trt[i, 2] * cond_prob_ctrl[j, 3])
        AIPW_mu_denom1_vec[j] <- 1
        AIPW_mu_denom2_vec[j] <- 1
      }
    }
    return(list(AIPW_num1 = AIPW_num1_vec, AIPW_denom1 = AIPW_denom1_vec, 
                AIPW_num2 = AIPW_num2_vec, AIPW_denom2 = AIPW_denom2_vec,
                AIPW_mu1 = AIPW_mu1_vec, AIPW_mu_denom1 = AIPW_mu_denom1_vec,
                AIPW_mu2 = AIPW_mu2_vec, AIPW_mu_denom2 = AIPW_mu_denom2_vec))
    
  }
  results_AIPW <- foreach(i=1:n_count) %dopar% {
    compare_rows_AIPW(i, df_comb, pi_func, cond_prob_trt, cond_prob_ctrl)
  }
  
  # Combine the results
  tau_AIPW_num1 <- do.call(rbind, lapply(results_AIPW, function(x) x$AIPW_num1))
  tau_AIPW_denom1 <- do.call(rbind, lapply(results_AIPW, function(x) x$AIPW_denom1))
  tau_AIPW_num2 <- do.call(rbind, lapply(results_AIPW, function(x) x$AIPW_num2))
  tau_AIPW_denom2 <- do.call(rbind, lapply(results_AIPW, function(x) x$AIPW_denom2))
  tau_AIPW_mu1 <- do.call(rbind, lapply(results_AIPW, function(x) x$AIPW_mu1))
  tau_AIPW_mu_denom1 <- do.call(rbind, lapply(results_AIPW, function(x) x$AIPW_mu_denom1))
  tau_AIPW_mu2 <- do.call(rbind, lapply(results_AIPW, function(x) x$AIPW_mu2))
  tau_AIPW_mu_denom2 <- do.call(rbind, lapply(results_AIPW, function(x) x$AIPW_mu_denom2))
  
  tau1_AIPW <- (sum(as.vector(tau_AIPW_num1), na.rm = T) / sum(as.vector(tau_AIPW_denom1), na.rm = T)) +
    (sum(as.vector(tau_AIPW_mu1), na.rm = T)/sum(as.vector(tau_AIPW_mu_denom1), na.rm = T))
  
  tau2_AIPW <- (sum(as.vector(tau_AIPW_num2), na.rm = T) / sum(as.vector(tau_AIPW_denom2), na.rm = T)) +
    (sum(as.vector(tau_AIPW_mu2), na.rm = T)/sum(as.vector(tau_AIPW_mu_denom2), na.rm = T))
  
  WP_trt_sim_AIPW[count_temp] <- tau1_AIPW
  WP_ctrl_sim_AIPW[count_temp] <- tau2_AIPW
  WR_sim_AIPW[count_temp] <- tau1_AIPW / tau2_AIPW
  WD_sim_AIPW[count_temp] <- tau1_AIPW - tau2_AIPW
  
  stopCluster(cl)
  
  ########### AIPW Theoretical Variance ###########
  mu_func_trt <- function(a, b, cond_prob_trt, cond_prob_ctrl) {
    mu_a_b <- (
      cond_prob_trt[a, 2] * cond_prob_ctrl[b, 1] + 
        cond_prob_trt[a, 3] * (cond_prob_ctrl[b, 1] + cond_prob_ctrl[b, 2])
    )
    return(mu_a_b)
  }
  
  mu_func_ctrl <- function(a, b, cond_prob_trt, cond_prob_ctrl) {
    mu_a_b <- (
      cond_prob_trt[a, 2] * cond_prob_ctrl[b, 3]+
        cond_prob_trt[a, 1] * (cond_prob_ctrl[b, 2] + cond_prob_ctrl[b, 3])
    )
    return(mu_a_b)
  }
  
  ###g1_AIPW_trt###
  n <- nrow(df_comb)
  outcomes = df_comb$outcomes_comb
  treatment = df_comb$treatment
  adjusted_sum_i_neq_j = sum(df_comb[,"treatment"]/pi_func) * sum((1-df_comb[,"treatment"])/(1-pi_func))
  g1_constant_denominator_AIPW <- (1/(n*(n-1))) * adjusted_sum_i_neq_j
  
  g1_AIPW_trt_fun <- function(i) {
    #n <- nrow(df_comb)
    pairwise_comparisons_AIPW_1 <- (1/2) * (
      ((treatment[i] * (1 - treatment) / (pi_func[i] * (1 - pi_func))) * 
         (as.numeric(outcomes[i] > outcomes) - sapply(1:n, function(b) mu_func_trt(i, b, cond_prob_trt, cond_prob_ctrl)))) +
        (((1 - treatment[i]) * treatment / (pi_func * (1 - pi_func[i]))) * 
           (as.numeric(outcomes > outcomes[i]) - sapply(1:n, function(a) mu_func_trt(a, i, cond_prob_trt, cond_prob_ctrl))))
    )
    pairwise_comparisons_AIPW_2 <- (1/2) * (sapply(1:n, function(b) mu_func_trt(i, b, cond_prob_trt, cond_prob_ctrl)) +
                                              sapply(1:n, function(a) mu_func_trt(a, i, cond_prob_trt, cond_prob_ctrl)))
    
    sum_g1_AIPW_trt <- (sum(pairwise_comparisons_AIPW_1) / g1_constant_denominator_AIPW) + 
      sum(pairwise_comparisons_AIPW_2)
    return((1/(n-1)) * sum_g1_AIPW_trt)
  }
  
  # Parallelize the computation for g1_AIPW_trt
  #results_g1_AIPW_trt <- mclapply(1:n, g1_AIPW_trt_fun, df_comb$outcomes_comb, df_comb$treatment, mc.cores = numCores - 1)
  results_g1_AIPW_trt <- apply(matrix(1:n,nrow=1),2,g1_AIPW_trt_fun)
  g1_AIPW_trt <- unlist(results_g1_AIPW_trt)
  
  
  ###g1_AIPW_ctrl###
  g1_AIPW_ctrl_fun <- function(i) {
    #n <- nrow(df_comb)
    pairwise_comparisons_AIPW_1 <- (1/2) * (
      ((treatment[i] * (1 - treatment) / (pi_func[i] * (1 - pi_func))) * 
         (as.numeric(outcomes[i] < outcomes) - sapply(1:n, function(b) mu_func_ctrl(i, b, cond_prob_trt, cond_prob_ctrl)))) +
        (((1 - treatment[i]) * treatment / (pi_func * (1 - pi_func[i]))) * 
           (as.numeric(outcomes < outcomes[i]) - sapply(1:n, function(a) mu_func_ctrl(a, i, cond_prob_trt, cond_prob_ctrl))))
    )
    pairwise_comparisons_AIPW_2 <- (1/2) * (sapply(1:n, function(b) mu_func_ctrl(i, b, cond_prob_trt, cond_prob_ctrl)) +
                                              sapply(1:n, function(a) mu_func_ctrl(a, i, cond_prob_trt, cond_prob_ctrl)))

    sum_g1_AIPW_ctrl <- (sum(pairwise_comparisons_AIPW_1) / g1_constant_denominator_AIPW) + 
      sum(pairwise_comparisons_AIPW_2)
    return((1/(n-1)) * sum_g1_AIPW_ctrl)
  }
  #results_g1_AIPW_ctrl <- mclapply(1:n, g1_AIPW_ctrl_fun, df_comb$outcomes_comb, df_comb$treatment, mc.cores = numCores - 1)
  results_g1_AIPW_ctrl <- apply(matrix(1:n,nrow=1),2,g1_AIPW_ctrl_fun)
  g1_AIPW_ctrl <- unlist(results_g1_AIPW_ctrl)
  
  ### B ###
  ### B1_trt & B1_ctrl
  B1_AIPW <- function(df, pi_func, cond_prob_trt, cond_prob_ctrl) {
    covariate_names <- c("x1", "x2", "x3", "x4", "x5", "x6")
    B1_value_f = function(j){
      Xj <- c(1,as.matrix(df[j, covariate_names]))
      numerator_vec1 <- numerator_vec2 <- numeric(length(covariate_names)+1)
      # Loop through all control indices and calculate contributions to numerator and denominator
      for (i in 1:n) {
        Xi <- c(1,as.matrix(df[i, covariate_names]))
        e_Xi <- pi_func[i]
        e_Xj <- pi_func[j]
        mu_ij_trt <- mu_func_trt(i, j, cond_prob_trt, cond_prob_ctrl)
        term1 <- (e_Xj * Xj - (1 - e_Xi) * Xi) * (df$treatment[i] * (1 - df$treatment[j])) *
          (as.numeric(df$outcomes[i] > df$outcomes[j]) - mu_ij_trt) / (e_Xi * (1 - e_Xj))
        numerator_vec1 <- numerator_vec1 + term1
        mu_ij_ctrl <- mu_func_ctrl(i, j, cond_prob_trt, cond_prob_ctrl)
        term2 <- (e_Xj * Xj - (1 - e_Xi) * Xi) * (df$treatment[i] * (1 - df$treatment[j])) *
          (as.numeric(df$outcomes[i] < df$outcomes[j]) - mu_ij_ctrl) / (e_Xi * (1 - e_Xj))
        numerator_vec2 <- numerator_vec2 + term2
      }
      return(list(numerator1 = numerator_vec1,numerator2 = numerator_vec2))
    }
    B1_values <- apply(matrix(1:n,nrow=1),2,B1_value_f)
    
    total_numerator1 <- Reduce("+", lapply(B1_values, `[[`, "numerator1"))
    total_numerator2 <- Reduce("+", lapply(B1_values, `[[`, "numerator2"))
    #denominator <- sum(sapply(B1_values, `[[`, "denominator"))
    denominator = adjusted_sum_i_neq_j
    B1_trt <- total_numerator1 / denominator
    B1_ctrl <- total_numerator2 / denominator
    return(list(B1_trt = B1_trt,B1_ctrl = B1_ctrl))
  }
  
  ### B2_trt & B2_ctrl ### (Dim: 1*1)
  B2_AIPW <- function(df, pi_func, cond_prob_trt, cond_prob_ctrl) {
    B2_value_f = function(j){
      numerator_vec1 <- numerator_vec2 <- 0
      # Loop through all control indices and calculate contributions to numerator and denominator
      for (i in 1:n) {
        e_Xi <- pi_func[i]
        e_Xj <- pi_func[j]
        mu_ij_trt <- mu_func_trt(i, j, cond_prob_trt, cond_prob_ctrl)
        term1 <- (df$treatment[i] * (1 - df$treatment[j])) * 
          (as.numeric(df$outcomes[i] > df$outcomes[j]) - mu_ij_trt) / (e_Xi * (1 - e_Xj))
        numerator_vec1 <- numerator_vec1 + term1
        mu_ij_ctrl <- mu_func_ctrl(i, j, cond_prob_trt, cond_prob_ctrl)
        term2 <- (df$treatment[i] * (1 - df$treatment[j])) * 
          (as.numeric(df$outcomes[i] < df$outcomes[j]) - mu_ij_ctrl) / (e_Xi * (1 - e_Xj))
        numerator_vec2 <- numerator_vec2 + term2
      }
      return(list(numerator1 = numerator_vec1, numerator2 = numerator_vec2))
    }
    B2_values <- apply(matrix(1:n,nrow=1),2,B2_value_f)
    
    total_numerator1 <- Reduce("+", lapply(B2_values, `[[`, "numerator1"))
    total_numerator2 <- Reduce("+", lapply(B2_values, `[[`, "numerator2"))
    denominator <- adjusted_sum_i_neq_j
    B2_trt <- total_numerator1 / denominator
    B2_ctrl <- total_numerator2 / denominator
    return(list(B2_trt = B2_trt,B2_ctrl = B2_ctrl))
  }
  
  ### B3 ### (Dim: 7*1)
  B3_AIPW_fun <- function(df, pi_func) {
    covariate_names <- c("x1", "x2", "x3", "x4", "x5", "x6")
    B3_term1_f = function(j){
      Xj <- c(1,as.matrix(df[j, covariate_names]))
      term1_vec <- numeric(length(covariate_names)+1)
      for (i in 1:n) {
        if (i != j) {
          Xi <- c(1,as.matrix(df[i, covariate_names]))
          term1_vec <- term1_vec + (-df$treatment[i]) * (1-pi_func[i]) * Xi / pi_func[i]
        }
      }
      term1_vec <- term1_vec * (1-df$treatment[j])/(1-pi_func[j])
      return(list(term1 = term1_vec))
    }
    B3_term1 <- apply(matrix(1:n,nrow=1),2,B3_term1_f)
    
    B3_term2_f = function(j){
      Xj <- c(1,as.matrix(df[j, covariate_names]))
      term2_vec <- numeric(length(covariate_names)+1)
      for (i in 1:n) {
        if (i != j) {
          term2_vec <- term2_vec + (df$treatment[i] / pi_func[i])
        }
      }
      term2_vec <- term2_vec * (1 - df$treatment[j]) * pi_func[j] * 
        Xj / (1 - pi_func[j])
      return(list(term2 = term2_vec))
    }
    B3_term2 <- apply(matrix(1:n,nrow=1),2,B3_term2_f)
    
    total_term1 <- Reduce("+", lapply(B3_term1, `[[`, "term1"))
    total_term2 <- Reduce("+", lapply(B3_term2, `[[`, "term2"))
    return(total_term1 + total_term2)  
  }
  
  ### B4 ### (Dim: 1*1)
  B4_AIPW_fun <- function(df, pi_func) {
    B4_term_f = function(j){
      term <- 0
      for (i in 1:n) {
        if (i != j) {
          term <- term + (df$treatment[i] / pi_func[i])
        }
      }
      term <- term * ((1 - df$treatment[j]) / (1 - pi_func[j]))
      return(list(B4term = term))
    }
    B4_term <- apply(matrix(1:n,nrow=1),2,B4_term_f)
    return(Reduce("+", lapply(B4_term, `[[`, "B4term")))
  }
  #note: B4 is equal to adjusted_sum_i_neq_j
  
  ### Bn ### (Dim: 6*1)
  B1 = B1_AIPW(df_comb, pi_func, cond_prob_trt, cond_prob_ctrl)
  B1_trt = matrix(B1$B1_trt,nrow=1)
  B1_ctrl = matrix(B1$B1_ctrl,nrow=1)
  B2 = B2_AIPW(df_comb, pi_func, cond_prob_trt, cond_prob_ctrl)
  B2_trt = B2$B2_trt
  B2_ctrl = B2$B2_ctrl
  B3 = B3_AIPW_fun(df_comb, pi_func)
  B4 = adjusted_sum_i_neq_j
  
  Bn_AIPW_trt_fun <- function(){
    Bn_trt <- B1_trt - B3*B2_trt/B4
    return(t(Bn_trt))
  }
  
  Bn_AIPW_ctrl_fun <- function(){
    Bn_ctrl <- B1_ctrl - B3*B2_ctrl/B4
    return(t(Bn_ctrl))
  }
  
  ### E_beta0beta0
  Ebeta0 <- function() {
    n <- nrow(df_comb)
    X <- as.matrix(cbind(1,df_comb[, c("x1", "x2", "x3", "x4", "x5", "x6")]))
    ps <- pi_func
    Ebeta = crossprod(sqrt(ps*(1-ps)) * X) / n
    return(Ebeta)
  }
  
  ### l_beta ### (Dim: 7*1)
  inv_Ebeta0 = solve(Ebeta0())
  l_beta_fun <- function(i) {
    X_i <- as.numeric(c(1,df_comb[i, c("x1", "x2", "x3", "x4", "x5", "x6")]))
    Z_i <- df_comb[i, "treatment"]
    pi_i <- pi_func[i]
    l_beta <- inv_Ebeta0 %*% X_i * (Z_i - pi_i)
    return(l_beta)
  }
  
  ### E_theta_theta
  thresholds_trt <- dr_trt$zeta
  thresholds_ctrl <- dr_ctrl$zeta
  
  E_thetatheta_fun <- function(df, cond_prob_trt, cond_prob_ctrl, thresholds_trt, thresholds_ctrl) {
    #n <- nrow(df)
    F_matrices <- lapply(1:n, function(i) {
      Xi0 <- as.matrix(df[i, c("x1", "x2", "x3", "x4", "x5", "x6")])
      Xi2 = Xi0^2  #quadratic terms
      Xi = t(cbind(Xi0,Xi2)) 
      thresholds <- if(df$treatment[i] == 1) thresholds_trt else thresholds_ctrl
      cond_prob <- if(df$treatment[i] == 1) cond_prob_trt else cond_prob_ctrl
      
      # Calculate delta values based on the outcome for individual i
      delta <- as.numeric(df[i, ]$outcomes_comb == 1:3)
      
      # Calculate the c values for individual i
      c1 <- cond_prob[i, 1]
      c2 <- cond_prob[i, 1] + cond_prob[i, 2]
      
      # Calculate the F matrix components for individual i
      F11 <- -(delta[1] + delta[2]) * (c1 - c1^2) - delta[2] * exp(thresholds[2] - thresholds[1]) *
        ((exp(thresholds[2] - thresholds[1]) - 1)^(-2))
      F12 <- delta[2] * exp(thresholds[2] - thresholds[1]) * ((exp(thresholds[2] - thresholds[1]) - 1)^(-2))
      F21 <- F12
      F13 <- (delta[1] + delta[2]) * (c1 - c1^2) * t(Xi)
      F31 <- t(F13)
      F22 <- -(delta[2] + delta[3]) * (c2 - c2^2) - delta[2] * exp(thresholds[1] - thresholds[2]) * 
        (1 - exp(thresholds[1] - thresholds[2]))^(-2)
      F23 <- (delta[2] + delta[3]) * (c2 - c2^2) * t(Xi)
      F32 <- t(F23)
      F33 <- -((delta[1] + delta[2]) * (c1 - c1^2) + (delta[2] + delta[3]) * (c2 - c2^2)) * (Xi %*% t(Xi))
      
      # Construct the F matrix for individual i
      top_left <- rbind(cbind(F11, F12), cbind(F21, F22))
      top_2row <- cbind(top_left, rbind(F13, F23))
      bottom_row <- cbind(cbind(F31, F32), F33)
      F_matrix <- as.matrix(rbind(top_2row, bottom_row))
      return(F_matrix)
    })
    # Compute the average of the F matrices across all individuals
    E_thetatheta <- Reduce("+", F_matrices) / n
    return(E_thetatheta)
  }
  E_thetatheta <- E_thetatheta_fun(df_comb, cond_prob_trt, cond_prob_ctrl, thresholds_trt, thresholds_ctrl)
  
  ### mu(Xi, Xj) 
  mu_tilde_AIPW_fun <- function(i, j, df, cond_prob_trt, cond_prob_ctrl) {
    # Extract Xi and Xj as row vectors
    Xi0 = df[i, c("x1", "x2", "x3", "x4", "x5", "x6")]
    Xi2 = Xi0^2  #quadratic terms
    Xi = t(cbind(Xi0,Xi2))
    Xj0 = df[j, c("x1", "x2", "x3", "x4", "x5", "x6")]
    Xj2 = Xj0^2
    Xj <- t(cbind(Xj0,Xj2))
    
    # Calculate the c values for individuals i and j
    ci_1 <- cond_prob_trt[i, 1]
    ci_2 <- cond_prob_trt[i, 1] + cond_prob_trt[i, 2]
    cj_1 <- cond_prob_ctrl[j, 1]
    cj_2 <- cond_prob_ctrl[j, 1] + cond_prob_ctrl[j, 2]
    
    # Construct the mu_tilde vector according to the formula
    mu_tilde_vector_1 <- -ci_1 * (1 - ci_1) * cj_1 + cj_1 * (1 - cj_1) * (ci_2 - ci_1)
    mu_tilde_vector_2 <- -ci_2 * (1 - ci_2) * (cj_2 - cj_1) + (1 - ci_2) * cj_2 * (1 - cj_2)
    mu_tilde_vector_3 <- (cj_1 * (ci_1 - ci_1^2) + (ci_2 - ci_2^2) * (cj_2 - cj_1)) * Xi - 
      ((ci_2 - ci_1) * (cj_1 - cj_2^2) + (1 - ci_2) * (cj_2 - cj_2^2)) * Xj
    mu_tilde_vector <- rbind(mu_tilde_vector_1, mu_tilde_vector_2,
                             mu_tilde_vector_3)
    return(mu_tilde_vector)
  }
  
  
  ### Cn ### (Dim: 8*1)
  Cn_AIPW_fun <- function(df, pi_func, cond_prob_trt, cond_prob_ctrl) {
    #n <- nrow(df)
    # Calculate the first term of Cn
    Cn_first_term_f = function(j){
      numerator_vec <- 0
      #denominator <- 0
      # Loop through all control indices and calculate contributions to numerator and denominator
      for (i in 1:n) {
        e_Xi <- pi_func[i]
        e_Xj <- pi_func[j]
        mu_tilde_ij <- mu_tilde_AIPW_fun(i, j, df, cond_prob_trt, cond_prob_ctrl)
        term <- (df$treatment[i] * (1 - df$treatment[j])) * (-mu_tilde_ij) /
          (e_Xi * (1 - e_Xj))
        numerator_vec <- numerator_vec + term
      }
      #denominator <- denominator * (1-df$treatment[j])/(1-pi_func[j])
      return(list(numerator = numerator_vec))
    }
    Cn_first_term <- apply(matrix(1:n,nrow=1),2,Cn_first_term_f)
    
    total_numerator <- Reduce("+", lapply(Cn_first_term, `[[`, "numerator"))
    #denominator <- sum(sapply(Cn_first_term, `[[`, "denominator"))
    denominator = adjusted_sum_i_neq_j
    first_term <- total_numerator / denominator
    
    # Second term of Cn
    second_term = matrix(0,14,1)
    for(i in 1:n){
      for(j in 1:n){
        if (i != j) {
          res = mu_tilde_AIPW_fun(i, j, df, cond_prob_trt, cond_prob_ctrl)
        }else{
          res = matrix(0,14,1)
        }
        second_term = second_term + res
      }
    }
    
    second_term <- second_term / (n * (n - 1))
    Cn_value <- first_term + second_term
    return(Cn_value)
  }
  
  # Use the function to calculate Cn
  Cn_value <- Cn_AIPW_fun(df_comb, pi_func, cond_prob_trt, cond_prob_ctrl)
  
  ### l_tilde_theta 
  E_thetatheta_inv <- ginv(E_thetatheta)
  
  l_theta_AIPW_fun <- function(i, df_comb, cond_prob_trt, cond_prob_ctrl, 
                               thresholds_trt, thresholds_ctrl, E_thetatheta_inv) {
    Xi0 <- as.matrix(df_comb[i, c("x1", "x2", "x3", "x4", "x5", "x6")])
    Xi2 = Xi0^2  #quadratic terms
    Xi = t(cbind(Xi0,Xi2)) 
    #Yi <- df_comb$outcomes[i]
    cond_prob <- if (df_comb$treatment[i] == 1) cond_prob_trt else cond_prob_ctrl
    
    c1 <- cond_prob[i, 1]
    c2 <- c1 + cond_prob[i, 2]
    delta <- as.numeric(df_comb[i, ]$outcomes_comb == 1:3)
    thresholds <- if (df_comb$treatment[i] == 1) thresholds_trt else thresholds_ctrl
    
    # Construct the l_tilde vector according to the formula
    l_tilde_vector_1 <- delta[1] - (delta[1] + delta[2]) * c1 - delta[2] * (exp(thresholds[2] - thresholds[1])-1)^(-1)
    l_tilde_vector_2 <- -(delta[2] + delta[3]) * c2 + delta[2] * (1 - exp(thresholds[1] - thresholds[2]))^(-1)
    l_tilde_vector_3 <- (-(delta[1] + delta[2]) + (delta[1] + delta[2]) * c1 + (delta[2] + delta[3]) * c2) * Xi
    l_tilde_vector <- rbind(l_tilde_vector_1, l_tilde_vector_2, l_tilde_vector_3)
    l_tilde_value <- E_thetatheta_inv %*% l_tilde_vector
    return(l_tilde_value)
  }
  
  AIPW_Variance_trt <- function() {
    var_sum <- 0
    Bn_t <- t(Bn_AIPW_trt_fun())
    Cn_t <- t(Cn_value)
    for (i in 1:n) {
      l_theta <- l_theta_AIPW_fun(i, df_comb, cond_prob_trt, cond_prob_ctrl, 
                                  thresholds_trt, thresholds_ctrl, E_thetatheta_inv)
      var_sum <- var_sum + (2 * (g1_AIPW_trt[i] - WP_trt_sim_AIPW[count_temp]) + 
                              Bn_t %*% l_beta_fun(i) +
                              Cn_t %*% l_theta)^2
    }
    AIPW_Var_trt <- var_sum / n
    return(AIPW_Var_trt)
  }
  
  AIPW_Variance_ctrl <- function() {
    var_sum <- 0
    Bn_t <- t(Bn_AIPW_ctrl_fun())
    Cn_t <- t(Cn_value)
    for (i in 1:n) {
      l_theta <- l_theta_AIPW_fun(i, df_comb, cond_prob_trt, cond_prob_ctrl, 
                                  thresholds_trt, thresholds_ctrl, E_thetatheta_inv)
      var_sum <- var_sum + (2 * (g1_AIPW_ctrl[i] - WP_ctrl_sim_AIPW[count_temp]) + 
                              Bn_t %*% l_beta_fun(i) +
                              Cn_t %*% l_theta)^2
    }
    AIPW_Var_ctrl <- var_sum / n
    return(AIPW_Var_ctrl)
  }
  
  #covariance 
  AIPW_cov = function(){
    cov_sum <- 0
    Bn_t_trt <- t(Bn_AIPW_trt_fun())
    Bn_t_ctrl <- t(Bn_AIPW_ctrl_fun())
    Cn_t <- t(Cn_value)
    for (i in 1:n) {
      l_theta <- l_theta_AIPW_fun(i, df_comb, cond_prob_trt, cond_prob_ctrl, 
                                  thresholds_trt, thresholds_ctrl, E_thetatheta_inv)
      l_beta <- l_beta_fun(i)
      temp_trt = (2 * (g1_AIPW_trt[i] - WP_trt_sim_AIPW[count_temp]) + 
                    Bn_t_trt %*% l_beta + Cn_t %*% l_theta)
      temp_ctrl = (2 * (g1_AIPW_ctrl[i] - WP_ctrl_sim_AIPW[count_temp]) + 
                     Bn_t_ctrl %*% l_beta + Cn_t %*% l_theta)
      cov_sum <- cov_sum + temp_trt*temp_ctrl
    }
    AIPW_cov_results <- cov_sum / n
    return(AIPW_cov_results)
  }
  
  
  Var_trt_AIPW_calculated <- AIPW_Variance_trt()/nrow(df_comb)
  Var_ctrl_AIPW_calculated <- AIPW_Variance_ctrl()/nrow(df_comb)
  cov_AIPW_calculated <- AIPW_cov()/nrow(df_comb)
  
  #################AIPW Estimator Output#########################
  ###Theoretical Variance###
  theory_var_AIPW_trt[count_temp] <- Var_trt_AIPW_calculated
  theory_var_AIPW_ctrl[count_temp] <- Var_ctrl_AIPW_calculated
  Cov_trt_ctrl[count_temp] <- cov_AIPW_calculated
  theory_var_WR[count_temp] <- var_ratio(WP_trt_sim_AIPW[count_temp],WP_ctrl_sim_AIPW[count_temp],
                                         theory_var_AIPW_trt[count_temp],theory_var_AIPW_ctrl[count_temp],
                                         Cov_trt_ctrl[count_temp])
  theory_var_WD[count_temp] <- Var_trt_AIPW_calculated+Var_ctrl_AIPW_calculated-2*cov_AIPW_calculated
}#open on line

end_time <- Sys.time()
run_time <- end_time - start_time
run_time
#Time difference of 1.709125 mins per 10 simulations under n_count = 60

theory_sd_AIPW_trt <- sqrt(theory_var_AIPW_trt)
theory_sd_AIPW_ctrl <- sqrt(theory_var_AIPW_ctrl)
theory_sd_WR <- sqrt(theory_var_WR)
theory_sd_WD = sqrt(theory_var_WD)

WP_trt_AIPW_coverage <- mean((WP_trt_sim_true > (WP_trt_sim_AIPW - qnorm(0.975)*theory_sd_AIPW_trt)) & 
                               (WP_trt_sim_true < (WP_trt_sim_AIPW + qnorm(0.975)*theory_sd_AIPW_trt)))
WP_ctrl_AIPW_coverage <- mean((WP_ctrl_sim_true > (WP_ctrl_sim_AIPW - qnorm(0.975)*theory_sd_AIPW_ctrl)) & 
                                (WP_ctrl_sim_true < (WP_ctrl_sim_AIPW + qnorm(0.975)*theory_sd_AIPW_ctrl)))
WR_coverage <- mean((WR_sim_true > (WR_sim_AIPW - qnorm(0.975)*theory_sd_WR)) & 
                      (WR_sim_true < (WR_sim_AIPW + qnorm(0.975)*theory_sd_WR)))
WD_coverage <- mean((WD_sim_true > (WD_sim_AIPW - qnorm(0.975)*theory_sd_WD)) & 
                      (WD_sim_true < (WD_sim_AIPW + qnorm(0.975)*theory_sd_WD)))

all_est = cbind(WP_trt_sim_AIPW,WP_ctrl_sim_AIPW,WR_sim_AIPW,WD_sim_AIPW)
wp_true = c(WP_trt_sim_true,WP_ctrl_sim_true,WR_sim_true,WD_sim_true)
wp_est = apply(all_est,2,mean)
wp_sd = apply(all_est,2,sd)
wp_se = sqrt(apply(cbind(theory_var_AIPW_trt,theory_var_AIPW_ctrl,theory_var_WR,theory_var_WD),2,mean))
wp_cr = c(WP_trt_AIPW_coverage,WP_ctrl_AIPW_coverage,WR_coverage,WD_coverage)

res = data.frame(wp_true,wp_est,wp_sd,wp_se,wp_cr)
res
write.csv(res,file="C:/Users/82655/Dropbox/research/covariate adjustment win ratio/AIPW_unb_200.csv")  


