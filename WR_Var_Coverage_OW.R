rm(list=ls())
library(MASS)
library(nnet)
library(ggplot2)
library(parallel)
library(foreach)
library(doParallel)
#numCores <- detectCores()

#############################################################
######### Variance Coverage of OW estimator ###########
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

WP_trt_sim_OW <- numeric(sim_num)
WP_ctrl_sim_OW <- numeric(sim_num)
WR_sim_OW <- numeric(sim_num)
WD_sim_OW <- numeric(sim_num)

theory_var_OW_trt <- numeric(sim_num)
theory_var_OW_ctrl <- numeric(sim_num)
theory_var_WR <- numeric(sim_num)
theory_var_WD <- numeric(sim_num)
Cov_trt_ctrl = numeric(sim_num)

#for balanced design
WP_trt_sim_true <- 0.2953701
WP_ctrl_sim_true <- 0.2050298
WR_sim_true <- 1.447066
WD_sim_true <- 0.0903403
#WP_trt_sim_true <-  0.2959808
#WP_ctrl_sim_true <- 0.2048627
#WR_sim_true <- 1.45124
#WD_sim_true <- 0.0911181

n_count <- 200

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
  treatment_assignment <- rbinom(n_count, 1, 0.5) #0.5 for balanced; 0.7 for unbalanced
  
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
  
  ##########################OW########################################
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
  
  # Calculate propensity scores; including intercept
  PropScore <- glm(treatment~x1 + x2 + x3 + x4 + x5 + x6, data = df_comb, family=binomial)
  pi_func <- fitted(PropScore)
  
  compare_rows_OW <- function(i) {
    g_hfunc_OW1_vec <- numeric(n_count)
    g_hfunc_OW2_vec <- numeric(n_count)
    
    for (j in 1:(n_count)) {
      if (df_comb[i,1] > df_comb[j,1]) {
        g_hfunc_OW1_vec[j] = (df_comb[i,"treatment"] * (1-df_comb[j,"treatment"]) * (1-pi_func[i]) * pi_func[j])
        g_hfunc_OW2_vec[j] = 0
      } else if (df_comb[i,1] < df_comb[j,1]) {
        g_hfunc_OW1_vec[j] = 0
        g_hfunc_OW2_vec[j] = (df_comb[i,"treatment"] * (1-df_comb[j,"treatment"]) * (1-pi_func[i]) * pi_func[j])
      }
    }
    return(list(OW1 = g_hfunc_OW1_vec, OW2 = g_hfunc_OW2_vec))
  }
  # Parallelize the outer loop
  results_OW = apply(matrix(1:n_count,nrow=1),2,compare_rows_OW)
  
  # Aggregate results
  g_hfunc_OW1 <- do.call(rbind, lapply(results_OW, function(x) x$OW1))
  g_hfunc_OW2 <- do.call(rbind, lapply(results_OW, function(x) x$OW2))
  
  # Calculate tau values
  tau1_OW <- sum(as.vector(g_hfunc_OW1), na.rm = TRUE) /
    (sum(df_comb[,"treatment"]*(1-pi_func), na.rm = TRUE) * sum((1-df_comb[,"treatment"])*pi_func, na.rm = TRUE))
  tau2_OW <- sum(as.vector(g_hfunc_OW2), na.rm = TRUE) /
    (sum(df_comb[,"treatment"]*(1-pi_func), na.rm = TRUE) * sum((1-df_comb[,"treatment"])*pi_func, na.rm = TRUE))
  
  # Update WD_sim_OW and WR_sim_OW
  WP_trt_sim_OW[count_temp] <- tau1_OW
  WP_ctrl_sim_OW[count_temp] <- tau2_OW
  WR_sim_OW[count_temp] <- tau1_OW / tau2_OW
  WD_sim_OW[count_temp] <- tau1_OW - tau2_OW
  
  ###########OW Theoretical Variance###########
  ###g1_OW_trt###
  n <- nrow(df_comb)
  #adjusted_sum_i_neq_j <- sum(sapply(1:n, function(j) {sum(df_comb$treatment[-j]*(1-pi_func[-j])) * 
  #    (1 - df_comb$treatment[j])*pi_func[j]}))
  adjusted_sum_i_neq_j = sum(df_comb[,"treatment"]*(1-pi_func))*sum((1-df_comb[,"treatment"])*pi_func)
  g1_constant_denominator_OW <- (1/(n*(n-1))) * adjusted_sum_i_neq_j
  
  outcomes = df_comb$outcomes_comb
  treatment = df_comb$treatment
  
  g1_OW_trt_fun_parallel <- function(i) {
    pairwise_comparisons_trt <- (1/2) * (
      (((1 - pi_func[i]) * pi_func) * treatment[i] * (1 - treatment) * as.numeric(outcomes[i] > outcomes)) +
        (((1 - pi_func) * pi_func[i]) * treatment * (1 - treatment[i]) * as.numeric(outcomes > outcomes[i]))
    )
    #pairwise_comparisons_trt[i] <- 0  # Avoid self-comparison
    sum_g1_OW_trt <- sum(pairwise_comparisons_trt) / g1_constant_denominator_OW
    return((1/(n-1)) * sum_g1_OW_trt)
  }
  #results_g1_OW_trt <- mclapply(1:n, g1_OW_trt_fun_parallel, df_comb$outcomes_comb, df_comb$treatment, pi_func, mc.cores = numCores - 1)
  results_g1_OW_trt = apply(matrix(1:n,nrow=1),2,g1_OW_trt_fun_parallel)
  
  g1_OW_trt <- unlist(results_g1_OW_trt)
  
  ###g1_OW_ctrl###
  g1_OW_ctrl_fun_parallel <- function(i) {
    pairwise_comparisons_ctrl <- (1/2) * (
      (((1 - pi_func[i]) * pi_func) * treatment[i] * (1 - treatment) * as.numeric(outcomes[i] < outcomes)) +
        (((1 - pi_func) * pi_func[i]) * treatment * (1 - treatment[i]) * as.numeric(outcomes < outcomes[i]))
    )
    #pairwise_comparisons_ctrl[i] <- 0  # Avoid self-comparison
    sum_g1_OW_ctrl <- sum(pairwise_comparisons_ctrl) / g1_constant_denominator_OW
    return((1/(n-1)) * sum_g1_OW_ctrl)
  }
  results_g1_OW_ctrl = apply(matrix(1:n,nrow=1),2, g1_OW_ctrl_fun_parallel)
  g1_OW_ctrl <- unlist(results_g1_OW_ctrl)
  

  ### B ### (Dim: 7*1)
  dBeta_B_OW <- function() {
    n <- nrow(df_comb)
    B <- matrix(0, nrow = 1, ncol = 7) # Initialize B as a zero matrix with 1 row and 7 columns
    # Pre-compute the values that are not dependent on i
    pi = pi_func
    Z = df_comb[, "treatment"]
    X = as.matrix(cbind(1,df_comb[, c("x1", "x2", "x3", "x4", "x5", "x6")]))
    
    # Compute sum_j1 and sum_j2 outside the loop
    sum_j1 <- sum((1 - Z) * pi)
    sum_j2 <- colSums((1 - Z) * pi * (1 - pi) * X)
    for (i in 1:n) {
      pi_i = pi[i]
      Z_i = Z[i]
      X_i = X[i, ]
      term1 <- (-Z_i * pi_i * (1-pi_i)) * X_i * sum_j1
      term2 <- (Z_i * (1 - pi_i)) * sum_j2
      B <- B + term1 + term2
    }
    B <- B / (n * (n - 1))
    return(B)
  }
  dBeta_B_OW_val <- dBeta_B_OW()

  ### C is a numerical
  dBeta_C_OW <- function() {
    n <- nrow(df_comb)
    pi <- pi_func
    Z <- df_comb[, "treatment"]
    sum_i <- sum(Z * (1 - pi))
    sum_j <- sum((1 - Z) * pi)
    C <- (1 / (n * (n - 1))) * sum_i * sum_j
    return(C)
  }
  dBeta_C_OW_val <- dBeta_C_OW()
  
  calculate_An_OW = function(df_comb,pi_func,adjusted_sum_i_neq_j,dBeta_B_OW_val,dBeta_C_OW_val){
    covariate_names <- c("x1", "x2", "x3", "x4", "x5", "x6")
    n <- nrow(df_comb)
    An_first_f = function(i){
      Xi <- c(1,as.matrix(df_comb[i, covariate_names]))
      numerator_vec1 <- numerator_vec2 <- numeric(length(covariate_names)+1)
      for(j in 1:n){
        Xj <- c(1,as.matrix(df_comb[j, covariate_names]))
        e_Xi <- pi_func[i]
        e_Xj <- pi_func[j]
        term1 <- ((1-e_Xj) * Xj - e_Xi * Xi) * (df_comb$treatment[i] * (1 - df_comb$treatment[j])) * 
          (e_Xj * (1 - e_Xi)) * as.numeric(df_comb$outcomes[i] > df_comb$outcomes[j]) 
        term2 <- ((1-e_Xj) * Xj - e_Xi * Xi) * (df_comb$treatment[i] * (1 - df_comb$treatment[j])) * 
          (e_Xj * (1 - e_Xi)) * as.numeric(df_comb$outcomes[i] < df_comb$outcomes[j]) 
        numerator_vec1 <- numerator_vec1 + term1
        numerator_vec2 <- numerator_vec2 + term2
      }
      return(list(numerator1 = numerator_vec1,numerator2 = numerator_vec2))
    }
    An_first_temp <- apply(matrix(1:n,nrow=1),2,An_first_f)
    total_numerator1 <- Reduce("+", lapply(An_first_temp, `[[`, "numerator1"))
    total_numerator2 <- Reduce("+", lapply(An_first_temp, `[[`, "numerator2"))
    An_first_trt = matrix(total_numerator1/adjusted_sum_i_neq_j,nrow=1)
    An_first_ctrl = matrix(total_numerator2/adjusted_sum_i_neq_j,nrow=1)
    An_trt = An_first_trt-dBeta_B_OW_val*(WP_trt_sim_OW[count_temp]/dBeta_C_OW_val)
    An_ctrl = An_first_ctrl-dBeta_B_OW_val*(WP_ctrl_sim_OW[count_temp]/dBeta_C_OW_val)
    return (list(An_trt=An_trt,An_ctrl=An_ctrl))
  }
  An_OW = calculate_An_OW(df_comb,pi_func,adjusted_sum_i_neq_j,dBeta_B_OW_val,dBeta_C_OW_val)
  
  An_trt_OW = matrix(An_OW$An_trt,ncol=1)
  An_ctrl_OW = matrix(An_OW$An_ctrl,ncol=1)
  
  ### E_beta0beta0 ### (Dim: 6*6)
  Ebeta0 <- function() {
    n <- nrow(df_comb)
    X <- as.matrix(cbind(1,df_comb[, c("x1", "x2", "x3", "x4", "x5", "x6")]))
    ps <- pi_func
    Ebeta = crossprod(sqrt(ps*(1-ps)) * X) / n
    return(Ebeta)
    
  }
  
  inv_Ebeta0 = solve(Ebeta0())
  ### l_beta ### (Dim: 6*1)
  l_beta_fun <- function(i) {
    X_i <- as.numeric(c(1,df_comb[i, c("x1", "x2", "x3", "x4", "x5", "x6")]))
    Z_i <- df_comb[i, "treatment"]
    pi_i <- pi_func[i]
    l_beta <- inv_Ebeta0 %*% X_i * (Z_i - pi_i)
    return(l_beta)
  }

  OW_Variance_trt <- function() {
    var_sum <- 0
    An_t <- t(An_trt_OW)
    for (i in 1:n) {
      var_sum <- var_sum + (2 * (g1_OW_trt[i] - WP_trt_sim_OW[count_temp]) + An_t %*% l_beta_fun(i))^2
    }
    OW_Var_trt <- var_sum / n
    return(OW_Var_trt)
  }
  
  OW_Variance_ctrl <- function() {
    var_sum <- 0
    An_t <- t(An_ctrl_OW)
    for (i in 1:n) {
      var_sum <- var_sum + (2 * (g1_OW_ctrl[i] - WP_ctrl_sim_OW[count_temp]) + An_t %*% l_beta_fun(i))^2
    }
    OW_Var_ctrl <- var_sum / n
    return(OW_Var_ctrl)
  }
  
  #covariance 
  OW_cov = function(){
    cov_sum <- 0
    An_t_trt <- t(An_trt_OW)
    An_t_ctrl <- t(An_ctrl_OW)
    for (i in 1:n) {
      l_beta_fun_i = l_beta_fun(i)
      temp_trt = (2 * (g1_OW_trt[i] - WP_trt_sim_OW[count_temp]) + An_t_trt %*% l_beta_fun_i)
      temp_ctrl = (2 * (g1_OW_ctrl[i] - WP_ctrl_sim_OW[count_temp]) + An_t_ctrl %*% l_beta_fun_i)
      cov_sum = cov_sum + temp_trt*temp_ctrl
    }
    OW_cov_results <- cov_sum / n
  }
  
  Var_trt_OW_calculated <- OW_Variance_trt()/nrow(df_comb)
  Var_ctrl_OW_calculated <- OW_Variance_ctrl()/nrow(df_comb)
  cov_OW_calculated <- OW_cov()/nrow(df_comb)
  
  #################OW Estimator Output#########################
  ###Theoretical Variance###
  theory_var_OW_trt[count_temp] <- Var_trt_OW_calculated
  theory_var_OW_ctrl[count_temp] <- Var_ctrl_OW_calculated
  Cov_trt_ctrl[count_temp] <- cov_OW_calculated
  theory_var_WR[count_temp] <- var_ratio(WP_trt_sim_OW[count_temp],WP_ctrl_sim_OW[count_temp],
                                         theory_var_OW_trt[count_temp],theory_var_OW_ctrl[count_temp],
                                         Cov_trt_ctrl[count_temp])
  theory_var_WD[count_temp] <- Var_trt_OW_calculated+Var_ctrl_OW_calculated-2*cov_OW_calculated
}#open on line

end_time <- Sys.time()
run_time <- end_time - start_time
run_time

######Coverage Output######
theory_sd_OW_trt <- sqrt(theory_var_OW_trt)
theory_sd_OW_ctrl <- sqrt(theory_var_OW_ctrl)
theory_sd_WR = sqrt(theory_var_WR)
theory_sd_WD = sqrt(theory_var_WD)

WP_trt_OW_coverage <- mean((WP_trt_sim_true > (WP_trt_sim_OW - qnorm(0.975)*theory_sd_OW_trt)) & 
                              (WP_trt_sim_true < (WP_trt_sim_OW + qnorm(0.975)*theory_sd_OW_trt)))
WP_ctrl_OW_coverage <- mean((WP_ctrl_sim_true > (WP_ctrl_sim_OW - qnorm(0.975)*theory_sd_OW_ctrl)) & 
                               (WP_ctrl_sim_true < (WP_ctrl_sim_OW + qnorm(0.975)*theory_sd_OW_ctrl)))
WR_coverage <- mean((WR_sim_true > (WR_sim_OW - qnorm(0.975)*theory_sd_WR)) & 
                      (WR_sim_true < (WR_sim_OW + qnorm(0.975)*theory_sd_WR)))
WD_coverage <- mean((WD_sim_true > (WD_sim_OW - qnorm(0.975)*theory_sd_WD)) & 
                      (WD_sim_true < (WD_sim_OW + qnorm(0.975)*theory_sd_WD)))

all_est = cbind(WP_trt_sim_OW,WP_ctrl_sim_OW,WR_sim_OW,WD_sim_OW)
wp_true = c(WP_trt_sim_true,WP_ctrl_sim_true,WR_sim_true,WD_sim_true)
wp_est = apply(all_est,2,mean)
wp_sd = apply(all_est,2,sd)
wp_se = sqrt(apply(cbind(theory_var_OW_trt,theory_var_OW_ctrl,theory_var_WR,theory_var_WD),2,mean))
wp_cr = c(WP_trt_OW_coverage,WP_ctrl_OW_coverage,WR_coverage,WD_coverage)
res = data.frame(wp_true,wp_est,wp_sd,wp_se,wp_cr)
res
write.csv(res,file="C:/Users/82655/Dropbox/research/covariate adjustment win ratio/OW_b_200.csv")  

