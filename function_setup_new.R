###unadjusted estimator,unadjusted estimator1 , IPW, OW, AIPW and AOW
###output: WP(win probability), LP(loss probability), WR(win ratio), WD(win difference)
###as well as their SEs and 95% confidence interval

#approximation variance of win ratio
var_ratio = function(uW,uL,varW,varL,cov_WL){
  WR_appro_var = (uW^2/uL^2)*(varW/uW^2-2*cov_WL/(uW*uL)+varL/uL^2)
  return(WR_appro_var)
}

#unadjusted estimation (based on multivariate multi-sample U-statistics in Bebu and M.Lachin (2016))
unadj_estimation = function(df_comb,outcomevar,treatment){
  #input treatment and outcome names in dataset df_comb
  treatments = df_comb[,treatment]
  outcomes = df_comb[,outcomevar]
  outcomes_trt = outcomes[treatments==1]
  outcomes_ctrl = outcomes[treatments==0]
  
  df_trt <- data.frame(outcomes_trt)
  df_ctrl <- data.frame(outcomes_ctrl)
  
  # Initialize matrices
  trtwin <- matrix(NA, nrow = nrow(df_trt), ncol = nrow(df_ctrl))
  ctrlwin <- matrix(NA, nrow = nrow(df_trt), ncol = nrow(df_ctrl))
  
  # Define the comparison function
  compare_rows_unadj <- function(i) {
    trt_row <- df_trt[i, 1]
    trt_vec <- numeric(nrow(df_ctrl))
    ctrl_vec <- numeric(nrow(df_ctrl))
    
    for (j in 1:nrow(df_ctrl)) {
      if (trt_row > df_ctrl[j, 1]) {
        trt_vec[j] = 1
        ctrl_vec[j] = 0
      } else if (trt_row < df_ctrl[j, 1]) {
        trt_vec[j] = 0
        ctrl_vec[j] = 1
      }
    }
    return(list(trt = trt_vec, ctrl = ctrl_vec))
  }
  
  results <- apply(matrix(1:nrow(df_trt),nrow=1),2,compare_rows_unadj)
  
  # Aggregate results
  for (i in 1:length(results)) {
    trtwin[i, ] <- results[[i]]$trt
    ctrlwin[i, ] <- results[[i]]$ctrl
  }
  # Calculate win proportions
  trt_winpr <- sum(as.vector(trtwin), na.rm = TRUE) / (nrow(df_trt) * nrow(df_ctrl))
  ctrl_winpr <- sum(as.vector(ctrlwin), na.rm = TRUE) / (nrow(df_trt) * nrow(df_ctrl))
  
  #win probability (WP_trt), loss probability (WP_ctrl), win ratio and win difference
  WP_trt_unadj <- trt_winpr
  WP_ctrl_unadj <- ctrl_winpr
  WR_unadj <- trt_winpr / ctrl_winpr
  WD_unadj <- trt_winpr-ctrl_winpr
  
  colnames(df_trt) <- "outcomes_comb"
  colnames(df_ctrl) <- "outcomes_comb"
  df_comb <- rbind(df_trt, df_ctrl)
  df_comb$treatment <- c(rep(1,nrow(df_trt)), rep(0,nrow(df_ctrl)))
  
  ###########Theoretical Variance###########
  n1 <- nrow(df_trt)
  n2 <- nrow(df_ctrl)
  sum_treatment <- sum(df_comb[,"treatment"])
  sum_control <- sum(1 - df_comb[,"treatment"])
  
  # Function to parallelize for the treatment group ijik
  p_ijik_trt_parallel <- function(i) {
    temp_inner <- numeric(n2)
    for (j in 1:n2) {
      wij_unadjusted <- 1
      temp_inner[j] <- wij_unadjusted * as.numeric(df_trt[i,"outcomes_comb"] > df_ctrl[j,"outcomes_comb"])
    }
    return((sum(temp_inner))^2)
  }
  
  # Parallel computation of the term1 for the treatment group
  term1_trt_list1 <- apply(matrix(1:(n1),nrow=1),2,p_ijik_trt_parallel)
  term1_trt1 <- sum(unlist(term1_trt_list1))/(n1 * n2 * (n2 - 1))
  term2_trt1 <- (1/(n2-1)) * WP_trt_unadj
  
  p_ijik_trt_results <- term1_trt1 - term2_trt1
  
  ###
  p_ijkj_trt_parallel <- function(j) {
    temp_inner <- numeric(n1)
    for (i in 1:n1) {
      wij_unadjusted <- 1
      temp_inner[i] <- wij_unadjusted * as.numeric(df_trt[i,"outcomes_comb"] > df_ctrl[j,"outcomes_comb"])
    }
    return((sum(temp_inner))^2)
  }
  
  # Parallel computation of the term1 for the treatment group
  term1_trt_list2 <- apply(matrix(1:(n2),nrow=1),2, p_ijkj_trt_parallel)
  term1_trt2 <- sum(unlist(term1_trt_list2))/(n2 * n1 * (n1 - 1))
  term2_trt2 <- (1/(n1-1)) * WP_trt_unadj
  
  p_ijkj_trt_results <- term1_trt2 - term2_trt2
  # Function for the control group
  p_ijik_ctrl_parallel <- function(i) {
    temp_inner <- numeric(n2)
    for (j in 1:n2) {
      wij_unadjusted <- 1
      temp_inner[j] <- wij_unadjusted * as.numeric(df_trt[i,"outcomes_comb"] < df_ctrl[j,"outcomes_comb"])
    }
    return((sum(temp_inner))^2)
  }

  # Parallel computation of the term1 for the control group
  term1_ctrl_list1 <- apply(matrix(1:(n1),nrow=1),2, p_ijik_ctrl_parallel)
  term1_ctrl1 <- sum(unlist(term1_ctrl_list1))/(n1 * n2 * (n2 - 1))
  term2_ctrl1 <- (1/(n2-1)) * WP_ctrl_unadj
  
  p_ijik_ctrl_results <- term1_ctrl1 - term2_ctrl1
  
  p_ijkj_ctrl_parallel <- function(j) {
    temp_inner <- numeric(n1)
    for (i in 1:n1) {
      wij_unadjusted <- 1
      temp_inner[i] <- wij_unadjusted * as.numeric(df_trt[i,"outcomes_comb"] < df_ctrl[j,"outcomes_comb"])
    }
    return((sum(temp_inner))^2)
  }
  
  # Parallel computation of the term1 for the control group
  #term1_ctrl_list <- mclapply(1:nrow(df_comb), p_ijik_ctrl_parallel, df_comb, mc.cores = numCores - 1)
  term1_ctrl_list2 <- apply(matrix(1:(n2),nrow=1),2, p_ijkj_ctrl_parallel)
  term1_ctrl2 <- sum(unlist(term1_ctrl_list2))/(n2 * n1 * (n1 - 1))
  term2_ctrl2 <- (1/(n1-1)) * WP_ctrl_unadj
  
  p_ijkj_ctrl_results <- term1_ctrl2 - term2_ctrl2
  
  #use no weight 
  Var_trt_unadj <- (1/n1)*p_ijik_trt_results + (1/n2)*p_ijkj_trt_results - (n1+n2)/(n1*n2)*(WP_trt_unadj)^2
  Var_ctrl_unadj <- (1/n1)*p_ijik_ctrl_results + (1/n2)*p_ijkj_ctrl_results - (n1+n2)/(n2*n1)*(WP_ctrl_unadj)^2 
  ###the following is to compute covariance
  p_ijik_cov <- function(i) {
    temp_inner1 <- temp_inner2 <- numeric(n2)
    for (j in 1:n2) {
      wij_unadjusted <- 1
      temp_inner1[j] <- wij_unadjusted * as.numeric(df_trt[i,"outcomes_comb"] > df_ctrl[j,"outcomes_comb"])
      temp_inner2[j] <- wij_unadjusted * as.numeric(df_trt[i,"outcomes_comb"] < df_ctrl[j,"outcomes_comb"])
    }
    return((sum(temp_inner1))*(sum(temp_inner2)))
  }
  term1_cov_list1 <- apply(matrix(1:(n1),nrow=1),2, p_ijik_cov)
  p_ijik_cov_results <- sum(unlist(term1_cov_list1))/(n1 * n2 * (n2 - 1))
  
  p_ijkj_cov <- function(j) {
    temp_inner1 <- temp_inner2 <- numeric(n1)
    for (i in 1:n1) {
      wij_unadjusted <- 1
      temp_inner1[i] <- wij_unadjusted * as.numeric(df_trt[i,"outcomes_comb"] > df_ctrl[j,"outcomes_comb"])
      temp_inner2[i] <- wij_unadjusted * as.numeric(df_trt[i,"outcomes_comb"] < df_ctrl[j,"outcomes_comb"])
    }
    return((sum(temp_inner1))*(sum(temp_inner2)))
  }
  term2_cov_list1 <- apply(matrix(1:(n2),nrow=1),2, p_ijkj_cov)
  p_ijkj_cov_results <- sum(unlist(term2_cov_list1))/(n2 * n1 * (n1 - 1))
  cov_unadj <- (1/n1)*p_ijik_cov_results + (1/n2)*p_ijkj_cov_results-(n1+n2)/(n1*n2)*(WP_trt_unadj*WP_ctrl_unadj)
  ###Variance of WR
  Var_WR_unadj <- var_ratio(WP_trt_unadj,WP_ctrl_unadj,Var_trt_unadj,
                            Var_ctrl_unadj,cov_unadj)
  Var_WD_unadj <- Var_trt_unadj+Var_ctrl_unadj-2*cov_unadj
  #################Unadjusted Estimator Output#########################
  Estimation = c(WP_trt_unadj,WP_ctrl_unadj,WR_unadj,WD_unadj)
  SE = c(sqrt(Var_trt_unadj),sqrt(Var_ctrl_unadj),sqrt(Var_WR_unadj),sqrt(Var_WD_unadj))
  conf_low = Estimation-qnorm(0.975)*SE
  conf_high = Estimation+qnorm(0.975)*SE
  res = data.frame(EST=Estimation,SE=SE,conf_low=conf_low,conf_high=conf_high)
  row.names(res)=c("WP","LP","WR","WD")
  return(res)
}


#variance estimation based on derived influence function
unadj_new_estimation = function(df_comb,outcomevar,treatment){
  #input treatment and outcome names in dataset df_comb
  treatments = df_comb[,treatment]
  outcomes = df_comb[,outcomevar]
  n = length(outcomes)
  
  compare_rows_unadj <- function(i) {
    g_hfunc_unadj1_vec <- numeric(n)
    g_hfunc_unadj2_vec <- numeric(n)
    for (j in 1:n) {
      if (outcomes[i] > outcomes[j]) {
        g_hfunc_unadj1_vec[j] = treatments[i] * (1-treatments[j])
        g_hfunc_unadj2_vec[j] = 0
      } else if (outcomes[i] < outcomes[j]) {
        g_hfunc_unadj1_vec[j] = 0
        g_hfunc_unadj2_vec[j] = treatments[i] * (1-treatments[j])
      }
    }
    return(list(unadj1 = g_hfunc_unadj1_vec, unadj2 = g_hfunc_unadj2_vec))
  }
  # Parallelize the outer loop
  results_unadj = apply(matrix(1:n,nrow=1),2,compare_rows_unadj)
  
  # Aggregate results
  g_hfunc_unadj1 <- do.call(rbind, lapply(results_unadj, function(x) x$unadj1))
  g_hfunc_unadj2 <- do.call(rbind, lapply(results_unadj, function(x) x$unadj2))
  
  adjusted_sum_i_neq_j = sum(treatments)*sum(1-treatments)
  g1_constant_denominator_unadj <- (1/(n*(n-1))) * adjusted_sum_i_neq_j
  
  # Calculate tau values
  tau1_unadj <- sum(as.vector(g_hfunc_unadj1), na.rm = TRUE)/adjusted_sum_i_neq_j
  tau2_unadj <- sum(as.vector(g_hfunc_unadj2), na.rm = TRUE)/adjusted_sum_i_neq_j
  
  #point estimate
  WP_trt_unadj <- tau1_unadj
  WP_ctrl_unadj <- tau2_unadj
  WR_unadj <- tau1_unadj / tau2_unadj
  WD_unadj <- tau1_unadj-tau2_unadj
  
  #variance estimation
  g1_unadj_trt_fun <- function(i) {
    # Vectorized operation to calculate pairwise comparisons
    pairwise_comparisons <- (1/2) * (
      (treatments[i] * (1 - treatments) * as.numeric(outcomes[i] > outcomes)) + 
        ((1-treatments[i]) * treatments * as.numeric(outcomes[i] < outcomes))
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
    pairwise_comparisons <- (1/2) * (
      (treatments[i] * (1 - treatments) * as.numeric(outcomes[i] < outcomes))+ 
        ((1-treatments[i]) * treatments * as.numeric(outcomes[i] > outcomes))
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
      var_sum <- var_sum + (2 * (g1_unadj_trt[i] - WP_trt_unadj))^2
    }
    unadj_Var_trt <- var_sum / n
    return(unadj_Var_trt)
  }
  
  unadj_Variance_ctrl <- function() {
    var_sum <- 0
    for (i in 1:n) {
      var_sum <- var_sum + (2 * (g1_unadj_ctrl[i] - WP_ctrl_unadj))^2
    }
    unadj_Var_ctrl <- var_sum / n
    return(unadj_Var_ctrl)
  }
  
  #covariance 
  unadj_cov = function(){
    cov_sum <- 0
    for (i in 1:n) {
      temp_trt = (2 * (g1_unadj_trt[i] - WP_trt_unadj))
      temp_ctrl = (2 * (g1_unadj_ctrl[i] - WP_ctrl_unadj))
      cov_sum = cov_sum + temp_trt*temp_ctrl
    }
    unadj_cov_results <- cov_sum / n
    return(unadj_cov_results)
  }
  
  Var_trt_unadj <- unadj_Variance_trt()/n
  Var_ctrl_unadj <- unadj_Variance_ctrl()/n
  cov_unadj <- unadj_cov()/n

  ###Variance of WR and WD
  Var_WR_unadj <- var_ratio(WP_trt_unadj,WP_ctrl_unadj,Var_trt_unadj,
                            Var_ctrl_unadj,cov_unadj)
  Var_WD_unadj <- Var_trt_unadj+Var_ctrl_unadj-2*cov_unadj
  #################Unadjusted Estimator Output#########################
  Estimation =c(WP_trt_unadj,WP_ctrl_unadj,WR_unadj,WD_unadj)
  SE = c(sqrt(Var_trt_unadj),sqrt(Var_ctrl_unadj),sqrt(Var_WR_unadj),sqrt(Var_WD_unadj))
  conf_low = Estimation-qnorm(0.975)*SE
  conf_high = Estimation+qnorm(0.975)*SE
  res = data.frame(EST=Estimation,SE=SE,conf_low=conf_low,conf_high=conf_high)
  row.names(res)=c("WP","LP","WR","WD")
  return(res)
}


#IPW estimation
ipw_estimation = function(df_comb,outcomevar,treatment,covariate){
  #df_comb is a data.frame including outcome, treatment, covariate in PS model
  n = dim(df_comb)[1]  #sample size
  ps.formula = as.formula(paste(treatment,"~",paste(covariate,collapse="+"),sep=""))
  # Calculate propensity scores
  PropScore <- glm(ps.formula,data = df_comb,family=binomial(link="logit"))
  pi_func <- fitted(PropScore)
  
  #simplify notations--outcomes, treatments, covariates
  outcomes = df_comb[,outcomevar]
  treatments = df_comb[,treatment]
  covariates = model.matrix(ps.formula,data=df_comb) #including intercept as covariate
  p = dim(covariates)[2]
  
  # Define the IPW comparison function
  compare_rows_IPW <- function(i) {
    g_hfunc_IPW1_vec <- numeric(n)
    g_hfunc_IPW2_vec <- numeric(n)
    for (j in 1:n) {
      if (outcomes[i] > outcomes[j]) {
        g_hfunc_IPW1_vec[j] = (treatments[i] * (1-treatments[j]) * 1 / (pi_func[i] * (1-pi_func[j])))
        g_hfunc_IPW2_vec[j] = 0
      } else if (outcomes[i] < outcomes[j]) {
        g_hfunc_IPW1_vec[j] = 0
        g_hfunc_IPW2_vec[j] = (treatments[i] * (1-treatments[j]) * 1 / (pi_func[i] * (1-pi_func[j])))
      }
    }
    return(list(IPW1 = g_hfunc_IPW1_vec, IPW2 = g_hfunc_IPW2_vec))
  }
  # Parallelize the outer loop
  results_IPW = apply(matrix(1:n,nrow=1),2,compare_rows_IPW)
  
  # Aggregate results
  g_hfunc_IPW1 <- do.call(rbind, lapply(results_IPW, function(x) x$IPW1))
  g_hfunc_IPW2 <- do.call(rbind, lapply(results_IPW, function(x) x$IPW2))
  
  # Calculate tau values
  tau1_IPW <- sum(as.vector(g_hfunc_IPW1), na.rm = TRUE) /
    (sum(treatments/pi_func, na.rm = TRUE) * sum((1-treatments)/(1-pi_func)))
  tau2_IPW <- sum(as.vector(g_hfunc_IPW2), na.rm = TRUE) /
    (sum(treatments/pi_func, na.rm = TRUE) * sum((1-treatments)/(1-pi_func)))
  
  # Update WD_sim_IPW and WR_sim_IPW
  WP_trt_ipw <- tau1_IPW
  WP_ctrl_ipw <- tau2_IPW
  WR_ipw <- tau1_IPW / tau2_IPW
  WD_ipw <-tau1_IPW-tau2_IPW
  
  ###########IPW Theoretical Variance###########
  ###g1_IPW_trt###
  # Pre-calculate constants
  adjusted_sum_i_neq_j = sum(treatments/pi_func)*sum((1-treatments)/(1-pi_func))
  g1_constant_denominator_IPW <- (1/(n*(n-1))) * adjusted_sum_i_neq_j
  
  g1_IPW_trt_fun <- function(i) {
    pairwise_comparisons <- (1/2) * (
      ((treatments[i] * (1 - treatments) * as.numeric(outcomes[i] > outcomes)) / 
         (pi_func[i] * (1 - pi_func))) + 
        (((1-treatments[i]) * treatments * as.numeric(outcomes[i] < outcomes)) / 
           (pi_func * (1 - pi_func[i])))
    )
    # Sum and adjust according to the formula
    sum_g1_IPW_trt <- sum(pairwise_comparisons) / g1_constant_denominator_IPW
    return((1/(n-1)) * sum_g1_IPW_trt)
  }
  results_g1_IPW_trt = apply(matrix(1:n,nrow=1),2,g1_IPW_trt_fun)
  # Aggregate results
  g1_IPW_trt <- unlist(results_g1_IPW_trt)
  
  ###g1_IPW_ctrl###
  g1_IPW_ctrl_fun <- function(i) {
    pairwise_comparisons <- (1/2) * (
      ((treatments[i] * (1 - treatments) * as.numeric(outcomes[i] < outcomes)) / 
         (pi_func[i] * (1 - pi_func))) + 
        (((1-treatments[i]) * treatments * as.numeric(outcomes[i] > outcomes)) / 
           (pi_func * (1 - pi_func[i])))
    )
    # Sum and adjust according to the formula
    sum_g1_IPW_ctrl <- sum(pairwise_comparisons) / g1_constant_denominator_IPW
    return((1/(n-1)) * sum_g1_IPW_ctrl)
  }
  results_g1_IPW_ctrl = apply(matrix(1:n,nrow=1),2,g1_IPW_ctrl_fun)
  # Aggregate results
  g1_IPW_ctrl <- unlist(results_g1_IPW_ctrl)
  
  ### B
  dBeta_B <- function() {
    X = as.matrix(covariates)
    # Initialize B as a zero matrix with 1 row and p columns
    B = matrix(0, nrow = 1, ncol = p) 
    # Pre-compute the values that are not dependent on i
    pi = pi_func
    Z = treatments
    # Compute sum_j1 and sum_j2 outside the loop
    sum_j1 <- sum((1 - Z) / (1 - pi))
    sum_j2 <- colSums((1 - Z) * pi * X / (1 - pi))
    for (i in 1:n) {
      pi_i = pi[i]
      Z_i = Z[i]
      X_i = X[i, ]
      term1 <- (-Z_i * (1 - pi_i) / pi_i) * X_i * sum_j1
      term2 <- (Z_i / pi_i) * sum_j2
      B <- B + term1 + term2
    }
    B <- B / (n * (n - 1))
    return(B)
  }
  dBeta_B_val <- dBeta_B()
  
  ### C ### (Dim: 1*1)
  dBeta_C <- function() {
    pi <- pi_func
    Z <- treatments
    sum_i <- sum(Z / pi)
    sum_j <- sum((1 - Z) / (1 - pi))
    C <- (1 / (n * (n - 1))) * sum_i * sum_j
    return(C)
  }
  dBeta_C_val <- dBeta_C()
  
  ###A Treatment### (Dim: 6*1)
  calculate_An_IPW = function(outcomes,treatments,covariates,pi_func,adjusted_sum_i_neq_j,dBeta_B_val,dBeta_C_val){
    An_first_f = function(i){
      Xi <- as.matrix(covariates[i,])
      numerator_vec1 <- numerator_vec2 <- numeric(p)
      for(j in 1:n){
        Xj <- as.matrix(covariates[j,])
        e_Xi <- pi_func[i]
        e_Xj <- pi_func[j]
        term1 <- (e_Xj * Xj - (1 - e_Xi) * Xi) * (treatments[i] * (1 - treatments[j])) *
          as.numeric(outcomes[i] > outcomes[j]) / (e_Xi * (1 - e_Xj))
        term2 <- (e_Xj * Xj - (1 - e_Xi) * Xi) * (treatments[i] * (1 - treatments[j])) *
          as.numeric(outcomes[i] < outcomes[j]) / (e_Xi * (1 - e_Xj))
        numerator_vec1 <- numerator_vec1 + term1
        numerator_vec2 <- numerator_vec2 + term2
      }
      return(list(numerator1 = numerator_vec1,numerator2 = numerator_vec2))
    }
    An_first_temp <- apply(matrix(1:n,nrow=1),2,An_first_f)
    total_numerator1 <- Reduce("+", lapply(An_first_temp, `[[`, "numerator1"))
    total_numerator2 <- Reduce("+", lapply(An_first_temp, `[[`, "numerator2"))
    An_first_trt = total_numerator1/adjusted_sum_i_neq_j
    An_first_ctrl = total_numerator2/adjusted_sum_i_neq_j
    An_trt = An_first_trt-t(dBeta_B_val)*(WP_trt_ipw/dBeta_C_val)
    An_ctrl = An_first_ctrl-t(dBeta_B_val)*(WP_ctrl_ipw/dBeta_C_val)
    return (list(An_trt=An_trt,An_ctrl=An_ctrl))
  }
  An_IPW = calculate_An_IPW(outcomes,treatments,covariates,pi_func,adjusted_sum_i_neq_j,dBeta_B_val,dBeta_C_val)
  
  An_trt_IPW = matrix(An_IPW$An_trt,ncol=1)
  An_ctrl_IPW = matrix(An_IPW$An_ctrl,ncol=1)
  
  ### E_beta0beta0 ### (Dim: 6*6)
  Ebeta0 <- function() {
    X <- as.matrix(covariates)
    ps <- pi_func
    Ebeta = crossprod(sqrt(ps*(1-ps)) * X) / n
    return(Ebeta)
  }
  
  inv_Ebeta0 = solve(Ebeta0())
  ### l_beta ### (Dim: 6*1)
  l_beta_fun <- function(i) {
    X_i <- as.numeric(covariates[i,])
    Z_i <- treatments[i]
    pi_i <- pi_func[i]
    l_beta <- inv_Ebeta0 %*% X_i * (Z_i - pi_i)
    return(l_beta)
  }
  
  IPW_Variance_trt <- function() {
    var_sum <- 0
    An_t <- t(An_trt_IPW)
    for (i in 1:n) {
      var_sum <- var_sum + (2 * (g1_IPW_trt[i] - WP_trt_ipw) + An_t %*% l_beta_fun(i))^2
    }
    IPW_Var_trt <- var_sum / n
    return(IPW_Var_trt)
  }
  
  IPW_Variance_ctrl <- function() {
    var_sum <- 0
    An_t <- t(An_ctrl_IPW)
    for (i in 1:n) {
      var_sum <- var_sum + (2 * (g1_IPW_ctrl[i] - WP_ctrl_ipw) + An_t %*% l_beta_fun(i))^2
    }
    IPW_Var_ctrl <- var_sum / n
    return(IPW_Var_ctrl)
  }
  
  #covariance 
  IPW_cov = function(){
    cov_sum <- 0
    An_t_trt <- t(An_trt_IPW)
    An_t_ctrl <- t(An_ctrl_IPW)
    for (i in 1:n) {
      l_beta_fun_i = l_beta_fun(i)
      temp_trt = (2 * (g1_IPW_trt[i] - WP_trt_ipw) + An_t_trt %*% l_beta_fun_i)
      temp_ctrl = (2 * (g1_IPW_ctrl[i] - WP_ctrl_ipw) + An_t_ctrl %*% l_beta_fun_i)
      cov_sum = cov_sum + temp_trt*temp_ctrl
    }
    IPW_cov_results <- cov_sum / n
  }
  
  Var_trt_ipw <- IPW_Variance_trt()/n
  Var_ctrl_ipw <- IPW_Variance_ctrl()/n
  cov_ipw <- IPW_cov()/n
  
  ###Theoretical Variance of WR and WD
  Var_WR_ipw <- var_ratio(WP_trt_ipw,WP_ctrl_ipw,Var_trt_ipw,Var_ctrl_ipw,cov_ipw)
  Var_WD_ipw <- Var_trt_ipw+Var_ctrl_ipw-2*cov_ipw
  #################IPW Estimator Output#########################
  Estimation=c(WP_trt_ipw,WP_ctrl_ipw,WR_ipw,WD_ipw)
  SE = c(sqrt(Var_trt_ipw),sqrt(Var_ctrl_ipw),sqrt(Var_WR_ipw),sqrt(Var_WD_ipw))
  conf_low = Estimation-qnorm(0.975)*SE
  conf_high = Estimation+qnorm(0.975)*SE
  res = data.frame(EST=Estimation,SE=SE,conf_low=conf_low,conf_high=conf_high)
  row.names(res)=c("WP","LP","WR","WD")
  return(res)
}

#OW estimation
ow_estimation = function(df_comb,outcomevar,treatment,covariate){
  #df_comb is a data.frame including outcome, treatment, covariate in PS model 
  n = dim(df_comb)[1]  #sample size
  ps.formula = as.formula(paste(treatment,"~",paste(covariate,collapse="+"),sep=""))
  # Calculate propensity scores
  PropScore <- glm(ps.formula,data = df_comb,family=binomial(link="logit"))
  pi_func <- fitted(PropScore)
  
  #simplify notations--outcomes, treatments, covariates
  outcomes = df_comb[,outcomevar]
  treatments = df_comb[,treatment]
  covariates = model.matrix(ps.formula,data=df_comb) #including intercept as covariate
  p = dim(covariates)[2]
  
  compare_rows_OW <- function(i) {
    g_hfunc_OW1_vec <- numeric(n)
    g_hfunc_OW2_vec <- numeric(n)
    
    for (j in 1:(n)) {
      if (outcomes[i] > outcomes[j]) {
        g_hfunc_OW1_vec[j] = (treatments[i] * (1-treatments[j]) * (1-pi_func[i]) * pi_func[j])
        g_hfunc_OW2_vec[j] = 0
      } else if (outcomes[i] < outcomes[j]) {
        g_hfunc_OW1_vec[j] = 0
        g_hfunc_OW2_vec[j] = (treatments[i] * (1-treatments[j]) * (1-pi_func[i]) * pi_func[j])
      }
    }
    return(list(OW1 = g_hfunc_OW1_vec, OW2 = g_hfunc_OW2_vec))
  }
  # Parallelize the outer loop
  results_OW = apply(matrix(1:n,nrow=1),2,compare_rows_OW)
  
  # Aggregate results
  g_hfunc_OW1 <- do.call(rbind, lapply(results_OW, function(x) x$OW1))
  g_hfunc_OW2 <- do.call(rbind, lapply(results_OW, function(x) x$OW2))
  
  # Calculate tau values
  tau1_OW <- sum(as.vector(g_hfunc_OW1), na.rm = TRUE) /
    (sum(treatments*(1-pi_func), na.rm = TRUE) * sum((1-treatments)*pi_func, na.rm = TRUE))
  tau2_OW <- sum(as.vector(g_hfunc_OW2), na.rm = TRUE) /
    (sum(treatments*(1-pi_func), na.rm = TRUE) * sum((1-treatments)*pi_func, na.rm = TRUE))
  
  # Update WD_sim_OW and WR_sim_OW
  WP_trt_ow <- tau1_OW
  WP_ctrl_ow <- tau2_OW
  WR_ow <- tau1_OW/tau2_OW
  WD_ow <- tau1_OW-tau2_OW
  
  ###########OW Theoretical Variance###########
  ###g1_OW_trt###
  adjusted_sum_i_neq_j = sum(treatments*(1-pi_func))*sum((1-treatments)*pi_func)
  g1_constant_denominator_OW <- (1/(n*(n-1))) * adjusted_sum_i_neq_j
  
  g1_OW_trt_fun_parallel <- function(i) {
    pairwise_comparisons_trt <- (1/2) * (
      (((1 - pi_func[i]) * pi_func) * treatments[i] * (1 - treatments) * as.numeric(outcomes[i] > outcomes)) +
        (((1 - pi_func) * pi_func[i]) * treatments * (1 - treatments[i]) * as.numeric(outcomes > outcomes[i]))
    )
    sum_g1_OW_trt <- sum(pairwise_comparisons_trt) / g1_constant_denominator_OW
    return((1/(n-1)) * sum_g1_OW_trt)
  }
  results_g1_OW_trt = apply(matrix(1:n,nrow=1),2,g1_OW_trt_fun_parallel)
  g1_OW_trt <- unlist(results_g1_OW_trt)
  
  ###g1_OW_ctrl###
  g1_OW_ctrl_fun_parallel <- function(i) {
    pairwise_comparisons_ctrl <- (1/2) * (
      (((1 - pi_func[i]) * pi_func) * treatments[i] * (1 - treatments) * as.numeric(outcomes[i] < outcomes)) +
        (((1 - pi_func) * pi_func[i]) * treatments * (1 - treatments[i]) * as.numeric(outcomes < outcomes[i]))
    )
    sum_g1_OW_ctrl <- sum(pairwise_comparisons_ctrl) / g1_constant_denominator_OW
    return((1/(n-1)) * sum_g1_OW_ctrl)
  }
  results_g1_OW_ctrl = apply(matrix(1:n,nrow=1),2, g1_OW_ctrl_fun_parallel)
  g1_OW_ctrl <- unlist(results_g1_OW_ctrl)
  
  ###OW A Treatment### (Dim: 6*1)
  dBeta_A_trt_OW <- function(i, j) {
    pi_i = pi_func[i]
    pi_j = pi_func[j]
    Z_i = treatments[i]
    Z_j = treatments[j]
    Y_i = outcomes[i]
    Y_j = outcomes[j]
    X_i = covariates[i,]
    X_j = covariates[j,]
    comparison_i_j = as.numeric(Y_i > Y_j)  # Convert logical TRUE/FALSE to 1/0
    comparison_j_i = as.numeric(Y_j > Y_i)  # Convert logical TRUE/FALSE to 1/0
    A = 0.5 * (
      ((1-pi_j) * X_j - pi_i * X_i) * (1-pi_i) * pi_j * Z_i * (1 - Z_j) * comparison_i_j  +
        ((1-pi_i) * X_i - pi_j * X_j) * (1-pi_j) * pi_i * Z_j * (1 - Z_i) * comparison_j_i
    )
    return(A)
  }
  
  ###OW A Control
  dBeta_A_ctrl_OW <- function(i, j) {
    pi_i = pi_func[i]
    pi_j = pi_func[j]
    Z_i = treatments[i]
    Z_j = treatments[j]
    Y_i = outcomes[i]
    Y_j = outcomes[j]
    X_i = covariates[i,]
    X_j = covariates[j,]
    
    comparison_i_j = as.numeric(Y_i < Y_j)  # Convert logical TRUE/FALSE to 1/0
    comparison_j_i = as.numeric(Y_j < Y_i)  # Convert logical TRUE/FALSE to 1/0
    A = 0.5 * (
      ((1-pi_j) * X_j - pi_i * X_i) * (1-pi_i) * pi_j * Z_i * (1 - Z_j) * comparison_i_j  +
        ((1-pi_i) * X_i - pi_j * X_j) * (1-pi_j) * pi_i * Z_j * (1 - Z_i) * comparison_j_i
    )
    return(A)
  }
  
  ### B ### (Dim: 6*1)
  dBeta_B_OW <- function() {
    X = as.matrix(covariates)
    p = dim(X)[2]
    B <- matrix(0, nrow = 1, ncol = p) # Initialize B as a zero matrix with 1 row and 6 columns
    # Pre-compute the values that are not dependent on i
    pi = pi_func
    Z = treatments
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
  
  ### C ### (Dim: 6*1)
  dBeta_C_OW <- function() {
    pi <- pi_func
    Z <- treatments
    sum_i <- sum(Z * (1 - pi))
    sum_j <- sum((1 - Z) * pi)
    C <- (1 / (n * (n - 1))) * sum_i * sum_j
    return(C)
  }
  dBeta_C_OW_val <- dBeta_C_OW()
  
  ###dBeta_gg_OW Treatment### (Dim: 6*1)
  dBeta_gg_OW_trt <- function(i,j) {
    temp_gg_OW_trt <- (1/2)*
      ((((1-pi_func[i]) * pi_func[j]) * 
          treatments[i]*(1-treatments[j])*
          as.numeric(outcomes[i] > outcomes[j]))
       + (((1-pi_func[j]) * pi_func[i]) * 
            treatments[j]*(1-treatments[i])*
            as.numeric(outcomes[j] > outcomes[]))
      ) /
      ((1/(n*(n-1)))*sum(treatments*(1-pi_func))*sum((1-treatments)*pi_func)
      )
    
    dBeta_gg_OW_trt_output <- dBeta_A_trt_OW(i,j)/dBeta_C_OW_val - temp_gg_OW_trt * dBeta_B_OW_val/dBeta_C_OW_val
    return(dBeta_gg_OW_trt_output)
  }
  
  ###dBeta_gg_OW Control### (Dim: 6*1)
  dBeta_gg_OW_ctrl <- function(i,j) {
    temp_gg_OW_ctrl <- (1/2)*
      ((((1-pi_func[i]) * pi_func[j]) * 
          treatments[i]*(1-treatments[j])*
          as.numeric(outcomes[i] < outcomes[j]))
       + (((1-pi_func[j]) * pi_func[i]) * 
            treatments[j]*(1-treatments[i])*
            as.numeric(outcomes[j] < outcomes[]))
      ) /
      ((1/(n*(n-1)))*sum(treatments*(1-pi_func))*sum((1-treatments)*pi_func)
      )
    dBeta_gg_OW_ctrl_output <- dBeta_A_ctrl_OW(i,j)/dBeta_C_OW_val - temp_gg_OW_ctrl * dBeta_B_OW_val/dBeta_C_OW_val
    return(dBeta_gg_OW_ctrl_output)
  }
  
  calculate_An_OW = function(outcomes,treatments,covariates,pi_func,adjusted_sum_i_neq_j,dBeta_B_OW_val,dBeta_C_OW_val){
    An_first_f = function(i){
      Xi <- as.matrix(covariates[i,])
      p = dim(Xi)[2]
      numerator_vec1 <- numerator_vec2 <- numeric(p)
      for(j in 1:n){
        Xj <- as.matrix(covariates[j,])
        e_Xi <- pi_func[i]
        e_Xj <- pi_func[j]
        term1 <- ((1-e_Xj) * Xj - e_Xi * Xi) * (treatments[i] * (1 - treatments[j])) * 
          (e_Xj * (1 - e_Xi)) * as.numeric(outcomes[i] > outcomes[j]) 
        term2 <- ((1-e_Xj) * Xj - e_Xi * Xi) * (treatments[i] * (1 - treatments[j])) * 
          (e_Xj * (1 - e_Xi)) * as.numeric(outcomes[i] < outcomes[j]) 
        numerator_vec1 <- numerator_vec1 + term1
        numerator_vec2 <- numerator_vec2 + term2
      }
      return(list(numerator1 = numerator_vec1,numerator2 = numerator_vec2))
    }
    An_first_temp <- apply(matrix(1:n,nrow=1),2,An_first_f)
    total_numerator1 <- Reduce("+", lapply(An_first_temp, `[[`, "numerator1"))
    total_numerator2 <- Reduce("+", lapply(An_first_temp, `[[`, "numerator2"))
    An_first_trt = total_numerator1/adjusted_sum_i_neq_j
    An_first_ctrl = total_numerator2/adjusted_sum_i_neq_j
    An_trt = An_first_trt-t(dBeta_B_OW_val)*(WP_trt_ow/dBeta_C_OW_val)
    An_ctrl = An_first_ctrl-t(dBeta_B_OW_val)*(WP_ctrl_ow/dBeta_C_OW_val)
    return (list(An_trt=An_trt,An_ctrl=An_ctrl))
  }
  An_OW = calculate_An_OW(outcomes,treatments,covariates,pi_func,adjusted_sum_i_neq_j,dBeta_B_OW_val,dBeta_C_OW_val)
  
  An_trt_OW = matrix(An_OW$An_trt,ncol=1)
  An_ctrl_OW = matrix(An_OW$An_ctrl,ncol=1)
  
  ### E_beta0beta0 ### (Dim: 6*6)
  Ebeta0 <- function() {
    X <- as.matrix(covariates)
    ps <- pi_func
    Ebeta = crossprod(sqrt(ps*(1-ps)) * X) / n
    return(Ebeta)
    
  }
  
  inv_Ebeta0 = solve(Ebeta0())
  ### l_beta 
  l_beta_fun <- function(i) {
    X_i <- as.numeric(covariates[i,])
    Z_i <- treatments[i]
    pi_i <- pi_func[i]
    l_beta <- inv_Ebeta0 %*% X_i * (Z_i - pi_i)
    return(l_beta)
  }
  
  OW_Variance_trt <- function() {
    var_sum <- 0
    An_t <- t(An_trt_OW)
    for (i in 1:n) {
      var_sum <- var_sum + (2 * (g1_OW_trt[i] - WP_trt_ow) + An_t %*% l_beta_fun(i))^2
    }
    OW_Var_trt <- var_sum / n
    return(OW_Var_trt)
  }
  
  OW_Variance_ctrl <- function() {
    var_sum <- 0
    An_t <- t(An_ctrl_OW)
    for (i in 1:n) {
      var_sum <- var_sum + (2 * (g1_OW_ctrl[i] - WP_ctrl_ow) + An_t %*% l_beta_fun(i))^2
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
      temp_trt = (2 * (g1_OW_trt[i] - WP_trt_ow) + An_t_trt %*% l_beta_fun_i)
      temp_ctrl = (2 * (g1_OW_ctrl[i] - WP_ctrl_ow) + An_t_ctrl %*% l_beta_fun_i)
      cov_sum = cov_sum + temp_trt*temp_ctrl
    }
    OW_cov_results <- cov_sum / n
  }
  
  Var_trt_ow <- OW_Variance_trt()/n
  Var_ctrl_ow <- OW_Variance_ctrl()/n
  cov_ow <- OW_cov()/n
  ###Theoretical Variance of WR and WD
  Var_WR_ow <- var_ratio(WP_trt_ow,WP_ctrl_ow,Var_trt_ow,Var_ctrl_ow,cov_ow)
  Var_WD_ow <- Var_trt_ow+Var_ctrl_ow-2*cov_ow
  #################OW Estimator Output#########################
  Estimation = c(WP_trt_ow,WP_ctrl_ow,WR_ow,WD_ow)
  SE = c(sqrt(Var_trt_ow),sqrt(Var_ctrl_ow),sqrt(Var_WR_ow),sqrt(Var_WD_ow))
  conf_low = Estimation-qnorm(0.975)*SE
  conf_high = Estimation+qnorm(0.975)*SE
  res = data.frame(EST=Estimation,SE=SE,conf_low=conf_low,conf_high=conf_high)
  row.names(res)=c("WP","LP","WR","WD")
  return(res)
}


###AIPW estimation
aipw_estimation = function(df_comb,outcomevar,treatment,covariate,covariate_reg){
  #df_comb is a data.frame including outcome, treatment, covariate used in PS model
  #and covariate_reg used in outcome model
  n = dim(df_comb)[1]  #sample size
  ps.formula = as.formula(paste(treatment,"~",paste(covariate,collapse="+"),sep=""))
  # Calculate propensity scores
  PropScore <- glm(ps.formula,data = df_comb,family=binomial(link="logit"))
  pi_func <- fitted(PropScore)
  
  #simplfy notations--outcomes, treatments, covariates
  outcomes = df_comb[,outcomevar]
  treatments = df_comb[,treatment]
  covariates = model.matrix(ps.formula,data=df_comb)
  covariates_reg = as.matrix(df_comb[,covariate_reg])
  p = dim(covariates)[2]
  po = dim(covariates_reg)[2]
  
  ### mu-function
  #actually, we do not use the start points when fitting polr
  #lp_start <- rep(0, po)
  #th_start <- c(-1, 0, 1)
  #start_values <- c(lp_start, th_start)
  model.formual = as.formula(paste("factor(kiss_outcome)","~",paste(covariate_reg,collapse="+"),sep=""))
  dr_trt <- polr(model.formual,data =df_comb[treatments==1,])
  cond_prob_trt <- predict(dr_trt, newdata = df_comb, type = "probs")
  
  dr_ctrl <- polr(model.formual, data = df_comb[treatments == 0,])
  cond_prob_ctrl <- predict(dr_ctrl, newdata = df_comb, type = "probs")
  
  numCores <- 8
  cl <- makeCluster(numCores - 1)  # Use all cores except one
  registerDoParallel(cl)
  
  #number of categories in ordinal outcome
  L = length(unique(outcomes))
  
  # Define the parallelized function
  compare_rows_AIPW <- function(i, outcomes, treatments, pi_func, cond_prob_trt, cond_prob_ctrl) {
    AIPW_num1_vec <- numeric(n)
    AIPW_denom1_vec <- numeric(n)
    AIPW_num2_vec <- numeric(n)
    AIPW_denom2_vec <- numeric(n)
    AIPW_mu1_vec <- numeric(n)
    AIPW_mu2_vec <- numeric(n)
    AIPW_mu_denom1_vec <- numeric(n)
    AIPW_mu_denom2_vec <- numeric(n)
    
    for(j in 1:n){
      if(i != j){
        mu_ij_trt = 0
        for(l in 2:L){
          temp1 = cond_prob_trt[i, l] * (sum(cond_prob_ctrl[j,1:(l-1)]))
          mu_ij_trt = mu_ij_trt+temp1
        }
        AIPW_num1_vec[j] = (treatments[i] * (1 - treatments[j])) * (1/(pi_func[i] * (1 - pi_func[j]))) *
          ((outcomes[i] > outcomes[j]) - mu_ij_trt)
        AIPW_denom1_vec[j] = (treatments[i] * (1 - treatments[j])) * (1/(pi_func[i] * (1 - pi_func[j])))
        
        #mu_j_i
        mu_ij_ctrl = 0
        for(l in 1:(L-1)){
          temp2 = cond_prob_trt[i, l] * (sum(cond_prob_ctrl[j,(l+1):L]))
          mu_ij_ctrl = mu_ij_ctrl+temp2
        }
        AIPW_num2_vec[j] = (treatments[i] * (1 - treatments[j])) * (1/(pi_func[i] * (1 - pi_func[j]))) *
          ((outcomes[i] < outcomes[j]) - mu_ij_ctrl)
        AIPW_denom2_vec[j] = (treatments[i] * (1 - treatments[j])) * (1/(pi_func[i] * (1 - pi_func[j])))
        
        AIPW_mu1_vec[j] = 1 * mu_ij_trt
        AIPW_mu2_vec[j] = 1 * mu_ij_ctrl
        
        AIPW_mu_denom1_vec[j] <- 1
        AIPW_mu_denom2_vec[j] <- 1
      }
    }
    return(list(AIPW_num1 = AIPW_num1_vec, AIPW_denom1 = AIPW_denom1_vec, 
                AIPW_num2 = AIPW_num2_vec, AIPW_denom2 = AIPW_denom2_vec,
                AIPW_mu1 = AIPW_mu1_vec, AIPW_mu_denom1 = AIPW_mu_denom1_vec,
                AIPW_mu2 = AIPW_mu2_vec, AIPW_mu_denom2 = AIPW_mu_denom2_vec))
    
  }
  results_AIPW <- foreach(i=1:n) %dopar% {
    compare_rows_AIPW(i, outcomes, treatments, pi_func, cond_prob_trt, cond_prob_ctrl)
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
  
  WP_trt_aipw <- tau1_AIPW
  WP_ctrl_aipw <- tau2_AIPW
  WR_aipw <- tau1_AIPW/tau2_AIPW
  WD_aipw <- tau1_AIPW-tau2_AIPW
  
  stopCluster(cl)
  
  ########### AIPW Theoretical Variance ###########
  mu_func_trt <- function(a, b, cond_prob_trt, cond_prob_ctrl) {
    mu_a_b = 0
    for(l in 2:L){
      temp1 = cond_prob_trt[a, l] * (sum(cond_prob_ctrl[b,1:(l-1)]))
      mu_a_b = mu_a_b+temp1
    }
    return(mu_a_b)
  }
  
  mu_func_ctrl <- function(a, b, cond_prob_trt, cond_prob_ctrl) {
    mu_a_b = 0
    for(l in 1:(L-1)){
      temp2 = cond_prob_trt[a, l] * (sum(cond_prob_ctrl[b,(l+1):L]))
      mu_a_b = mu_a_b+temp2
    }
    return(mu_a_b)
  }
  
  ###g1_AIPW_trt###
  adjusted_sum_i_neq_j = sum(treatments/pi_func) * sum((1-treatments)/(1-pi_func))
  g1_constant_denominator_AIPW <- (1/(n*(n-1))) * adjusted_sum_i_neq_j
  
  g1_AIPW_trt_fun <- function(i) {
    pairwise_comparisons_AIPW_1 <- (1/2) * (
      ((treatments[i] * (1 - treatments) / (pi_func[i] * (1 - pi_func))) * 
         (as.numeric(outcomes[i] > outcomes) - sapply(1:n, function(b) mu_func_trt(i, b, cond_prob_trt, cond_prob_ctrl)))) +
        (((1 - treatments[i]) * treatments / (pi_func * (1 - pi_func[i]))) * 
           (as.numeric(outcomes > outcomes[i]) - sapply(1:n, function(a) mu_func_trt(a, i, cond_prob_trt, cond_prob_ctrl))))
    )
    pairwise_comparisons_AIPW_2 <- (1/2) * (sapply(1:n, function(b) mu_func_trt(i, b, cond_prob_trt, cond_prob_ctrl)) +
                                              sapply(1:n, function(a) mu_func_trt(a, i, cond_prob_trt, cond_prob_ctrl)))
    
    sum_g1_AIPW_trt <- (sum(pairwise_comparisons_AIPW_1) / g1_constant_denominator_AIPW) + 
      sum(pairwise_comparisons_AIPW_2)
    return((1/(n-1)) * sum_g1_AIPW_trt)
  }
  
  # Parallelize the computation for g1_AIPW_trt
  results_g1_AIPW_trt <- apply(matrix(1:n,nrow=1),2,g1_AIPW_trt_fun)
  g1_AIPW_trt <- unlist(results_g1_AIPW_trt)
  
  
  ###g1_AIPW_ctrl###
  g1_AIPW_ctrl_fun <- function(i) {
    pairwise_comparisons_AIPW_1 <- (1/2) * (
      ((treatments[i] * (1 - treatments) / (pi_func[i] * (1 - pi_func))) * 
         (as.numeric(outcomes[i] < outcomes) - sapply(1:n, function(b) mu_func_ctrl(i, b, cond_prob_trt, cond_prob_ctrl)))) +
        (((1 - treatments[i]) * treatments / (pi_func * (1 - pi_func[i]))) * 
           (as.numeric(outcomes < outcomes[i]) - sapply(1:n, function(a) mu_func_ctrl(a, i, cond_prob_trt, cond_prob_ctrl))))
    )
    pairwise_comparisons_AIPW_2 <- (1/2) * (sapply(1:n, function(b) mu_func_ctrl(i, b, cond_prob_trt, cond_prob_ctrl)) +
                                              sapply(1:n, function(a) mu_func_ctrl(a, i, cond_prob_trt, cond_prob_ctrl)))
    sum_g1_AIPW_ctrl <- (sum(pairwise_comparisons_AIPW_1) / g1_constant_denominator_AIPW) + 
      sum(pairwise_comparisons_AIPW_2)
    return((1/(n-1)) * sum_g1_AIPW_ctrl)
  }
  results_g1_AIPW_ctrl <- apply(matrix(1:n,nrow=1),2,g1_AIPW_ctrl_fun)
  g1_AIPW_ctrl <- unlist(results_g1_AIPW_ctrl)
  
  ### B ###
  ### B1_trt & B1_ctrl is a p*1 vector
  B1_AIPW <- function(outcomes,treatments,covariates,pi_func,cond_prob_trt,cond_prob_ctrl) {
    B1_value_f = function(j){
      Xj <- matrix(covariates[j,],nrow=1)
      numerator_vec1 <- numerator_vec2 <- numeric(p)
      # Loop through all control indices and calculate contributions to numerator and denominator
      for (i in 1:n) {
        Xi <- matrix(covariates[i,],nrow=1)
        e_Xi <- pi_func[i]
        e_Xj <- pi_func[j]
        mu_ij_trt <- mu_func_trt(i, j, cond_prob_trt, cond_prob_ctrl)
        term1 <- (e_Xj * Xj - (1 - e_Xi) * Xi) * (treatments[i] * (1 - treatments[j])) *
          (as.numeric(outcomes[i] > outcomes[j]) - mu_ij_trt) / (e_Xi * (1 - e_Xj))
        numerator_vec1 <- numerator_vec1 + term1
        mu_ij_ctrl <- mu_func_ctrl(i, j, cond_prob_trt, cond_prob_ctrl)
        term2 <- (e_Xj * Xj - (1 - e_Xi) * Xi) * (treatments[i] * (1 - treatments[j])) *
          (as.numeric(outcomes[i] < outcomes[j]) - mu_ij_ctrl) / (e_Xi * (1 - e_Xj))
        numerator_vec2 <- numerator_vec2 + term2
      }
      return(list(numerator1 = numerator_vec1,numerator2 = numerator_vec2))
    }
    B1_values <- apply(matrix(1:n,nrow=1),2,B1_value_f)
    
    total_numerator1 <- Reduce("+", lapply(B1_values, `[[`, "numerator1"))
    total_numerator2 <- Reduce("+", lapply(B1_values, `[[`, "numerator2"))
    denominator = adjusted_sum_i_neq_j
    B1_trt <- total_numerator1 / denominator
    B1_ctrl <- total_numerator2 / denominator
    return(list(B1_trt = B1_trt,B1_ctrl = B1_ctrl))
  }
  
  ### B2_trt & B2_ctrl is numerical
  B2_AIPW <- function(outcomes,treatments,pi_func, cond_prob_trt, cond_prob_ctrl) {
    B2_value_f = function(j){
      numerator_vec1 <- numerator_vec2 <- 0
      # Loop through all control indices and calculate contributions to numerator and denominator
      for (i in 1:n) {
        e_Xi <- pi_func[i]
        e_Xj <- pi_func[j]
        mu_ij_trt <- mu_func_trt(i, j, cond_prob_trt, cond_prob_ctrl)
        term1 <- (treatments[i] * (1 - treatments[j])) * 
          (as.numeric(outcomes[i] > outcomes[j]) - mu_ij_trt) / (e_Xi * (1 - e_Xj))
        numerator_vec1 <- numerator_vec1 + term1
        mu_ij_ctrl <- mu_func_ctrl(i, j, cond_prob_trt, cond_prob_ctrl)
        term2 <- (treatments[i] * (1 - treatments[j])) * 
          (as.numeric(outcomes[i] < outcomes[j]) - mu_ij_ctrl) / (e_Xi * (1 - e_Xj))
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
  
  ### B3 is a p*1 vector
  B3_AIPW_fun <- function(treatments,covariates,pi_func) {
    B3_term1_f = function(j){
      Xj <- matrix(covariates[j,],nrow=1)
      term1_vec <- numeric(p)
      for (i in 1:n) {
        if (i != j) {
          Xi <- matrix(covariates[i,],nrow=1)
          term1_vec <- term1_vec + (-treatments[i]) * (1-pi_func[i]) * Xi / pi_func[i]
        }
      }
      term1_vec <- term1_vec * (1-treatments[j])/(1-pi_func[j])
      return(list(term1 = term1_vec))
    }
    B3_term1 <- apply(matrix(1:n,nrow=1),2,B3_term1_f)
    
    B3_term2_f = function(j){
      Xj <- matrix(covariates[j,],nrow=1)
      term2_vec <- numeric(p)
      for (i in 1:n) {
        if (i != j) {
          term2_vec <- term2_vec + (treatments[i] / pi_func[i])
        }
      }
      term2_vec <- term2_vec * (1 - treatments[j]) * pi_func[j] * Xj / (1 - pi_func[j])
      return(list(term2 = term2_vec))
    }
    B3_term2 <- apply(matrix(1:n,nrow=1),2,B3_term2_f)
    
    total_term1 <- Reduce("+", lapply(B3_term1, `[[`, "term1"))
    total_term2 <- Reduce("+", lapply(B3_term2, `[[`, "term2"))
    return(total_term1 + total_term2)  
  }
  
  ### B4 is a numerical
  B4_AIPW_fun <- function(treatments, pi_func) {
    B4_term_f = function(j){
      term <- 0
      for (i in 1:n) {
        if (i != j) {
          term <- term + (treatments[i] / pi_func[i])
        }
      }
      term <- term * ((1 - treatments[j]) / (1 - pi_func[j]))
      return(list(B4term = term))
    }
    B4_term <- apply(matrix(1:n,nrow=1),2,B4_term_f)
    return(Reduce("+", lapply(B4_term, `[[`, "B4term")))
  }
  #note: B4 is equal to adjusted_sum_i_neq_j
  
  ### Bn ### 
  B1 = B1_AIPW(outcomes,treatments,covariates, pi_func, cond_prob_trt, cond_prob_ctrl)
  B1_trt = matrix(B1$B1_trt,nrow=1)
  B1_ctrl = matrix(B1$B1_ctrl,nrow=1)
  B2 = B2_AIPW(outcomes,treatments, pi_func, cond_prob_trt, cond_prob_ctrl)
  B2_trt = B2$B2_trt
  B2_ctrl = B2$B2_ctrl
  B3 = B3_AIPW_fun(treatments,covariates,pi_func)
  B4 = adjusted_sum_i_neq_j
  
  Bn_AIPW_trt_fun <- function(){
    Bn_trt <- B1_trt - B3*B2_trt/B4
    return(t(Bn_trt))
  }
  
  Bn_AIPW_ctrl_fun <- function(){
    Bn_ctrl <- B1_ctrl - B3*B2_ctrl/B4
    return(t(Bn_ctrl))
  }
  
  ### E_beta0beta0 is p*p matrix
  Ebeta0 <- function() {
    X <- as.matrix(covariates)
    ps <- pi_func
    Ebeta = crossprod(sqrt(ps*(1-ps)) * X) / n
    return(Ebeta)
  }
  
  ### l_beta
  inv_Ebeta0 = solve(Ebeta0())
  l_beta_fun <- function(i) {
    X_i <- as.matrix(covariates[i,])
    Z_i <- treatments[i]
    pi_i <- pi_func[i]
    l_beta <- inv_Ebeta0 %*% X_i * (Z_i - pi_i)
    return(l_beta)
  }
  
  ### E_theta_theta
  thresholds_trt <- dr_trt$zeta
  thresholds_ctrl <- dr_ctrl$zeta
  
  E_thetatheta_fun <- function(outcomes, treatments, covariates_reg, cond_prob_trt, cond_prob_ctrl, thresholds_trt, thresholds_ctrl) {
    F_matrices <- lapply(1:n, function(i) {
      Xi = as.matrix(covariates_reg[i,])
      thresholds <- if(treatments[i] == 1) thresholds_trt else thresholds_ctrl
      cond_prob <- if(treatments[i] == 1) cond_prob_trt else cond_prob_ctrl
      # Calculate delta values based on the outcome for individual i
      delta <- as.numeric(outcomes[i] == 1:L) #which is a 1*L vector
      if(L==3){
        # Calculate the c values for individual i
        c1 <- cond_prob[i, 1]
        c2 <- cond_prob[i, 1] + cond_prob[i,2]
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
      }else{ #in this case L>3
        k = L-1+po
        FM = matrix(0,k,k)
        #the first row
        c1 = cond_prob[i,1]
        f11 = -sum(delta[1:2])*(c1-c1^2)-delta[2]*exp(thresholds[2]-thresholds[1])*((exp(thresholds[2]-thresholds[1])-1)^(-2))
        f12 = delta[2]*exp(thresholds[2]-thresholds[1])*((exp(thresholds[2]-thresholds[1])-1)^(-2))
        #3:(L-1) columns
        f1lm1 = rep(0,L-3)
        f1l = sum(delta[1:2])*(c1 - c1^2)*t(Xi) #1*po
        FM[1,] = c(f11,f12,f1lm1,f1l)
        #2:(L-2) row
        for(l in 2:(L-2)){
          for(j in 1:k){
            if (j < l) FM[l,j] = FM[j,l]
            else if (j == l){
              cl = sum(cond_prob[i,1:l])
              temp1 = sum(delta[l:(l+1)])*(cl-cl^2)
              temp2 = delta[l]*exp(thresholds[l-1]-thresholds[l])*(1-exp(thresholds[l-1]-thresholds[l]))^(-2)
              temp3 = delta[l+1]*exp(thresholds[l+1]-thresholds[l])*((exp(thresholds[l+1]-thresholds[l])-1)^(-2))
              FM[l,j] = -temp1-temp2-temp3
            }
            else if (j==(l+1)) FM[l,j] = temp3
            else if (j<=(L-1)) FM[l,j] = 0
            else FM[l,L:k] = temp1*t(Xi)
          }
        }
        #the L-1 row
        for(j in 1:(L-2)){
          FM[L-1,j] = FM[j,L-1]
        }
        cLm1 = sum(cond_prob[i,1:(L-1)])
        temp1 = sum(delta[(L-1):L])*(cLm1-cLm1^2)
        temp2 = delta[L-1]*exp(thresholds[L-2]-thresholds[L-1])*(1-exp(thresholds[L-2]-thresholds[L-1]))^(-2)
        FM[L-1,L-1] = -temp1-temp2
        FM[L-1,L:k] = temp1*t(Xi)
        #row--L:k, 1:(L-1) columns
        for(l in L:k){
          for(j in 1:(L-1)){
            FM[l,j] = FM[j,l]
          }
        }
        #row--L:k, L:k columns
        factor_val = 0
        for(j in 1:(L-1)){
          cj = sum(cond_prob[i,1:j])
          temp = sum(delta[j:(j+1)])*(cj-cj^2)
          factor_val = factor_val+temp
        }
        FM[L:k,L:k] = -factor_val* (Xi %*% t(Xi))
        F_matrix = FM
      }
      return(F_matrix)
    })
    # Compute the average of the F matrices across all individuals
    E_thetatheta <- Reduce("+", F_matrices) / n
    return(E_thetatheta)
  }
  E_thetatheta <- E_thetatheta_fun(outcomes, treatments, covariates_reg, cond_prob_trt, cond_prob_ctrl, thresholds_trt, thresholds_ctrl)
  
  ### tilde(\mu)(Xi, Xj) 
  mu_tilde_AIPW_fun <- function(i, j, covariates_reg, cond_prob_trt, cond_prob_ctrl) {
    # Extract Xi and Xj as row vectors
    Xi <- as.matrix(covariates_reg[i,])
    Xj <- as.matrix(covariates_reg[j,])
    if(L==3){
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
    }else{
      mu_tilde_vector_l = matrix(0,L-1,1)
      #the first element
      ci1 = cond_prob_trt[i,1]; ci2 = sum(cond_prob_trt[i,1:2])
      cj1 = cond_prob_ctrl[j,1];
      mu_tilde_vector_l[1] = -(ci1-ci1^2)*cj1+(ci2-ci1)*(cj1-cj1^2)
      #the 2 to (L-2) element
      for(l in 2:(L-2)){
        cil = sum(cond_prob_trt[i,1:l]); cilp1 = sum(cond_prob_trt[i,1:(l+1)])
        cjlm1 = sum(cond_prob_ctrl[j,1:(l-1)]); cjl = sum(cond_prob_ctrl[j,1:l])
        mu_tilde_vector_l[l] = (cil-cil^2)*(cjlm1-cjl)+(cilp1-cil)*(cjl-cjl^2)
      }
      #the (L-1)th element
      ciLm1 = sum(cond_prob_trt[i,1:(L-1)]); cjLm2 = sum(cond_prob_ctrl[j,1:(L-2)]) 
      cjLm1 = sum(cond_prob_ctrl[j,1:(L-1)])
      mu_tilde_vector_l[L-1] = (ciLm1-ciLm1^2)*(cjLm2-cjLm1)+(1-ciLm1)*(cjLm1-cjLm1^2)
      #calculate the final term
      factor1 = factor2 = 0
      for(l in 1:(L-2)){
        cil = sum(cond_prob_trt[i,1:l]); cilp1 = sum(cond_prob_trt[i,1:(l+1)])
        cjl = sum(cond_prob_ctrl[j,1:l])
        temp1 = ((cil-cil^2)-(cilp1-cilp1^2))*cjl
        factor1 = factor1 + temp1
        temp2 = (cilp1-cil)*(cjl-cjl^2)
        factor2 = factor2 + temp2
      }
      ciLm1 = sum(cond_prob_trt[i,1:(L-1)]) 
      cjLm1 = sum(cond_prob_ctrl[j,1:(L-1)])
      factor1 = factor1 + (ciLm1-ciLm1^2)*cjLm1
      factor2 = factor2 + (1-ciLm1)*(cjLm1-cjLm1^2)
      mu_tilde_vector_final = factor1*Xi - factor2*Xj
      mu_tilde_vector = rbind(mu_tilde_vector_l,mu_tilde_vector_final)
    }
    return(mu_tilde_vector)
  }
  
  ### Cn 
  Cn_AIPW_fun <- function(treatments, covariates_reg, pi_func, cond_prob_trt, cond_prob_ctrl) {
    # Calculate the first term of Cn
    Cn_first_term_f = function(j){
      numerator_vec <- 0
      # Loop through all control indices and calculate contributions to numerator and denominator
      for (i in 1:n) {
        e_Xi <- pi_func[i]
        e_Xj <- pi_func[j]
        mu_tilde_ij <- mu_tilde_AIPW_fun(i, j, covariates_reg, cond_prob_trt, cond_prob_ctrl)
        term <- (treatments[i] * (1 - treatments[j])) * (-mu_tilde_ij)/(e_Xi * (1 - e_Xj))
        numerator_vec <- numerator_vec + term
      }
      #denominator <- denominator * (1-df$treatment[j])/(1-pi_func[j])
      return(list(numerator = numerator_vec))
    }
    Cn_first_term <- apply(matrix(1:n,nrow=1),2,Cn_first_term_f)
    
    total_numerator <- Reduce("+", lapply(Cn_first_term, `[[`, "numerator"))
    denominator = adjusted_sum_i_neq_j
    first_term <- total_numerator / denominator
    
    # Second term of Cn
    second_term = matrix(0,L-1+po,1)
    for(i in 1:n){
      for(j in 1:n){
        if (i != j) {
          res = mu_tilde_AIPW_fun(i,j, covariates_reg, cond_prob_trt, cond_prob_ctrl)
        }else{
          res = matrix(0,L-1+po,1)
        }
        second_term = second_term + res
      }
    }
    second_term <- second_term/(n*(n-1))
    Cn_value <- first_term + second_term
    return(Cn_value)
  }
  
  # Use the function to calculate Cn
  Cn_value <- Cn_AIPW_fun(treatments, covariates_reg, pi_func, cond_prob_trt, cond_prob_ctrl)
  
  ### l_tilde_theta ###
  E_thetatheta_inv <- ginv(E_thetatheta)
  
  l_theta_AIPW_fun <- function(i,outcomes,treatments,covariates_reg,cond_prob_trt, cond_prob_ctrl, 
                               thresholds_trt, thresholds_ctrl, E_thetatheta_inv) {
    Xi <- as.matrix(covariates_reg[i,])
    cond_prob <- if(treatments[i] == 1) cond_prob_trt else cond_prob_ctrl
    delta <- as.numeric(outcomes[i] == 1:L)
    thresholds <- if(treatments[i] == 1) thresholds_trt else thresholds_ctrl
    if(L==3){
      c1 <- cond_prob[i, 1]
      c2 <- c1 + cond_prob[i, 2]
      # Construct the l_tilde vector according to the formula
      l_tilde_vector_1 <- delta[1] - (delta[1] + delta[2]) * c1 - delta[2] * (exp(thresholds[2] - thresholds[1])-1)^(-1)
      l_tilde_vector_2 <- -(delta[2] + delta[3]) * c2 + delta[2] * (1 - exp(thresholds[1] - thresholds[2]))^(-1)
      l_tilde_vector_3 <- (-(delta[1] + delta[2]) + (delta[1] + delta[2]) * c1 + (delta[2] + delta[3]) * c2) * Xi
      l_tilde_vector <- rbind(l_tilde_vector_1, l_tilde_vector_2, l_tilde_vector_3)
      
    }else{
      l_tilde_vector_l = matrix(0,L-1,1)
      c1 = cond_prob[i,1]
      l_tilde_vector_l[1] = delta[1]-sum(delta[1:2])*c1-delta[2]*(exp(thresholds[2]-thresholds[1])-1)^(-1)
      for(l in 2:(L-2)){
        cl = sum(cond_prob[i,1:l])
        l_tilde_vector_l[l] = -sum(delta[l:(l+1)])*cl+delta[l]*1/(1-exp(thresholds[l-1]-thresholds[l]))-delta[l+1]*1/(exp(thresholds[l+1]-thresholds[l])-1)
      }
      cLm1 = sum(cond_prob[i,1:(L-1)])
      l_tilde_vector_l[L-1] = -sum(delta[(L-1):L])*cLm1+delta[L-1]*1/(1-exp(thresholds[L-2]-thresholds[L-1]))
      factor_val = 0
      for(l in 1:(L-1)){
        cl = sum(cond_prob[i,1:l])
        temp = sum(delta[l:(l+1)])*cl
        factor_val = factor_val+temp
      }
      factor_val = factor_val-sum(delta[1:(L-1)])
      l_tilde_vector_finall = factor_val*Xi
      l_tilde_vector = rbind(l_tilde_vector_l,l_tilde_vector_finall)
    }
    l_tilde_value <- E_thetatheta_inv %*% l_tilde_vector
    return(l_tilde_value)
  }
  
  AIPW_Variance_trt <- function() {
    var_sum <- 0
    Bn_t <- t(Bn_AIPW_trt_fun())
    Cn_t <- t(Cn_value)
    for (i in 1:n) {
      l_theta <- l_theta_AIPW_fun(i, outcomes,treatments,covariates_reg,cond_prob_trt, cond_prob_ctrl, 
                                  thresholds_trt, thresholds_ctrl, E_thetatheta_inv)
      var_sum <- var_sum+(2*(g1_AIPW_trt[i]-WP_trt_aipw)+Bn_t%*%l_beta_fun(i)+Cn_t%*%l_theta)^2
    }
    AIPW_Var_trt <- var_sum / n
    return(AIPW_Var_trt)
  }
  
  AIPW_Variance_ctrl <- function() {
    var_sum <- 0
    Bn_t <- t(Bn_AIPW_ctrl_fun())
    Cn_t <- t(Cn_value)
    for (i in 1:n) {
      l_theta <- l_theta_AIPW_fun(i,outcomes,treatments,covariates_reg,cond_prob_trt, cond_prob_ctrl, 
                                  thresholds_trt, thresholds_ctrl, E_thetatheta_inv)
      var_sum <- var_sum+(2*(g1_AIPW_ctrl[i]-WP_ctrl_aipw)+Bn_t%*%l_beta_fun(i)+Cn_t%*%l_theta)^2
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
      l_theta <- l_theta_AIPW_fun(i,outcomes,treatments,covariates_reg,cond_prob_trt, cond_prob_ctrl, 
                                  thresholds_trt, thresholds_ctrl, E_thetatheta_inv)
      l_beta <- l_beta_fun(i)
      temp_trt = (2*(g1_AIPW_trt[i] - WP_trt_aipw) + Bn_t_trt %*% l_beta + Cn_t %*% l_theta)
      temp_ctrl = (2*(g1_AIPW_ctrl[i] - WP_ctrl_aipw) + Bn_t_ctrl %*% l_beta + Cn_t %*% l_theta)
      cov_sum <- cov_sum + temp_trt*temp_ctrl
    }
    AIPW_cov_results <- cov_sum / n
    return(AIPW_cov_results)
  }
  
  
  Var_trt_aipw <- AIPW_Variance_trt()/n
  Var_ctrl_aipw <- AIPW_Variance_ctrl()/n
  cov_aipw <- AIPW_cov()/n
  ###Theoretical Variance of WR and WD
  Var_WR_aipw <- var_ratio(WP_trt_aipw,WP_ctrl_aipw,Var_trt_aipw,Var_ctrl_aipw,cov_aipw)
  Var_WD_aipw <- Var_trt_aipw+Var_ctrl_aipw-2*cov_aipw
  #################AIPW Estimator Output#########################
  Estimation = c(WP_trt_aipw,WP_ctrl_aipw,WR_aipw,WD_aipw)
  SE = c(sqrt(Var_trt_aipw),sqrt(Var_ctrl_aipw),sqrt(Var_WR_aipw),sqrt(Var_WD_aipw))
  conf_low = Estimation-qnorm(0.975)*SE
  conf_high = Estimation+qnorm(0.975)*SE
  res = data.frame(EST=Estimation,SE=SE,conf_low=conf_low,conf_high=conf_high)
  row.names(res)=c("WP","LP","WR","WD")
  return(res)
}

###AOW estimation
aow_estimation = function(df_comb,outcomevar,treatment,covariate,covariate_reg){
  #df_comb is a data.frame including outcome, treatment, covariate in PS model
  #and covariate_reg in outcome model
  n = dim(df_comb)[1]  #sample size
  ps.formula = as.formula(paste(treatment,"~",paste(covariate,collapse="+"),sep=""))
  # Calculate propensity scores
  PropScore <- glm(ps.formula,data = df_comb,family=binomial(link="logit"))
  pi_func <- fitted(PropScore)
  
  #simplfy notations--outcomes, treatments, covariates
  outcomes = df_comb[,outcomevar]
  treatments = df_comb[,treatment]
  covariates = model.matrix(ps.formula,data=df_comb) #including intercept as covariate
  covariates_reg = as.matrix(df_comb[,covariate_reg])
  p = dim(covariates)[2]
  po = dim(covariates_reg)[2]
  
  ### mu-function
  #actually, we do not use the start points when fitting polr
  #lp_start <- rep(0, po)
  #th_start <- c(-1, 0, 1)
  #start_values <- c(lp_start, th_start)
  model.formual = as.formula(paste("factor(kiss_outcome)","~",paste(covariate_reg,collapse="+"),sep=""))
  dr_trt <- polr(model.formual,data =df_comb[treatments==1,])
  cond_prob_trt <- predict(dr_trt, newdata = df_comb, type = "probs")
  
  dr_ctrl <- polr(model.formual, data = df_comb[treatments == 0,])
  cond_prob_ctrl <- predict(dr_ctrl, newdata = df_comb, type = "probs")
  
  numCores <- 8
  cl <- makeCluster(numCores - 1)  # Use all cores except one
  registerDoParallel(cl)
  
  #number of categories in ordinal outcome
  L = length(unique(outcomes))
  
  compare_rows_AOW <- function(i,outcomes,treatments,pi_func,cond_prob_trt,cond_prob_ctrl) {
    AOW_num1_vec <- numeric(n)
    AOW_denom1_vec <- numeric(n)
    AOW_num2_vec <- numeric(n)
    AOW_denom2_vec <- numeric(n)
    AOW_mu1_vec <- numeric(n)
    AOW_mu2_vec <- numeric(n)
    AOW_mu_denom1_vec <- numeric(n)
    AOW_mu_denom2_vec <- numeric(n)
    
    for(j in 1:n){
      if(i != j){
        AOW_num1_vec[j] = (treatments[i] * (1 - treatments[j])) * ((1-pi_func[i]) * pi_func[j]) *
          ((outcomes[i] > outcomes[j]) - (cond_prob_trt[i, 2] * cond_prob_ctrl[j, 1] + 
                                            cond_prob_trt[i, 3] * (cond_prob_ctrl[j, 1] + cond_prob_ctrl[j, 2])))
        AOW_denom1_vec[j] = (treatments[i] * (1 - treatments[j])) * ((1-pi_func[i]) * pi_func[j])
        
        AOW_num2_vec[j] = (treatments[i] * (1 - treatments[j])) * ((1-pi_func[i]) * pi_func[j]) *
          ((outcomes[i] < outcomes[j]) - (cond_prob_trt[i, 1] * cond_prob_ctrl[j, 2] + 
                                            cond_prob_trt[i, 1] * cond_prob_ctrl[j, 3] +
                                            cond_prob_trt[i, 2] * cond_prob_ctrl[j, 3]))
        AOW_denom2_vec[j] = (treatments[i] * (1 - treatments[j])) * ((1-pi_func[i]) * pi_func[j])
        
        AOW_mu1_vec[j] = (pi_func[i] * pi_func[j] * (1-pi_func[i]) *(1-pi_func[j])) * 
          (cond_prob_trt[i, 2] * cond_prob_ctrl[j, 1] + 
             cond_prob_trt[i, 3] * cond_prob_ctrl[j, 1] + 
             cond_prob_trt[i, 3] * cond_prob_ctrl[j, 2])
        
        AOW_mu2_vec[j] = (pi_func[i] * pi_func[j] * (1-pi_func[i]) *(1-pi_func[j])) * 
          (cond_prob_trt[i, 1] * cond_prob_ctrl[j, 2] + 
             cond_prob_trt[i, 1] * cond_prob_ctrl[j, 3] + 
             cond_prob_trt[i, 2] * cond_prob_ctrl[j, 3])
        
        AOW_mu_denom1_vec[j] = pi_func[i] * pi_func[j] * (1-pi_func[i]) *(1-pi_func[j])
        AOW_mu_denom2_vec[j] = pi_func[i] * pi_func[j] * (1-pi_func[i]) *(1-pi_func[j])
      }
    }
    
    return(list(AOW_num1 = AOW_num1_vec, AOW_denom1 = AOW_denom1_vec, 
                AOW_num2 = AOW_num2_vec, AOW_denom2 = AOW_denom2_vec, 
                AOW_mu1 = AOW_mu1_vec, AOW_mu_denom1 = AOW_mu_denom1_vec,
                AOW_mu2 = AOW_mu2_vec, AOW_mu_denom2 = AOW_mu_denom2_vec))
  }
  
  results_AOW <- foreach(i=1:n) %dopar% {
    compare_rows_AOW(i,outcomes,treatments,pi_func, cond_prob_trt, cond_prob_ctrl)
  }
  
  # Combine the results
  tau_AOW_num1 <- do.call(rbind, lapply(results_AOW, function(x) x$AOW_num1))
  tau_AOW_denom1 <- do.call(rbind, lapply(results_AOW, function(x) x$AOW_denom1))
  tau_AOW_num2 <- do.call(rbind, lapply(results_AOW, function(x) x$AOW_num2))
  tau_AOW_denom2 <- do.call(rbind, lapply(results_AOW, function(x) x$AOW_denom2))
  tau_AOW_mu1 <- do.call(rbind, lapply(results_AOW, function(x) x$AOW_mu1))
  tau_AOW_mu_denom1 <- do.call(rbind, lapply(results_AOW, function(x) x$AOW_mu_denom1))
  tau_AOW_mu2 <- do.call(rbind, lapply(results_AOW, function(x) x$AOW_mu2))
  tau_AOW_mu_denom2 <- do.call(rbind, lapply(results_AOW, function(x) x$AOW_mu_denom2))
  
  tau1_AOW <- (sum(as.vector(tau_AOW_num1), na.rm = T) / sum(as.vector(tau_AOW_denom1), na.rm = T)) +
    (sum(as.vector(tau_AOW_mu1), na.rm = T)/sum(as.vector(tau_AOW_mu_denom1), na.rm = T))
  
  tau2_AOW <- (sum(as.vector(tau_AOW_num2), na.rm = T) / sum(as.vector(tau_AOW_denom2), na.rm = T)) +
    (sum(as.vector(tau_AOW_mu2), na.rm = T)/sum(as.vector(tau_AOW_mu_denom2), na.rm = T))
  
  #point estimation
  WP_trt_aow <- tau1_AOW
  WP_ctrl_aow <- tau2_AOW
  WR_aow <- tau1_AOW / tau2_AOW
  WD_aow <- tau1_AOW-tau2_AOW
  stopCluster(cl)
  
  ########### AOW Theoretical Variance###########
  mu_func_trt <- function(a, b, cond_prob_trt, cond_prob_ctrl) {
    mu_a_b = 0
    for(l in 2:L){
      temp1 = cond_prob_trt[a, l] * (sum(cond_prob_ctrl[b,1:(l-1)]))
      mu_a_b = mu_a_b+temp1
    }
    return(mu_a_b)
  }
  
  mu_func_ctrl <- function(a, b, cond_prob_trt, cond_prob_ctrl) {
    mu_a_b = 0
    for(l in 1:(L-1)){
      temp2 = cond_prob_trt[a, l] * (sum(cond_prob_ctrl[b,(l+1):L]))
      mu_a_b = mu_a_b+temp2
    }
    return(mu_a_b)
  }
  
  
  ###g1_AOW_trt###
  adjusted_denom1_term = sum(treatments*(1-pi_func))*sum((1-treatments)*pi_func)
  g1_denominator1_AOW <- (1/(n*(n-1))) * adjusted_denom1_term
  adjusted_denom2_term <- sum(sapply(1:n, function(j) {sum(pi_func[-j]*(1-pi_func[-j])) * 
      (pi_func[j]) * (1-pi_func[j])}))
  g1_denominator2_AOW <- (1/(n*(n-1))) * adjusted_denom2_term
  
  
  g1_AOW_trt_fun <- function(i) {
    pairwise_comparisons_AOW_term1 <- (1/2) * (
      ((treatments[i] * (1 - treatments) * (1 - pi_func[i]) * pi_func) * 
         (as.numeric(outcomes[i] > outcomes) - sapply(1:n, function(b) mu_func_trt(i, b, cond_prob_trt, cond_prob_ctrl)))) +
        (((1 - treatments[i]) * treatments * pi_func[i] * (1 - pi_func)) * 
           (as.numeric(outcomes > outcomes[i]) - sapply(1:n, function(a) mu_func_trt(a, i, cond_prob_trt, cond_prob_ctrl))))
    )
    pairwise_comparisons_AOW_term2 <- (1/2) * pi_func[i] * (1 - pi_func[i]) *
      pi_func * (1 - pi_func) * (sapply(1:n, function(b) mu_func_trt(i, b, cond_prob_trt, cond_prob_ctrl)) +
                                   sapply(1:n, function(a) mu_func_trt(a, i, cond_prob_trt, cond_prob_ctrl)))
    # Sum and adjust according to the formula
    sum_g1_AOW_trt <- (sum(pairwise_comparisons_AOW_term1) / g1_denominator1_AOW) +
      (sum(pairwise_comparisons_AOW_term2) / g1_denominator2_AOW)
    return((1/(n-1)) * sum_g1_AOW_trt)
    
  }
  
  # Parallelize the computation for g1_AIPW_trt
  results_g1_AOW_trt = apply(matrix(1:n,nrow=1),2,g1_AOW_trt_fun)
  g1_AOW_trt <- unlist(results_g1_AOW_trt)
  
  g1_AOW_ctrl_fun <- function(i) {
    pairwise_comparisons_AOW_term1 <- (1/2) * (
      ((treatments[i] * (1 - treatments) * (1 - pi_func[i]) * pi_func) * 
         (as.numeric(outcomes[i] < outcomes) - sapply(1:n, function(b) mu_func_ctrl(i, b, cond_prob_trt, cond_prob_ctrl)))) +
        (((1 - treatments[i]) * treatments * pi_func[i] * (1 - pi_func)) * 
           (as.numeric(outcomes < outcomes[i]) - sapply(1:n, function(a) mu_func_ctrl(a, i, cond_prob_trt, cond_prob_ctrl))))
    )
    pairwise_comparisons_AOW_term2 <- (1/2) * pi_func[i] * (1 - pi_func[i]) *
      pi_func * (1 - pi_func) * (sapply(1:n, function(b) mu_func_ctrl(i, b, cond_prob_trt, cond_prob_ctrl)) +
                                   sapply(1:n, function(a) mu_func_ctrl(a, i, cond_prob_trt, cond_prob_ctrl)))
    
    sum_g1_AOW_ctrl <- (sum(pairwise_comparisons_AOW_term1) / g1_denominator1_AOW) +
      (sum(pairwise_comparisons_AOW_term2) / g1_denominator2_AOW)
    return((1/(n-1)) * sum_g1_AOW_ctrl)
    
  }
  results_g1_AOW_ctrl = apply(matrix(1:n,nrow=1),2,g1_AOW_ctrl_fun)
  g1_AOW_ctrl <- unlist(results_g1_AOW_ctrl)
  
  ### B ###
  ### B1_trt & B1_ctrl is a p*1 matrix
  B1_AOW <- function(outcomes,treatments,covariates,pi_func, cond_prob_trt, cond_prob_ctrl) {
    B1_value_f = function(j){
      Xj <- as.matrix(covariates[j,])
      numerator_vec1 <- numerator_vec2 <- numeric(p)
      # Loop through all control indices and calculate contributions to numerator and denominator
      for (i in 1:n) {
        Xi <- as.matrix(covariates[i,])
        e_Xi <- pi_func[i]
        e_Xj <- pi_func[j]
        mu_ij_trt <- mu_func_trt(i, j, cond_prob_trt, cond_prob_ctrl)
        term1 <- (treatments[i] * (1 - treatments[j])) * (1 - e_Xi) * e_Xj * 
          ((1-e_Xj) * Xj - e_Xi * Xi) * (as.numeric(outcomes[i] > outcomes[j]) - mu_ij_trt)
        numerator_vec1 <- numerator_vec1 + term1
        mu_ij_ctrl <- mu_func_ctrl(i, j, cond_prob_trt, cond_prob_ctrl)
        term2 <- (treatments[i] * (1 - treatments[j])) * (1 - e_Xi) * e_Xj * 
          ((1-e_Xj) * Xj - e_Xi * Xi) * (as.numeric(outcomes[i] < outcomes[j]) - mu_ij_ctrl)
        numerator_vec2 <- numerator_vec2 + term2
      }
      return(list(numerator1 = numerator_vec1, numerator2 = numerator_vec2))
    }
    B1_values <- apply(matrix(1:n,nrow=1),2,B1_value_f)
    
    total_numerator1 <- Reduce("+", lapply(B1_values, `[[`, "numerator1"))
    total_numerator2 <- Reduce("+", lapply(B1_values, `[[`, "numerator2"))
    denominator <- adjusted_denom1_term 
    B1_trt <- total_numerator1 / denominator
    B1_ctrl <- total_numerator2 / denominator
    return(list(B1_trt = B1_trt,B1_ctrl = B1_ctrl))
  }
  B1 = B1_AOW(outcomes,treatments,covariates,pi_func, cond_prob_trt, cond_prob_ctrl) 
  
  ### B2_trt & B2_ctrl 
  B2_AOW <- function(outcomes,treatments,pi_func, cond_prob_trt, cond_prob_ctrl) {
    B2_value_f = function(j){
      numerator_vec1 <- numerator_vec2 <- 0
      # Loop through all control indices and calculate contributions to numerator and denominator
      for (i in 1:n) {
        e_Xi <- pi_func[i]
        e_Xj <- pi_func[j]
        mu_ij_trt <- mu_func_trt(i, j, cond_prob_trt, cond_prob_ctrl)
        term1 <- (treatments[i] * (1 - treatments[j])) * (1-e_Xi) * e_Xj *
          (as.numeric(outcomes[i] > outcomes[j]) - mu_ij_trt)
        numerator_vec1 <- numerator_vec1 + term1
        mu_ij_ctrl <- mu_func_ctrl(i, j, cond_prob_trt, cond_prob_ctrl)
        term2 <- (treatments[i] * (1 - treatments[j])) * (1-e_Xi) * e_Xj *
          (as.numeric(outcomes[i] < outcomes[j]) - mu_ij_ctrl)
        numerator_vec2 <- numerator_vec2 + term2
      }
      return(list(numerator1 = numerator_vec1, numerator2 = numerator_vec2))
    }
    
    B2_values <- apply(matrix(1:n,nrow=1),2,B2_value_f )
    total_numerator1 <- Reduce("+", lapply(B2_values, `[[`, "numerator1"))
    total_numerator2 <- Reduce("+", lapply(B2_values, `[[`, "numerator2"))
    denominator <- adjusted_denom1_term 
    B2_trt <- total_numerator1 / denominator
    B2_ctrl <- total_numerator2 / denominator
    return(list(B2_trt = B2_trt,B2_ctrl = B2_ctrl))
  }
  B2 = B2_AOW(outcomes,treatments,pi_func, cond_prob_trt, cond_prob_ctrl) 
  
  ### B3 ###
  B3_AOW_fun <- function(treatments,covariates,pi_func) {
    B3_value_f = function(j){
      Xj <- as.matrix(covariates[j,])
      numerator_vec <- numeric(p)
      # Loop through all control indices and calculate contributions to numerator and denominator
      for (i in 1:n) {
        Xi <- as.matrix(covariates[i,])
        e_Xi <- pi_func[i]
        e_Xj <- pi_func[j]
        term <- (treatments[i] * (1 - treatments[j])) * (1 - e_Xi) * e_Xj * 
          ((1-e_Xj) * Xj - e_Xi * Xi)
        numerator_vec <- numerator_vec + term
      }
      return(list(numerator = numerator_vec))
    }
    B3_values <- apply(matrix(1:n,nrow=1),2,B3_value_f)
    
    total_numerator <- Reduce("+", lapply(B3_values, `[[`, "numerator"))
    return(total_numerator)
  }
  B3 = B3_AOW_fun(treatments,covariates,pi_func)
  
  ### B4 is a numeric
  B4_AOW_fun <- function(treatments,pi_func) {
    B4_term_f = function(j){
      term_val <- 0
      for (i in 1:n) {
        if (i != j) {
          term_val <- term_val + treatments[i] * (1-pi_func[i])
        }
      }
      term_val <- term_val * (1 - treatments[j]) * pi_func[j]
      return(list(B4_term_val = term_val))
    }
    B4_term <- apply(matrix(1:n,nrow=1),2,B4_term_f)
    B4_AOW <- Reduce("+", lapply(B4_term, `[[`, "B4_term_val"))
    return(B4_AOW)  
  }
  ##note: B4=adjusted_denom1_term
  B4 = adjusted_denom1_term
  
  ### B8 ### (Dim: 1*1), since B8 is the denominator of B5 and B6
  B8_AOW_fun <- function(pi_func) {
    B8_term_f = function(j){
      term_val <- 0
      for (i in 1:n) {
        if (i != j) {
          term_val <- term_val + pi_func[i] * (1-pi_func[i])
        }
      }
      term_val <- term_val * pi_func[j] * (1-pi_func[j])
      return(list(B8_term_val = term_val))
    }
    B8_term <- apply(matrix(1:n,nrow=1),2,B8_term_f)
    B8_AOW <- Reduce("+", lapply(B8_term, `[[`, "B8_term_val"))
    return(B8_AOW)  
  }
  B8 = B8_AOW_fun(pi_func)
  
  ### B5 ### (Dim: 6*1)
  B5_AOW <- function(covariates, pi_func, cond_prob_trt, cond_prob_ctrl) {
    B5_value_f = function(j){
      numerator_vec1 <- numerator_vec2 <- numeric(p)
      # Loop through all control indices and calculate contributions to numerator and denominator
      for (i in 1:n) {
        if (i != j) {
          Xi <- as.matrix(covariates[i,])
          e_Xi <- pi_func[i]
          e_Xj <- pi_func[j]
          h_ij <- e_Xi * (1-e_Xi) * e_Xj * (1-e_Xj)
          mu_ij_trt <- mu_func_trt(i, j, cond_prob_trt, cond_prob_ctrl)
          term1 <- h_ij * ((1 - 2 * e_Xi) * Xi) * mu_ij_trt
          numerator_vec1 <- numerator_vec1 + term1
          mu_ij_ctrl <- mu_func_ctrl(i, j, cond_prob_trt, cond_prob_ctrl)
          term2 <- h_ij * ((1 - 2 * e_Xi) * Xi) * mu_ij_ctrl
          numerator_vec2 <- numerator_vec2 + term2
        }
      }
      return(list(numerator1 = numerator_vec1, numerator2 = numerator_vec2))
    }
    B5_values <- apply(matrix(1:n,nrow=1),2,B5_value_f)
    
    total_numerator1 <- Reduce("+", lapply(B5_values, `[[`, "numerator1"))
    total_numerator2 <- Reduce("+", lapply(B5_values, `[[`, "numerator2"))
    denominator <-  B8
    B5_trt <- total_numerator1 / denominator
    B5_ctrl <- total_numerator2 / denominator
    return(list(B5_trt=B5_trt,B5_ctrl=B5_ctrl))
  }
  B5 = B5_AOW(covariates, pi_func, cond_prob_trt, cond_prob_ctrl)
  
  ### B6 ### (Dim: 1*1)
  B6_AOW <- function(pi_func, cond_prob_trt, cond_prob_ctrl) {
    B6_value_f = function(j){
      numerator_val1 <- numerator_val2 <- 0
      # Loop through all control indices and calculate contributions to numerator and denominator
      for (i in 1:n) {
        if (i != j) {
          e_Xi <- pi_func[i]
          e_Xj <- pi_func[j]
          h_ij <- e_Xi * (1-e_Xi) * e_Xj * (1-e_Xj)
          mu_ij_trt <- mu_func_trt(i, j, cond_prob_trt, cond_prob_ctrl)
          term1 <- h_ij * mu_ij_trt
          numerator_val1 <- numerator_val1 + term1
          mu_ij_ctrl <- mu_func_ctrl(i, j, cond_prob_trt, cond_prob_ctrl)
          term2 <- h_ij * mu_ij_ctrl
          numerator_val2 <- numerator_val2 + term2
        }
      }
      return(list(numerator1 = numerator_val1, numerator2 = numerator_val2))
    }
    B6_values <- apply(matrix(1:n,nrow=1),2,B6_value_f)
    
    total_numerator1 <- Reduce("+", lapply(B6_values, `[[`, "numerator1"))
    total_numerator2 <- Reduce("+", lapply(B6_values, `[[`, "numerator2"))
    denominator <- B8
    B6_trt <- total_numerator1 / denominator
    B6_ctrl <- total_numerator2 / denominator
    return(list(B6_trt=B6_trt,B6_ctrl=B6_ctrl))
  }
  B6 = B6_AOW(pi_func, cond_prob_trt, cond_prob_ctrl)
  
  ### B7 ### (Dim: 6*1)
  B7_AOW_fun <- function(covariates, pi_func) {
    B7_term1_f = function(j){
      term1_vec <- numeric(p)
      for (i in 1:n) {
        if (i != j) {
          Xi <- as.matrix(covariates[i,])
          term1_vec <- term1_vec + pi_func[i] * (1-pi_func[i]) * ((1-2*pi_func[i]) * Xi)
        }
      }
      term1_vec <- term1_vec * pi_func[j] * (1 - pi_func[j])
      return(list(term1 = term1_vec))
    }
    B7_term1 <- apply(matrix(1:n,nrow=1),2,B7_term1_f)
    
    B7_term2_f = function(j){
      Xj <- as.matrix(covariates[j,])
      term2_vec <- numeric(p)
      for (i in 1:n) {
        if (i != j) {
          term2_vec <- term2_vec + (pi_func[i] * (1-pi_func[i]))
        }
      }
      term2_vec <- term2_vec * pi_func[j] * (1 - pi_func[j]) * (1 - 2*pi_func[j]) * Xj
      return(list(term2 = term2_vec))
    }
    B7_term2 <- apply(matrix(1:n,nrow=1),2,B7_term2_f)
    
    total_term1 <- Reduce("+", lapply(B7_term1, `[[`, "term1"))
    total_term2 <- Reduce("+", lapply(B7_term2, `[[`, "term2"))
    return(total_term1 + total_term2)  
  }
  B7 = B7_AOW_fun(covariates,pi_func)
  
  ### Bn ### (Dim: 6*1)
  Bn_AOW_trt_fun <- function(){
    B1_trt = B1$B1_trt; B2_trt = B2$B2_trt
    B5_trt = B5$B5_trt; B6_trt = B6$B6_trt
    Bn_trt <- B1_trt-B3*(B2_trt/B4)+B5_trt-B7*(B6_trt/B8)
    return(t(Bn_trt))
  }
  
  Bn_AOW_ctrl_fun <- function(){
    B1_ctrl = B1$B1_ctrl; B2_ctrl = B2$B2_ctrl
    B5_ctrl = B5$B5_ctrl; B6_ctrl = B6$B6_ctrl
    Bn_ctrl <- B1_ctrl-B3*(B2_ctrl/B4)+B5_ctrl-B7*(B6_ctrl/B8)
    return(t(Bn_ctrl))
  }
  ### E_beta0beta0 ### (Dim: 6*6)
  Ebeta0 <- function() {
    X <- as.matrix(covariates)
    ps <- pi_func
    Ebeta = crossprod(sqrt(ps*(1-ps)) * X) / n
    return(Ebeta)
    
  }
  
  ### l_beta ### (Dim: 6*1)
  inv_Ebeta0 = solve(Ebeta0())
  l_beta_fun <- function(i) {
    X_i <- as.matrix(covariates[i,])
    Z_i <- treatments[i]
    pi_i <- pi_func[i]
    l_beta <- inv_Ebeta0 %*% X_i * (Z_i - pi_i)
    return(l_beta)
  }
  
  ### E_theta_theta
  thresholds_trt <- dr_trt$zeta
  thresholds_ctrl <- dr_ctrl$zeta
  
  E_thetatheta_fun <- function(outcomes, treatments, covariates_reg, cond_prob_trt, cond_prob_ctrl, thresholds_trt, thresholds_ctrl) {
    F_matrices <- lapply(1:n, function(i) {
      Xi = as.matrix(covariates_reg[i,])
      thresholds <- if(treatments[i] == 1) thresholds_trt else thresholds_ctrl
      cond_prob <- if(treatments[i] == 1) cond_prob_trt else cond_prob_ctrl
      # Calculate delta values based on the outcome for individual i
      delta <- as.numeric(outcomes[i] == 1:L) #which is a 1*L vector
      if(L==3){
        # Calculate the c values for individual i
        c1 <- cond_prob[i, 1]
        c2 <- cond_prob[i, 1] + cond_prob[i,2]
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
      }else{ #in this case L>3
        k = L-1+po
        FM = matrix(0,k,k)
        #the first row
        c1 = cond_prob[i,1]
        f11 = -sum(delta[1:2])*(c1-c1^2)-delta[2]*exp(thresholds[2]-thresholds[1])*((exp(thresholds[2]-thresholds[1])-1)^(-2))
        f12 = delta[2]*exp(thresholds[2]-thresholds[1])*((exp(thresholds[2]-thresholds[1])-1)^(-2))
        #3:(L-1) columns
        f1lm1 = rep(0,L-3)
        f1l = sum(delta[1:2])*(c1 - c1^2)*t(Xi) #1*po
        FM[1,] = c(f11,f12,f1lm1,f1l)
        #2:(L-2) row
        for(l in 2:(L-2)){
          for(j in 1:k){
            if (j < l) FM[l,j] = FM[j,l]
            else if (j == l){
              cl = sum(cond_prob[i,1:l])
              temp1 = sum(delta[l:(l+1)])*(cl-cl^2)
              temp2 = delta[l]*exp(thresholds[l-1]-thresholds[l])*(1-exp(thresholds[l-1]-thresholds[l]))^(-2)
              temp3 = delta[l+1]*exp(thresholds[l+1]-thresholds[l])*((exp(thresholds[l+1]-thresholds[l])-1)^(-2))
              FM[l,j] = -temp1-temp2-temp3
            }
            else if (j==(l+1)) FM[l,j] = temp3
            else if (j<=(L-1)) FM[l,j] = 0
            else FM[l,L:k] = temp1*t(Xi)
          }
        }
        #the L-1 row
        for(j in 1:(L-2)){
          FM[L-1,j] = FM[j,L-1]
        }
        cLm1 = sum(cond_prob[i,1:(L-1)])
        temp1 = sum(delta[(L-1):L])*(cLm1-cLm1^2)
        temp2 = delta[L-1]*exp(thresholds[L-2]-thresholds[L-1])*(1-exp(thresholds[L-2]-thresholds[L-1]))^(-2)
        FM[L-1,L-1] = -temp1-temp2
        FM[L-1,L:k] = temp1*t(Xi)
        #row--L:k, 1:(L-1) columns
        for(l in L:k){
          for(j in 1:(L-1)){
            FM[l,j] = FM[j,l]
          }
        }
        #row--L:k, L:k columns
        factor_val = 0
        for(j in 1:(L-1)){
          cj = sum(cond_prob[i,1:j])
          temp = sum(delta[j:(j+1)])*(cj-cj^2)
          factor_val = factor_val+temp
        }
        FM[L:k,L:k] = -factor_val* (Xi %*% t(Xi))
        F_matrix = FM
      }
      return(F_matrix)
    })
    # Compute the average of the F matrices across all individuals
    E_thetatheta <- Reduce("+", F_matrices) / n
    return(E_thetatheta)
  }
  E_thetatheta <- E_thetatheta_fun(outcomes, treatments, covariates_reg, cond_prob_trt, cond_prob_ctrl, thresholds_trt, thresholds_ctrl)
  
  
  ### mu(Xi, Xj) ### (Dim: 8*1)
  mu_tilde_AOW_fun <- function(i, j, covariates_reg, cond_prob_trt, cond_prob_ctrl) {
    # Extract Xi and Xj as row vectors
    Xi <- as.matrix(covariates_reg[i,])
    Xj <- as.matrix(covariates_reg[j,])
    if(L==3){
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
    }else{
      mu_tilde_vector_l = matrix(0,L-1,1)
      #the first element
      ci1 = cond_prob_trt[i,1]; ci2 = sum(cond_prob_trt[i,1:2])
      cj1 = cond_prob_ctrl[j,1];
      mu_tilde_vector_l[1] = -(ci1-ci1^2)*cj1+(ci2-ci1)*(cj1-cj1^2)
      #the 2 to (L-2) element
      for(l in 2:(L-2)){
        cil = sum(cond_prob_trt[i,1:l]); cilp1 = sum(cond_prob_trt[i,1:(l+1)])
        cjlm1 = sum(cond_prob_ctrl[j,1:(l-1)]); cjl = sum(cond_prob_ctrl[j,1:l])
        mu_tilde_vector_l[l] = (cil-cil^2)*(cjlm1-cjl)+(cilp1-cil)*(cjl-cjl^2)
      }
      #the (L-1)th element
      ciLm1 = sum(cond_prob_trt[i,1:(L-1)]); cjLm2 = sum(cond_prob_ctrl[j,1:(L-2)]) 
      cjLm1 = sum(cond_prob_ctrl[j,1:(L-1)])
      mu_tilde_vector_l[L-1] = (ciLm1-ciLm1^2)*(cjLm2-cjLm1)+(1-ciLm1)*(cjLm1-cjLm1^2)
      #calculate the final term
      factor1 = factor2 = 0
      for(l in 1:(L-2)){
        cil = sum(cond_prob_trt[i,1:l]); cilp1 = sum(cond_prob_trt[i,1:(l+1)])
        cjl = sum(cond_prob_ctrl[j,1:l])
        temp1 = ((cil-cil^2)-(cilp1-cilp1^2))*cjl
        factor1 = factor1 + temp1
        temp2 = (cilp1-cil)*(cjl-cjl^2)
        factor2 = factor2 + temp2
      }
      ciLm1 = sum(cond_prob_trt[i,1:(L-1)]) 
      cjLm1 = sum(cond_prob_ctrl[j,1:(L-1)])
      factor1 = factor1 + (ciLm1-ciLm1^2)*cjLm1
      factor2 = factor2 + (1-ciLm1)*(cjLm1-cjLm1^2)
      mu_tilde_vector_final = factor1*Xi - factor2*Xj
      mu_tilde_vector = rbind(mu_tilde_vector_l,mu_tilde_vector_final)
    }
    return(mu_tilde_vector)
  }
  
  
  ### Cn ### (Dim: 8*1)
  Cn_AOW_fun <- function(treatments, covariates_reg, pi_func, cond_prob_trt, cond_prob_ctrl) {
    # Calculate the first term of Cn
    Cn_first_term_f = function(j){
      numerator_vec <- 0
      # Loop through all control indices and calculate contributions to numerator and denominator
      for (i in 1:n) {
        e_Xi <- pi_func[i]
        e_Xj <- pi_func[j]
        mu_tilde_ij <- mu_tilde_AOW_fun(i, j, covariates_reg, cond_prob_trt, cond_prob_ctrl)
        term <- (treatments[i] * (1 - treatments[j])) * (1-e_Xi) * e_Xj * (-mu_tilde_ij)
        numerator_vec <- numerator_vec + term
      }
      return(list(numerator = numerator_vec))
    }
    Cn_first_term <- apply(matrix(1:n,nrow=1),2,Cn_first_term_f)
    
    term1_numerator <- Reduce("+", lapply(Cn_first_term, `[[`, "numerator"))
    term1_denominator <- adjusted_denom1_term
    first_term <- term1_numerator / term1_denominator
    
    # Second term of Cn
    Cn_second_term_f = function(j){
      numerator_vec <- 0
      # Loop through all control indices and calculate contributions to numerator and denominator
      for (i in 1:n) {
        if (i != j) {
          e_Xi <- pi_func[i]
          e_Xj <- pi_func[j]
          h_ij <- e_Xi * (1-e_Xi) * e_Xj * (1-e_Xj)
          mu_tilde_ij <- mu_tilde_AOW_fun(i, j, covariates_reg, cond_prob_trt, cond_prob_ctrl)
          term <- h_ij * (mu_tilde_ij)
          numerator_vec <- numerator_vec + term
        }
      }
      return(list(numerator = numerator_vec))
    }
    Cn_second_term <- apply(matrix(1:n,nrow=1),2,Cn_second_term_f)
    
    term2_numerator <- Reduce("+", lapply(Cn_second_term, `[[`, "numerator"))
    term2_denominator <- B8
    second_term <- term2_numerator / term2_denominator
    Cn_value <- first_term + second_term
    return(Cn_value)
  }
  # Use the function to calculate Cn
  Cn_value <- Cn_AOW_fun(treatments, covariates_reg, pi_func, cond_prob_trt, cond_prob_ctrl)
  
  E_thetatheta_inv <- ginv(E_thetatheta)
  ### l_tilde_theta ### (Dim: 8*1)
  l_theta_AOW_fun <- function(i,outcomes,treatments,covariates_reg, cond_prob_trt, cond_prob_ctrl, 
                              thresholds_trt, thresholds_ctrl, E_thetatheta_inv) {
    
    Xi <- as.matrix(covariates_reg[i,])
    cond_prob <- if(treatments[i] == 1) cond_prob_trt else cond_prob_ctrl
    delta <- as.numeric(outcomes[i] == 1:L)
    thresholds <- if (treatments[i] == 1) thresholds_trt else thresholds_ctrl
    if(L==3){
      c1 <- cond_prob[i, 1]
      c2 <- c1 + cond_prob[i, 2]
      # Construct the l_tilde vector according to the formula
      l_tilde_vector_1 <- delta[1] - (delta[1] + delta[2]) * c1 - delta[2] * (exp(thresholds[2] - thresholds[1])-1)^(-1)
      l_tilde_vector_2 <- -(delta[2] + delta[3]) * c2 + delta[2] * (1 - exp(thresholds[1] - thresholds[2]))^(-1)
      l_tilde_vector_3 <- (-(delta[1] + delta[2]) + (delta[1] + delta[2]) * c1 + (delta[2] + delta[3]) * c2) * Xi
      l_tilde_vector <- rbind(l_tilde_vector_1, l_tilde_vector_2, l_tilde_vector_3)
      
    }else{
      l_tilde_vector_l = matrix(0,L-1,1)
      c1 = cond_prob[i,1]
      l_tilde_vector_l[1] = delta[1]-sum(delta[1:2])*c1-delta[2]*(exp(thresholds[2]-thresholds[1])-1)^(-1)
      for(l in 2:(L-2)){
        cl = sum(cond_prob[i,1:l])
        l_tilde_vector_l[l] = -sum(delta[l:(l+1)])*cl+delta[l]*1/(1-exp(thresholds[l-1]-thresholds[l]))-delta[l+1]*1/(exp(thresholds[l+1]-thresholds[l])-1)
      }
      cLm1 = sum(cond_prob[i,1:(L-1)])
      l_tilde_vector_l[L-1] = -sum(delta[(L-1):L])*cLm1+delta[L-1]*1/(1-exp(thresholds[L-2]-thresholds[L-1]))
      factor_val = 0
      for(l in 1:(L-1)){
        cl = sum(cond_prob[i,1:l])
        temp = sum(delta[l:(l+1)])*cl
        factor_val = factor_val+temp
      }
      factor_val = factor_val-sum(delta[1:(L-1)])
      l_tilde_vector_finall = factor_val*Xi
      l_tilde_vector = rbind(l_tilde_vector_l,l_tilde_vector_finall)
    }
    l_tilde_value <- E_thetatheta_inv %*% l_tilde_vector
    return(l_tilde_value)
  }
  
  AOW_Variance_trt <- function() {
    var_sum <- 0
    Bn_t <- Bn_AOW_trt_fun()
    Cn_t <- t(Cn_value)
    for (i in 1:n) {
      l_theta <- l_theta_AOW_fun(i,outcomes,treatments,covariates_reg,cond_prob_trt, cond_prob_ctrl, 
                                 thresholds_trt, thresholds_ctrl, E_thetatheta_inv)
      var_sum <- var_sum+(2*(g1_AOW_trt[i]-WP_trt_aow)+Bn_t%*%l_beta_fun(i)+Cn_t%*%l_theta)^2
    }
    AOW_Var_trt <- var_sum/n
    return(AOW_Var_trt)
  }
  
  AOW_Variance_ctrl <- function() {
    var_sum <- 0
    Bn_t <- Bn_AOW_ctrl_fun()
    Cn_t <- t(Cn_value)
    for (i in 1:n) {
      l_theta <- l_theta_AOW_fun(i,outcomes,treatments,covariates_reg,cond_prob_trt, cond_prob_ctrl, 
                                 thresholds_trt, thresholds_ctrl, E_thetatheta_inv)
      var_sum <- var_sum+(2*(g1_AOW_ctrl[i]-WP_ctrl_aow)+Bn_t %*% l_beta_fun(i)+Cn_t %*% l_theta)^2
    }
    AOW_Var_ctrl <- var_sum/n
    return(AOW_Var_ctrl)
  }
  
  #covariance 
  AOW_cov = function(){
    cov_sum <- 0
    Bn_t_trt <- Bn_AOW_trt_fun()
    Bn_t_ctrl <- Bn_AOW_ctrl_fun()
    Cn_t <- t(Cn_value)
    for (i in 1:n) {
      l_theta <- l_theta_AOW_fun(i,outcomes,treatments,covariates_reg,cond_prob_trt, cond_prob_ctrl, 
                                 thresholds_trt, thresholds_ctrl, E_thetatheta_inv)
      l_beta <- l_beta_fun(i)
      temp_trt = (2*(g1_AOW_trt[i]-WP_trt_aow) + Bn_t_trt %*% l_beta + Cn_t %*% l_theta)
      temp_ctrl = (2*(g1_AOW_ctrl[i]-WP_ctrl_aow) + Bn_t_ctrl %*% l_beta + Cn_t %*% l_theta)
      cov_sum <- cov_sum + temp_trt*temp_ctrl
    }
    AOW_cov_results <- cov_sum / n
    return(AOW_cov_results)
  }
  
  Var_trt_aow <- AOW_Variance_trt()/n
  Var_ctrl_aow <- AOW_Variance_ctrl()/n
  cov_aow <- AOW_cov()/n
  ###Theoretical Variance of WR and WD
  Var_WR_aow <- var_ratio(WP_trt_aow,WP_ctrl_aow,Var_trt_aow,Var_ctrl_aow,cov_aow)
  Var_WD_aow <- Var_trt_aow+Var_ctrl_aow-2*cov_aow
  #################AOW Estimator Output#########################
  Estimation = c(WP_trt_aow,WP_ctrl_aow,WR_aow,WD_aow)
  SE = c(sqrt(Var_trt_aow),sqrt(Var_ctrl_aow),sqrt(Var_WR_aow),sqrt(Var_WD_aow))
  conf_low = Estimation-qnorm(0.975)*SE
  conf_high = Estimation+qnorm(0.975)*SE
  res = data.frame(EST=Estimation,SE=SE,conf_low=conf_low,conf_high=conf_high)
  row.names(res)=c("WP","LP","WR","WD")
  return(res)
}

