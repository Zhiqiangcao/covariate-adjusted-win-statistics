library(MASS)
library(nnet)
library(ggplot2)
library(parallel)
library(foreach)
library(doParallel)

######################################################
########## True Win Probability/Win Ratio ############
######################################################
args<-commandArgs(trailingOnly = TRUE)
k <- as.integer(args[1])
if (is.na(k)) k = 1
paste("Scenario:",k)

numCores<-as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK",8))
if (is.na(numCores)) numCores = 8

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
re_count <- 3 #sample size count
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

WP_trt_sim_true <- numeric(sim_num)
WP_ctrl_sim_true <- numeric(sim_num)
WR_sim_true <- numeric(sim_num)

sample_size <- 1000

for (count_temp in 1:sim_num){
  set.seed(count_temp)
  print(count_temp)
  n_count <- sample_size
  # covariates
  x1 <- rnorm(n_count, mean = 1, sd = sds[1])
  x2 <- rnorm(n_count, mean = 0.9, sd = sds[2])
  x3 <- rnorm(n_count, mean = 0.8, sd = sds[3])
  
  x4 <- rbinom(n_count,1,bern_param[1])
  x5 <- rbinom(n_count,1,bern_param[2])
  x6 <- rbinom(n_count,1,bern_param[3])
  
  df_cov <- data.frame(x1, x2, x3, x4, x5, x6)
  treatment_assignment <- rbinom(n_count, 1, 0.5) #0.5 is for balanced design; 0.7 is for unbalanced
  
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

  ############################True Win Ratio##############################
  # ctrl but assigned in trt group
  logodds1_trt_assigned <- b01 + as.matrix(combined_ctrl_quad) %*% c(bi_trt, bi_trt_quad) + trt_eff1
  logodds2_trt_assigned <- b02 + as.matrix(combined_ctrl_quad) %*% c(bi_trt, bi_trt_quad) + trt_eff1
  
  # trt but assigned in ctrl group
  logodds1_ctrl_assigned <- b01 + as.matrix(combined_trt_quad) %*% c(bi_ctrl, bi_ctrl_quad)
  logodds2_ctrl_assigned <- b02 + as.matrix(combined_trt_quad) %*% c(bi_ctrl, bi_ctrl_quad)
  
  ## Probability trt
  prob_2to3_trt_assigned <- inv_logit(logodds1_trt_assigned)
  prob_3_trt_assigned <- inv_logit(logodds2_trt_assigned)
  prob_1_trt_assigned <- 1 - prob_2to3_trt_assigned
  prob_2_trt_assigned <- prob_2to3_trt_assigned - prob_3_trt_assigned
  
  ## Probability ctrl
  prob_2to3_ctrl_assigned <- inv_logit(logodds1_ctrl_assigned)
  prob_3_ctrl_assigned <- inv_logit(logodds2_ctrl_assigned)
  prob_1_ctrl_assigned <- 1 - prob_2to3_ctrl_assigned
  prob_2_ctrl_assigned <- prob_2to3_ctrl_assigned - prob_3_ctrl_assigned
  
  ## outcomes
  outcomes_trt_assigned <- c()
  for (i in 1:nrow(ctrl_cov)) {
    outcomes_trt_assigned[i] <- sample(
      outcomes_3lvl, 
      size = 1,
      prob = c(prob_1_trt_assigned[i], prob_2_trt_assigned[i], prob_3_trt_assigned[i])
    )
  }
  
  outcomes_ctrl_assigned <- c()
  for (i in 1:nrow(trt_cov)) {
    outcomes_ctrl_assigned[i] <- sample(
      outcomes_3lvl, 
      size = 1,
      prob = c(prob_1_ctrl_assigned[i], prob_2_ctrl_assigned[i], prob_3_ctrl_assigned[i])
    )
  }
  
  df_trt_true <- data.frame(c(outcomes_trt, outcomes_trt_assigned))
  df_ctrl_true <- data.frame(c(outcomes_ctrl_assigned, outcomes_ctrl))
  
  # Initialize matrices
  trtwin_true <- matrix(NA, nrow = nrow(df_trt_true), ncol = nrow(df_ctrl_true))
  ctrlwin_true <- matrix(NA, nrow = nrow(df_trt_true), ncol = nrow(df_ctrl_true))
  
  # Define the comparison function for true values
  compare_rows_true <- function(i) {
    trt_row_true <- df_trt_true[i, 1]
    trt_vec_true <- numeric(nrow(df_ctrl_true))
    ctrl_vec_true <- numeric(nrow(df_ctrl_true))
    
    for (j in 1:nrow(df_ctrl_true)) {
      if (trt_row_true > df_ctrl_true[j, 1]) {
        trt_vec_true[j] = 1
        ctrl_vec_true[j] = 0
      } else if (trt_row_true < df_ctrl_true[j, 1]) {
        trt_vec_true[j] = 0
        ctrl_vec_true[j] = 1
      }
    }
    return(list(trt = trt_vec_true, ctrl = ctrl_vec_true))
  }
  # Parallelize the outer loop for true values
  #results_true <- mclapply(1:nrow(df_trt_true), compare_rows_true, mc.cores = numCores)
  results_true <- apply(matrix(1:nrow(df_trt_true),nrow=1),2,compare_rows_true)
  
  # Aggregate results for true values
  for (i in 1:length(results_true)) {
    trtwin_true[i, ] <- results_true[[i]]$trt
    ctrlwin_true[i, ] <- results_true[[i]]$ctrl
  }
  # Calculate win proportions for true values
  trt_winpr_true <- sum(as.vector(trtwin_true), na.rm = TRUE) / (nrow(df_trt_true) * nrow(df_ctrl_true))
  ctrl_winpr_true <- sum(as.vector(ctrlwin_true), na.rm = TRUE) / (nrow(df_trt_true) * nrow(df_ctrl_true))

  WP_trt_sim_true[count_temp] <- trt_winpr_true
  WP_ctrl_sim_true[count_temp] <- ctrl_winpr_true
  WR_sim_true[count_temp] <- trt_winpr_true / ctrl_winpr_true
}


end_time <- Sys.time()
run_time <- end_time - start_time
run_time
#16.78479 mins per 1000 sim_num

######Coverage Output######
WP_trt_trueVal <- mean(WP_trt_sim_true)
WP_ctrl_trueVal <- mean(WP_ctrl_sim_true)
WR_trueVal <- mean(WR_sim_true)
WPWR_TrueValue <- as.data.frame(cbind(WP_trt_trueVal,
                                      WP_ctrl_trueVal,
                                      WR_trueVal))

#balanced design (pi=0.5) with sample size=500
# WP_trt_trueVal WP_ctrl_trueVal WR_trueVal
#      0.2953701       0.2050298   1.447066
#balanced design (pi=0.5) with sample size=1000
#WP_trt_trueVal WP_ctrl_trueVal WR_trueVal
#      0.2957303       0.2046075   1.448601

#unbalanced design (pi=0.7) with sample size=500
#WP_trt_trueVal WP_ctrl_trueVal WR_trueVal
#1      0.2959808       0.2048627    1.45124



















