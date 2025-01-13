# covariate-adjusted-win-statistics

The folders and the R code help to reproducing Tables 1-2 and Figures 1-4 in the article "Covariate-adjusted win statistics in randomized clinical trials with ordinal outcomes" by Zhiqiang Cao, Scott Zuo, Kendra Plourde, Mary Ryan, Patrick Heagerty, Tong, Guangyu and Fan Li (under review)

For questions or comments about the code, please contact Zhiqiang Cao zcaoae@connect.ust.hk. You will need to change the directory to use the example code in script R code. This folder includes the following functions:

1.WR_TureWP.R: generate true win statistics in Table 1;

2.WR_Var_Coverage_UNADJ.R: simulation in Table 1-2 under unadjusted estimator (variance based on multi-sample U-statistics theory derived in Bebu and M.Lachin (2016));

3.WR_Var_Coverage_UNADJ_new.R: simulation in Table 1-2 under unadjusted estimator (variance based on influence function);

4.WR_Var_Coverage_IPW.R: simulation in Table 1-2 under IPW method;

5.WR_Var_Coverage_OW.R: simulation in Table 1-2 under OW method;

6.WR_Var_Coverage_AIPW.R: simulation in Table 1-2 under AIPW method;

7.WR_Var_Coverage_AIPW_mis.R: simulation in Table 1-2 under AIPW method with misspecified outcome regression model;

8.WR_Var_Coverage_AOW.R: simulation in Table 1-2 under AOW method;

9.WR_Var_Coverage_AOW_mis.R: simulation in Table 1-2 under AOW method with misspecified outcome regression model;

10.function_setup_new.R: estimate WP(win probability), LP(loss probability), WR(win ratio), WD(win difference) as well as their SEs and 95% confidence interval
   through unadjusted estimator (two versions of variance estimation), IPW, OW, AIPW and AOW.
   
