rm(list = ls())

load("truth_fin_res_estimand_GLMM_nclock_1000_rangeN_1000_1200_M_1000")

fin_res$res_rec

#   alpha        DE SEavg_log_mu1_diff SEavg_log_mu0_diff OEavg_log_mu_diff TEavg_log_mu10_diff
# 1   0.2 0.1740945          1.1583838          1.1014316         1.2265457          0.19175319
# 2   0.3 0.1655351          1.0000000          1.0000000         1.0000000          0.16553511
# 3   0.4 0.1572585          0.8607729          0.9060757         0.8012189          0.14248814
# 4   0.5 0.1494499          0.7394463          0.8190323         0.6279141          0.12240433
# 5   0.6 0.1420051          0.6338031          0.7388231         0.4781967          0.10491667
# 6   0.7 0.1348482          0.5418589          0.6651675         0.3499513          0.08969667
# 7   0.8 0.1278432          0.4620847          0.5983207         0.2412624          0.07649125

load("MC_results_nclock_1000_rangeN_1000_1200_M_100")

MC_results[[1]]$phi_hat_TND$coefs

# true value, realized one
#parameter_vec = c(alpha=NA,beta_0=-0.6,beta_1=0.3,beta_2=0.5,include_RE=T, sigma=0.5)
# beta_0=-0.6,beta_1=0.3,beta_2=0.5, sigma=0.5

mean(sapply(MC_results, function(tbl) tbl$phi_hat_TND$coefs[1]))
# -0.6017117
mean(sapply(MC_results, function(tbl) tbl$phi_hat_TND$coefs[2]))
# 0.2982706
mean(sapply(MC_results, function(tbl) tbl$phi_hat_TND$coefs[3]))
# 0.4997246
mean(sapply(MC_results, function(tbl) tbl$phi_hat_TND$re_var))
# variance 0.2476494

# true values = estimated values

rowMeans(sapply(MC_results, function(tbl) tbl$GM_DE[,2]))
# 0.72728781 0.26323129 0.12422652 0.06890951 0.04772321 0.05527685 0.12270826
# SRS full sample, TND 

#rowMeans(sapply(MC_results, function(tbl) tbl$GM_IEv1[8:14,2]))
#  5.815220e+00 1.000000e+00 1.124599e+00 5.375551e+00 9.680496e+01 4.388561e+03 2.014532e+05

GM_IEv1_res<-rowMeans(sapply(MC_results, function(tbl) tbl$GM_IEv1[,2]))
GM_IEv0_res<-rowMeans(sapply(MC_results, function(tbl) tbl$GM_IEv0[,2]))
# IEv1, IEv0, not good. 
# GM_DE, seems to be okay. 
# check GM_IEv1 and GM_IEv0 

temp_res<-cbind(MC_results[[1]]$GM_IEv1[,c(1,6)],GM_IEv1_res)
cbind(MC_results[[1]]$GM_IEv1[,c(1,6)],GM_IEv1_res)
temp_res[which(temp_res$alpha1==0.3),]

temp_res<-cbind(MC_results[[1]]$GM_IEv0[,c(1,6)],GM_IEv0_res)
cbind(MC_results[[1]]$GM_IEv0[,c(1,6)],GM_IEv0_res)
temp_res[which(temp_res$alpha1==0.3),]


#IEv1 and IEv0 need to check the function.

MC_results[[1]]$GM_IEv0[8:14,]
rowMeans(sapply(MC_results, function(tbl) tbl$GM_IEv0[8:14,2]))
# 2.107131e+00 1.000000e+00 2.378121e+00 2.045495e+01 5.309221e+02 2.072430e+04 4.266074e+05


#   alpha        DE SEavg_log_mu1_diff SEavg_log_mu0_diff OEavg_log_mu_diff TEavg_log_mu10_diff
# 1   0.2 0.1740945          1.1583838          1.1014316         1.2265457          0.19175319
# 2   0.3 0.1655351          1.0000000          1.0000000         1.0000000          0.16553511
# 3   0.4 0.1572585          0.8607729          0.9060757         0.8012189          0.14248814
# 4   0.5 0.1494499          0.7394463          0.8190323         0.6279141          0.12240433
# 5   0.6 0.1420051          0.6338031          0.7388231         0.4781967          0.10491667
# 6   0.7 0.1348482          0.5418589          0.6651675         0.3499513          0.08969667
# 7   0.8 0.1278432          0.4620847          0.5983207         0.2412624          0.07649125


