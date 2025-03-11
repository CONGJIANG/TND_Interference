#' TND Partial Interference DATA generation
#' 
#'selection probability
# N = sample(x = 50:60, size = nblocks, replace = T)
# tune parameters
# N = sample(x = 20:30, size = nblocks, replace = T)
# num_sub=1000
# OR_1=1, OR_2 = 2
# OR_1=2, OR_2 = 3
# both work
#' 
#' Realized: mixed random effect model with parameters (so only set the value, don't need to resolve it), not necessary keep
#' the overal level as alpha.
#' hypotehtical: fixed effect model or type B 
#' so estimand, full data, based on fixed effect model or type B 
#' realized: TND sample with random effect model. 
#' seperate to generate. 

# change the corresponding lines in this file
# change the corresponding lines function class 
# diagnolize the problem in estimators
# run the simulations 
# list which sets of simulations we need
# real data examples?
# up to function_class.R right now.

rm(list = ls())

library(truncnorm)
library(Matrix)
library(lme4)
library(dplyr)
#*** key for lme4 works.
#utils::install.packages("lme4", type = "source")

source("TNDintIPWestimator.R")
source("function_class.R")

# step 1: generate estimand based on type B or fixed effect model. (super population perspective)
# use this to generate any datatype of data, can generate full, use for estimand 
# can also generate tnd sample, use for estimation.

#treat_type="GLMM2"
#parameter_vec = c(alpha=NA,beta_0=-0.6,beta_1=0.3,beta_2=0.5,include_RE=T, sigma=0.5)
#treat_type="GLMM1"
#parameter_vec=c(alpha= 0.3,beta_0=NA, beta_1=0.3, beta_2=0.5,include_RE=F, sigma=0.5)
# return_full=F
# rangeN = 400:500
# nblocks=1000
# treat_type<-"GLMM2"
# parameter_vec = c(alpha=NA,beta_0=-0.6,beta_1=0.3,beta_2=0.5,include_RE=T, sigma=0.5)

num_sub<-1000
alpha_vec<-seq(0.2,0.8,by=0.1)
nblocks<-2000
rangeN<-5000:8000
estimands<-matrix(0,nrow=nblocks,ncol=4)
estimands<-data.frame(estimands)
colnames(estimands)<-c("mu","mu_1","mu_0","block")
estimands$block<-1:nblocks
estimands_list<-lapply(1:length(alpha_vec), function(i) estimands )

treat_type<-"GLMM"


for(i in 1:length(alpha_vec)){
  #"GLMM1"
  parameter_vec<-c(alpha=alpha_vec[i],beta_0=NA, beta_1=0.3, beta_2=0.5,include_RE=F, sigma=0.5)
  
  # the size of clusters is random 
  for(m in 1:num_sub){
    # for each repetition, 
    # we have 2000 clusters, each of them randomly have 5000-8000 units 
    
    set.seed(m)
    if (m %% 10 == 0) {
      print(paste("Iteration:", m))
    }
    
    
    data <- datagen_int(rangeN = rangeN, nblocks = nblocks, parameter_vec = parameter_vec)
    full_data <-data$data_full
    
    N<-as.vector(table(full_data$block))
    
    mu<- mu_1<- mu_0<-rep(0, length(N) )
    
    # mu   = mean(H*Infec_COVID, na.rm = T)
    # mu_1 = mean(H*Infec_COVID * V / rep(alphas[i], num_sub), na.rm = T)
    # mu_0 = mean(H*Infec_COVID * (1-V) / rep(1-alphas[i], num_sub), na.rm = T)
    mu <-mu+(full_data %>%group_by(block) %>% summarise(mean_value = mean(H * Y, na.rm = TRUE)))[,2]
    
    mu_1 <-mu_1+(full_data %>%group_by(block) %>% summarise(mean_value = mean( (H * Y * V)/parameter_vec["alpha"] , na.rm = TRUE)))[,2]
    
    mu_0 <-mu_0+(full_data %>%group_by(block) %>% summarise(mean_value = mean( (H * Y * (1-V) )/(1-parameter_vec["alpha"]) , na.rm = TRUE)))[,2]
    
    mu<-mu/num_sub
    mu_1<-mu_1/num_sub
    mu_0<-mu_0/num_sub
    
    estimands_list[[i]] <- estimands_list[[i]][,1:3] + cbind(mu, mu_1, mu_0)
    
  }
  
}

file_path<-paste0("estimand_GLMM_RE_0_alpha_vec_M_1000")

#save(estimands_list, file=file_path)

load("estimand_GLMM_RE_0_alpha_vec_M_1000")
estimands.list<-estimands_list

################################################################## estimand DE

estimands.list <- lapply(estimands.list, function(df) {
  # Identify rows where mu_1 or mu_0 is exactly 0
  # avoid -Inf- (-Inf), which produces NaN
  # here it just said mu_1 and mu_0 =0
  
  # rare_case_rec <- which(abs(df$mu_1) == 0 | abs(df$mu_0) == 0)
  # 
  # # Remove those rows
  # if (length(rare_case_rec) > 0) {
  #   df <- df[-rare_case_rec, ]
  # }
  
  
  rare_case_rec_0 <- which(abs(df$mu_0) == 0)
  
  # Remove those rows, then not have that alpha's value
  # having dataset all rows have 0 exists.
  # if something equal to 0, then change to a pretty small number 
  if(length(rare_case_rec_0) > 0 & length(rare_case_rec_0)<length(df$mu_0) ) {
    df$mu_0[rare_case_rec_0]<-10^(-6)
    #df <- df[-rare_case_rec, ]
    df$log_mu0 <- log(df$mu_0)
  }else if(length(rare_case_rec_0)==length(df$mu_0)){
    
    # if mu_1 all =0 across alphas 
    df$mu_0<-10^(-6)
    df$log_mu0 <- log(df$mu_0)
    
  }else{
    df$log_mu0 <- log(df$mu_0)
  }
  
  rare_case_rec_1 <- which(abs(df$mu_1) == 0)
  
  # Remove those rows, then not have that alpha's value
  # having dataset all rows have 0 exists.
  # if something equal to 0, then change to a pretty small number 
  if(length(rare_case_rec_1) > 0 & length(rare_case_rec_1)<length(df$mu_1) ) {
    df$mu_1[rare_case_rec_1]<-10^(-6)
    #df <- df[-rare_case_rec, ]
    df$log_mu1 <- log(df$mu_1)
  }else if(length(rare_case_rec_1)==length(df$mu_1)){
    
    # if mu_1 all =0 across alphas 
    df$mu_1<-10^(-6)
    df$log_mu1 <- log(df$mu_1)
    
  }else{
    df$log_mu1 <- log(df$mu_1)
  }
  
  return(df)
  
})


######### ********
estimands.list <- lapply(estimands.list, function(df) {
  # Identify rows where mu_1 or mu_0 is exactly 0
  # avoid -Inf- (-Inf), which produces NaN
  # here it just said mu_1 and mu_0 =0
  
  # rare_case_rec <- which(abs(df$mu_1) == 0 | abs(df$mu_0) == 0)
  # 
  # # Remove those rows
  # if (length(rare_case_rec) > 0) {
  #   df <- df[-rare_case_rec, ]
  # }
  
  # Compute log difference
  df$log_diff <- log(df$mu_1) - log(df$mu_0)
  
  return(df)
})
############## *********


# dim(combined_df)= 11000 and 5 
# length(alphas)*length(N) (the number of blocks)
#combined_df <- combined_df[order(combined_df$alpha), ]
#combined_df[432,]
# mu_1=0, mu_0=0.0348

res_rec<-data.frame(alpha=alpha_vec, DE=rep(0,length(alpha_vec)))

res_rec$DE<- unlist(lapply(estimands.list, function(df) {
  # Identify rows where mu_1 or mu_0 is exactly 0
  # avoid -Inf- (-Inf), which produces NaN
  # here it just said mu_1 and mu_0 =0
  
  # rare_case_rec <- which(abs(df$mu_1) == 0 | abs(df$mu_0) == 0)
  # 
  # # Remove those rows
  # if (length(rare_case_rec) > 0) {
  #   df <- df[-rare_case_rec, ]
  # }
  
  # Compute log difference
  res_rec<-exp(mean(df$log_diff, na.rm = TRUE))
  
  return(res_rec)
}))

estimands.list <- lapply(estimands.list, function(df) {
  # if df$mu_1=0, then log(0)=-Inf, change to a fairly small number.
  # then when put back to exp, it will close to 0
  #df$log_mu1[which(df$log_mu1==-Inf)]<--
  # get the one which simulated as alpha super close to 0.3.
  # why set as 0.2 does not work???
  standard_log_mu1 <- estimands.list[[which(abs(alpha_vec-0.3)< 1e-6) ]]$log_mu1 # df$alpha == 0.30
  df$log_mu1_diff <- df$log_mu1 - standard_log_mu1
  return(df)
})

res_rec$SEavg_log_mu1_diff<- unlist(lapply(estimands.list, function(df) {
  # Identify rows where mu_1 or mu_0 is exactly 0
  # avoid -Inf- (-Inf), which produces NaN
  # here it just said mu_1 and mu_0 =0
  
  # rare_case_rec <- which(abs(df$mu_1) == 0 | abs(df$mu_0) == 0)
  # 
  # # Remove those rows
  # if (length(rare_case_rec) > 0) {
  #   df <- df[-rare_case_rec, ]
  # }
  
  # Compute log difference
  res_rec<-exp(mean(df$log_mu1_diff, na.rm = TRUE))
  
  return(res_rec)
}))


estimands.list <- lapply(estimands.list, function(df) {
  standard_log_mu0 <- estimands.list[[which(abs(alpha_vec-0.3)< 1e-6) ]]$log_mu0 # df$alpha == 0.30
  df$log_mu0_diff <- df$log_mu0 - standard_log_mu0
  return(df)
})

res_rec$SEavg_log_mu0_diff<- unlist(lapply(estimands.list, function(df) {
  # Identify rows where mu_1 or mu_0 is exactly 0
  # avoid -Inf- (-Inf), which produces NaN
  # here it just said mu_1 and mu_0 =0
  
  # rare_case_rec <- which(abs(df$mu_1) == 0 | abs(df$mu_0) == 0)
  # 
  # # Remove those rows
  # if (length(rare_case_rec) > 0) {
  #   df <- df[-rare_case_rec, ]
  # }
  
  # Compute log difference
  res_rec<-exp(mean(df$log_mu0_diff, na.rm = TRUE))
  
  return(res_rec)
}))

estimands.list <- lapply(estimands.list, function(df) {
  df$log_mu <- log(df$mu) # df$alpha == 0.30
  return(df)
})


#OE: Compute the log differences and combine the results
estimands.list<- lapply(estimands.list, function(df) {
  standard_log_mu <-  (estimands.list[[which(abs(alpha_vec-0.3)< 1e-6) ]])$log_mu
  df$log_mu_diff <- df$log_mu - standard_log_mu
  return(df)
})


# Calculate the average log difference for each alpha value (?? Difference between OE and TE?)
res_rec$OEavg_log_mu_diff <- unlist(lapply(estimands.list, function(df) {
  # Identify rows where mu_1 or mu_0 is exactly 0
  # avoid -Inf- (-Inf), which produces NaN
  # here it just said mu_1 and mu_0 =0
  
  # rare_case_rec <- which(abs(df$mu_1) == 0 | abs(df$mu_0) == 0)
  # 
  # # Remove those rows
  # if (length(rare_case_rec) > 0) {
  #   df <- df[-rare_case_rec, ]
  # }
  
  # Compute log difference
  res_rec<-exp(mean(df$log_mu_diff, na.rm = TRUE))
  
  return(res_rec)
}))

#OE_log_diff_df %>%
#group_by(alpha) %>%
#summarise(OE = exp(mean(log_mu_diff, na.rm = TRUE)))
#OEavg_log_mu_diff

#TE: Compute the log differences and combine the results (total effect)
estimands.list <- lapply(estimands.list, function(df) {
  df$log_mu10_diff <- df$log_mu1 - (estimands.list[[which(abs(alpha_vec-0.3)< 1e-6) ]])$log_mu0
  return(df)
})

# Calculate the average log difference for each alpha value
res_rec$TEavg_log_mu10_diff <- unlist(lapply(estimands.list, function(df) {
  
  res_rec<-exp(mean(df$log_mu10_diff, na.rm = TRUE))
  
  return(res_rec)
}))


fin_res<-list(estimands.list=estimands.list,res_rec=res_rec )
file_path<-paste0("truth_fin_res_estimand_GLMM_RE_0_alpha_vec_M_1000")
#save(fin_res,file=file_path)

load("truth_fin_res_estimand_GLMM_RE_0_alpha_vec_M_1000")

fin_res<-fin_res$res_rec[,-c(3:5)]


#####################################################################

##### estimation part ###############

# if sigma pretty small, then it will have problem 
# need large sample size.
#Error in if (re_sd > 0) { : missing value where TRUE/FALSE needed
#In addition: There were 50 or more warnings (use warnings() to see the first 50)

set.seed(2025)
treat_type<-"GLMM2"
#parameter_vec=c(alpha=NA,beta_0=-0.25,beta_1=0.3,include_RE=T, sigma=0.3)
parameter_vec = c(alpha=NA,beta_0=-0.6,beta_1=0.3,beta_2=0.5,include_RE=T, sigma=0.5)
# treat_type="Type-B"
#parameter_vec=c(alpha=0.5,beta_0=NA,beta_1=0.3,include_RE=F, sigma=0.3)
# ***the number of cluster larger
# *** need return full, estimation of beta_0, beta_1, sigma close as TND sample return =F
# *** the number of clsters increases, sigma, larger. Estimation sigma...
# nblocks<-1000
nblocks<-2000
rangeN<-5000:8000

#datTND<-datagen_int(OR_1=2, OR_2 =3,nblocks=1000,return_full = "BOTH",treat_type=treat_type,parameter_vec=parameter_vec)

# Load required packages
library(ggplot2)
library(dplyr)

# Number of Monte Carlo replications
#MC_reps <- 2
# Function to generate one MC replication

# record: 
# random effect based on TND sample, correct or not 
# full sample, Barkley, estimate should be correct 
# estimator: hypothetical directly oracle Georgia's paper 
# realized: Barkley's model.
# should get back to hypothetical theoretical value. 

run_MC_replication <- function() {
  
  # estimation part, use Barkeley's model. 
  data <- datagen_int(rangeN = rangeN, nblocks = nblocks, treat_type="GLMM", parameter_vec = parameter_vec)
  
  datTND<-data$data_TND
  datFULL <-data$data_full
  datSRS <-data$data_SRS
  
  # ----------- PS estimates ----------- #
  # what's the difference between phi_hat and gamma_numer? 
  cov_names <- c('C1', 'C2')
  cov_cols <- which(names(datTND) %in% cov_names)
  cov_names <- names(datTND)[cov_cols]
  #
  glm_form <- paste('V ~ (1 | block) +', paste(cov_names, collapse = ' + '))
  # GLMM with random effect, full sample seems okay.
  # TND sample bad, but possibly small clusters
  # large sample seems okay.
  # -0.25, 0.3, sd=0.3
  # ***Tnd, Y=0, TND_IPW
  # what return is the variance 
  phi_hat_TND <- PS_model(datTND, glm_form = glm_form, method = 'TND_IPW')
  # ***full data, ori_IPW
  phi_hat_FULL <- PS_model(datFULL, glm_form = glm_form, method = 'ORG_IPW')
  # ***tnd all data, traditional_IPW
  phi_hat_TRAD <- PS_model(datTND, glm_form = glm_form, method = 'ORG_IPW')
  # ***srs data, srs_IPW
  phi_hat_SRS<- PS_model(datSRS, glm_form = glm_form, method = 'ORG_IPW')
  
  # ---------- Coefficients of counterfactual treatment allocation ----------- #
  #gamma_numer <- Policy_coef(datTND, gamma_form= glm_form, method = 'TND_IPW')
  #hypo_beta_1=0.3
  #parameter_vec<-c(alpha=alpha_vec[i],beta_0=NA,beta_1=0.3,include_RE=F, sigma=0.1)
  gamma_numer<-c(NA,0.3)
  
  # ----------- Calculating the IPW ----------- #
  # when averaging across TND sample, the alpha should be biased 
  # but when checking the trueth, it should keep the proportion as alpha. 
  # Yes. Checked. It is the case.
  #obs_alpha <- aggregate(V ~ block, data = datTND, FUN = function(x) mean(x == 1))
  # alpha_range <- quantile(obs_alpha$V, probs = c(0.2, 0.9))
  # alpha <- seq(alpha_range[1], alpha_range[2], length.out = 5)
  # alpha <- unique(sort(round(c(alpha, 0.34, 0.5), 2)))
  alpha<-alpha_vec
  
  # alpha: the hypothetical policies we are interested in. 
  # output: cluster mu_i(1), mu_i(0) under specific alpha.
  # verify: data generating alpha = estimator alpha, verify if it is the same.
  # dim(datTND)[1]/dim(datFULL)[1] cluster level 
  # ??e.g., numerator georgia's one, only need to select estimand = 'GLMM'?
  # change to type B to estimand = 'GLMM', also fine.
  # realized not need to specify 
  # ??denomenator do not need specify here ?
  
  # dta: TND sample generated by realized treatment rule. Unknown, need to be estimated.
  # 
  resB <- GroupIPW(dta = datTND, cov_cols = cov_cols, phi_hat = phi_hat_TND,
                  alpha = alpha, trt_col = which(names(datTND) == 'V'), out_col = which(names(datTND) == 'Y'),
                  keep_re_alpha = F, estimand = 'GLMM',
                  gamma_numer = gamma_numer)
  # rho_i * hat{mu} = mu

  ygroup = resB$yhat_group
  
  # ----------- Estimates and asymptotic variance of the population average potential----------- #
  Score.est <- CalcScore(dta = datTND, neigh_ind = datTND$block, phi_hat = phi_hat_TND, cov_cols = cov_cols,
                         trt_name = 'V', integral_bound = 10)
  
  # Compute GM_DE results
  GM_DE_result <- GM_DE(ygroup = ygroup, boots = NULL, alpha = alpha, alpha_level = 0.05, scores = Score.est, dta = datTND)
  
  # Compute GM_IE results
  GM_IE_v1 <- GM_IE(ygroupM = ygroup[, 1, ], scores = Score.est)
  GM_IE_v0 <- GM_IE(ygroupM = ygroup[, 2, ], scores = Score.est)

  # Extract estimates and confidence intervals
  GM_DE_summary <- data.frame(alpha = colnames(GM_DE_result),
                              est = GM_DE_result["est", ],
                              var = GM_DE_result["var", ],
                              low_int = GM_DE_result["low_int", ],
                              high_int = GM_DE_result["high_int", ])
  
  GM_IE_sumv1 <- list(); GM_IE_sumv0 <- list()
  for (alpha2 in dimnames(GM_IE_v1)[[3]]) {
    temp <- data.frame(alpha1 = colnames(GM_IE_v1[, , alpha2]),
                       est = GM_IE_v1["est", , alpha2],
                       var = GM_IE_v1["var", , alpha2],
                       LB = GM_IE_v1["LB", , alpha2],
                       UB = GM_IE_v1["UB", , alpha2],
                       alpha2 = alpha2)
    GM_IE_sumv1[[alpha2]] <- temp
  }
  
  for (alpha2 in dimnames(GM_IE_v0)[[3]]) {
    temp <- data.frame(alpha1 = colnames(GM_IE_v0[, , alpha2]),
                       est = GM_IE_v0["est", , alpha2],
                       var = GM_IE_v0["var", , alpha2],
                       LB = GM_IE_v0["LB", , alpha2],
                       UB = GM_IE_v0["UB", , alpha2],
                       alpha2 = alpha2)
    GM_IE_sumv0[[alpha2]] <- temp
  }
  
  GM_IE_sumv1 <- bind_rows(GM_IE_sumv1);  GM_IE_sumv0 <- bind_rows(GM_IE_sumv0)
  
  return(list(GM_DE = GM_DE_summary, GM_IEv1 = GM_IE_sumv1, GM_IEv0 = GM_IE_sumv0,phi_hat_TND=phi_hat_TND,phi_hat_FULL=phi_hat_FULL))
}

set.seed(1)
# Run MC replications and store results
MC_reps<-100
MC_results <- vector("list", MC_reps)

#MC_reps <- 1000  # Example
#results <- vector("list", MC_reps)

for (i in 68:MC_reps) {
  MC_results[[i]] <- run_MC_replication()
  
  # Print every 100 iterations
  if (i %% 10 == 0) {
    print(paste("Iteration:", i))
  }
}

save(MC_results,file="MC_results_M_100")

# using the list of mu, mu_1, mu_0, to take the average on each cluster, to obtain the true value approximation.

############################################################################




















datTND[which(datTND$block == 50),]

# Aggregating the data: Calculate the mean of V being 1 for each 'block'
obs_alpha <- aggregate(V ~ block, data = datTND, FUN = function(x) mean(x == 1))
hist(obs_alpha$V, breaks = 100)
print(paste(sum(obs_alpha$V %in% c(0, 1)), 'all treated/control'))


# Identifying blocks with all 0s or all 1s
blocks_all_treated_control <- obs_alpha$block[obs_alpha$V %in% c(0, 1)]
# Removing rows in datTND where block is all treated or all control
datTND <- datTND[!datTND$block %in% blocks_all_treated_control, ]
obs_alpha_fil <- aggregate(V ~ block, data = datTND, FUN = function(x) mean(x == 1))
hist(obs_alpha_fil$V, breaks = 100)
print(paste(sum(obs_alpha_fil$V %in% c(0, 1)), 'all treated/control'))

# ----------- PS estimates ----------- #
cov_names <- c('C')
cov_cols <- which(names(datTND) %in% cov_names)
cov_names <- names(datTND)[cov_cols]

glm_form <- paste('V ~ (1 | block) +', paste(cov_names, collapse = ' + '))
phi_hat <- PS_model(datTND, glm_form = glm_form, method = 'TND_IPW')
# ---------- Coefficients of counterfactual treatment allocation policy Q----------- #
gamma_numer <- Policy_coef(datTND, gamma_form= glm_form, method = 'TND_IPW')
# ----------- Calculating the IPW ----------- #
# GLMM: estimand = 'GLMM'
obs_alpha <- aggregate(V ~ block, data = datTND, FUN = function(x) mean(x == 1))
alpha_range <- quantile(obs_alpha$V, probs = c(0.2, 0.9))
alpha <- seq(alpha_range[1], alpha_range[2], length.out = 5)
alpha <- unique(sort(round(c(alpha, 0.34, 0.5), 2)))

resB <-GroupIPW(dta = datTND, cov_cols = cov_cols, phi_hat = phi_hat,
                alpha = alpha, trt_col = which(names(datTND) == 'V'), out_col = which(names(datTND) == 'Y'),
                estimand = 'GLMM',
                gamma_numer = gamma_numer)
ygroup = resB$yhat_group
apply(ygroup, c(2, 3), mean, na.rm = TRUE)

ypop <- apply(ygroup, 2, mean)
exp(apply( (log(ygroup[, 2, ]) - log(ygroup[, 1, ])), 2, mean))

# ----------- Estimates and asymptotic variance of the population average potential----------- #
Score.est <- CalcScore(dta = datTND, neigh_ind = datTND$block, phi_hat = phi_hat, cov_cols = cov_cols,
                       trt_name = 'V', integral_bound = 10)

DEvar(ygroup = ygroup, scores = Score.est, dta = datTND)
GM_DE(ygroup = ygroup, boots = NULL, alpha = alpha, alpha_level = 0.05, scores = Score.est, dta = datTND)

(se0 <- GM_IE(ygroupM = ygroup[, 1, ], scores = Score.est))
(se1 <- GM_IE(ygroupM = ygroup[, 2, ], scores = Score.est))

# change: (1) the way to generate treatment, can specify the policy. 
# (2) share the data generating each time, then, select sample to do estimation, whole sample to do 
# the truth. 


# Load required packages
library(ggplot2)
library(dplyr)

# Number of Monte Carlo replications
MC_reps <- 2
# Function to generate one MC replication
run_MC_replication <- function() {
  
  datTND <- datagen_int(nblocks=1000)
  # ----------- PS estimates ----------- #
  cov_names <- c('C')
  cov_cols <- which(names(datTND) %in% cov_names)
  cov_names <- names(datTND)[cov_cols]
  glm_form <- paste('V ~ (1 | block) +', paste(cov_names, collapse = ' + '))
  phi_hat <- PS_model(datTND, glm_form = glm_form, method = 'TND_IPW')
  
  
  # ---------- Coefficients of counterfactual treatment allocation ----------- #
  gamma_numer <- Policy_coef(datTND, gamma_form= glm_form, method = 'TND_IPW')
  
  # ----------- Calculating the IPW ----------- #
  # Type B: estimand = '1'
  obs_alpha <- aggregate(V ~ block, data = datTND, FUN = function(x) mean(x == 1))
  alpha_range <- quantile(obs_alpha$V, probs = c(0.2, 0.9))
  alpha <- seq(alpha_range[1], alpha_range[2], length.out = 5)
  alpha <- unique(sort(round(c(alpha, 0.34, 0.5), 2)))
  
  resB <-GroupIPW(dta = datTND, cov_cols = cov_cols, phi_hat = phi_hat,
                  alpha = alpha, trt_col = which(names(datTND) == 'V'), out_col = which(names(datTND) == 'Y'),
                  estimand = 'GLMM',
                  gamma_numer = gamma_numer)
  ygroup = resB$yhat_group
  
  
  
  # ----------- Estimates and asymptotic variance of the population average potential----------- #
  Score.est <- CalcScore(dta = datTND, neigh_ind = datTND$block, phi_hat = phi_hat, cov_cols = cov_cols,
                         trt_name = 'V', integral_bound = 10)
  
  # Compute GM_DE results
  GM_DE_result <- GM_DE(ygroup = ygroup, boots = NULL, alpha = alpha, alpha_level = 0.05, scores = Score.est, dta = datTND)
  
  # Compute GM_IE results
  GM_IE_v1 <- GM_IE(ygroupM = ygroup[, 1, ], scores = Score.est)
  GM_IE_v0 <- GM_IE(ygroupM = ygroup[, 2, ], scores = Score.est)
  
  # Extract estimates and confidence intervals
  GM_DE_summary <- data.frame(alpha = colnames(GM_DE_result),
                              est = GM_DE_result["est", ],
                              var = GM_DE_result["var", ],
                              low_int = GM_DE_result["low_int", ],
                              high_int = GM_DE_result["high_int", ])
  
  GM_IE_sumv1 <- list(); GM_IE_sumv0 <- list()
  for (alpha2 in dimnames(GM_IE_v1)[[3]]) {
    temp <- data.frame(alpha1 = colnames(GM_IE_v1[, , alpha2]),
                       est = GM_IE_v1["est", , alpha2],
                       var = GM_IE_v1["var", , alpha2],
                       LB = GM_IE_v1["LB", , alpha2],
                       UB = GM_IE_v1["UB", , alpha2],
                       alpha2 = alpha2)
    GM_IE_sumv1[[alpha2]] <- temp
  }
  
  for (alpha2 in dimnames(GM_IE_v0)[[3]]) {
    temp <- data.frame(alpha1 = colnames(GM_IE_v0[, , alpha2]),
                       est = GM_IE_v0["est", , alpha2],
                       var = GM_IE_v0["var", , alpha2],
                       LB = GM_IE_v0["LB", , alpha2],
                       UB = GM_IE_v0["UB", , alpha2],
                       alpha2 = alpha2)
    GM_IE_sumv0[[alpha2]] <- temp
  }
  
  GM_IE_sumv1 <- bind_rows(GM_IE_sumv1);  GM_IE_sumv0 <- bind_rows(GM_IE_sumv0)
  
  return(list(GM_DE = GM_DE_summary, GM_IEv1 = GM_IE_sumv1, GM_IEv0 = GM_IE_sumv0))
}

set.seed(1)
# Run MC replications and store results
MC_reps<-2
MC_results <- replicate(MC_reps, run_MC_replication(), simplify = FALSE)

# Aggregate GM_DE results
GM_DE_all <- bind_rows(lapply(MC_results, `[[`, "GM_DE"))
GM_DE_summary <- GM_DE_all %>%
  group_by(alpha) %>%
  summarise(mean_est = mean(est), sd_est = sd(est), 
            mean_low_int = mean(low_int), mean_high_int = mean(high_int))

# Aggregate GM_IE results
GM_IE_all <- bind_rows(lapply(MC_results, `[[`, "GM_IE"))
GM_IE_sumv1 <- GM_IE_all %>%
  group_by(alpha1, alpha2) %>%
  summarise(mean_est = mean(est), sd_est = sd(est), 
            mean_LB = mean(LB), mean_UB = mean(UB))

# Plot GM_DE estimates with confidence intervals
ggplot(GM_DE_summary, aes(x = as.numeric(alpha), y = mean_est)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean_low_int, ymax = mean_high_int), width = 0.02) +
  labs(title = "GM_DE Estimates with Confidence Intervals",
       x = "Alpha", y = "Estimate") +
  theme_minimal()

# Plot GM_IE estimates with confidence intervals
ggplot(GM_IE_sumv1, aes(x = as.numeric(alpha1), y = mean_est, color = as.factor(alpha2))) +
  geom_point() +
  geom_errorbar(aes(ymin = mean_LB, ymax = mean_UB), width = 0.02) +
  labs(title = "GM_IE Estimates with Confidence Intervals",
       x = "Alpha1", y = "Estimate", color = "Alpha2") +
  theme_minimal()


















r <- 3
# Initialize a list to store results
resDE.org <- vector("list", length = r)
resSE0.org <- vector("list", length = r)
resSE1.org <- vector("list", length = r)

resDE.tnd <- vector("list", length = r)
resSE0.tnd <- vector("list", length = r)
resSE1.tnd <- vector("list", length = r)

# Run the function r times and store the results
for (i in 1:r) {
  tryCatch({
    datTND <- datagen_int(nblocks=1000)
    # ----------- PS estimates ----------- #
    cov_names <- c('C')
    cov_cols <- which(names(datTND) %in% cov_names)
    cov_names <- names(datTND)[cov_cols]
    glm_form <- paste('V ~ (1 | block) +', paste(cov_names, collapse = ' + '))
    phi_hat <- PS_model(datTND, glm_form = glm_form, method = 'TND_IPW')
    
    
    # ---------- Coefficients of counterfactual treatment allocation ----------- #
    gamma_numer <- Policy_coef(datTND, gamma_form= glm_form, method = 'TND_IPW')
    
    # ----------- Calculating the IPW ----------- #
    # Type B: estimand = '1'
    obs_alpha <- aggregate(V ~ block, data = datTND, FUN = function(x) mean(x == 1))
    alpha_range <- quantile(obs_alpha$V, probs = c(0.2, 0.9))
    alpha <- seq(alpha_range[1], alpha_range[2], length.out = 5)
    alpha <- unique(sort(round(c(alpha, 0.34, 0.5), 2)))
    
    resB <-GroupIPW(dta = datTND, cov_cols = cov_cols, phi_hat = phi_hat,
                    alpha = alpha, trt_col = which(names(datTND) == 'V'), out_col = which(names(datTND) == 'Y'),
                    estimand = 'GLMM',
                    gamma_numer = gamma_numer)
    ygroup = resB$yhat_group
    
    # ----------- Estimates and asymptotic variance of the population average potential----------- #
    Score.est <- CalcScore(dta = datTND, neigh_ind = datTND$block, phi_hat = phi_hat, cov_cols = cov_cols,
                           trt_name = 'V', integral_bound = 10)
    resDE.tnd[[i]] <- GM_DE(ygroup = ygroup, boots = NULL, alpha = alpha, alpha_level = 0.05, scores = Score.est, dta = datTND)
    
    #. indices corresponding to treatment = 0,  and [ygroup[, 2, ] is for treatment = 1]
    resSE0.tnd[[i]] <- GM_IE(ygroupM = ygroup[, 1, ], scores = Score.est)
    resSE1.tnd[[i]] <- GM_IE(ygroupM = ygroup[, 2, ], scores = Score.est)
  }, error = function(e) {
    # Print the error message
    cat("Error in iteration", i, ":", conditionMessage(e), "\n")
    # Optionally, you can print more details about the error using conditionCall(e) or traceback()
  })
}






# Initialize lists to store results
resDE.tnd <- vector("list", r)
resSE0.tnd <- vector("list", r)
resSE1.tnd <- vector("list", r)

# Run the function r times and store the results
for (i in 1:r) {
  tryCatch({
    cat("Running iteration", i, "...\n")  # Progress tracking
    
    # Data generation
    datTND <- datagen_int(nblocks = 1000)
    
    # ----------- PS estimates ----------- #
    cov_names <- c('C')
    cov_cols <- which(names(datTND) %in% cov_names)
    cov_names <- names(datTND)[cov_cols]
    glm_form <- paste('V ~ (1 | block) +', paste(cov_names, collapse = ' + '))
    phi_hat <- PS_model(datTND, glm_form = glm_form, method = 'TND_IPW')
    
    # ---------- Coefficients of counterfactual treatment allocation ----------- #
    gamma_numer <- Policy_coef(datTND, gamma_form = glm_form, method = 'TND_IPW')
    
    # ----------- Calculating the IPW ----------- #
    obs_alpha <- aggregate(V ~ block, data = datTND, FUN = function(x) mean(x == 1))
    alpha_range <- quantile(obs_alpha$V, probs = c(0.2, 0.9))
    alpha <- seq(alpha_range[1], alpha_range[2], length.out = 5)
    alpha <- unique(sort(round(c(alpha, 0.34, 0.5), 2)))  # Ensure unique, sorted alphas
    
    # Compute Group IPW
    resB <- GroupIPW(
      dta = datTND, cov_cols = cov_cols, phi_hat = phi_hat,
      alpha = alpha, trt_col = which(names(datTND) == 'V'), 
      out_col = which(names(datTND) == 'Y'),
      estimand = 'GLMM', gamma_numer = gamma_numer
    )
    ygroup <- resB$yhat_group
    
    # ----------- Compute Estimates ----------- #
    Score.est <- CalcScore(
      dta = datTND, neigh_ind = datTND$block, phi_hat = phi_hat, cov_cols = cov_cols,
      trt_name = 'V', integral_bound = 10
    )
    
    # Store results
    resDE.tnd[[i]] <- GM_DE(ygroup = ygroup, boots = NULL, alpha = alpha, alpha_level = 0.05, scores = Score.est, dta = datTND)
    resSE0.tnd[[i]] <- GM_IE(ygroupM = ygroup[, 1, ], scores = Score.est)
    resSE1.tnd[[i]] <- GM_IE(ygroupM = ygroup[, 2, ], scores = Score.est)
    
  }, error = function(e) {
    cat("Error in iteration", i, ":", conditionMessage(e), "\n")
  })
}

cat("Monte Carlo simulations completed!\n")


# Convert lists to data frames
resDE.df <- do.call(rbind, resDE.tnd)
resSE0.df <- do.call(rbind, resSE0.tnd)
resSE1.df <- do.call(rbind, resSE1.tnd)


library(dplyr)

# Function to summarize results
summarize_results <- function(df) {
  df %>%
    group_by(alpha) %>%  # Group by alpha values
    summarise(
      mean_est = mean(estimate, na.rm = TRUE),
      sd_est = sd(estimate, na.rm = TRUE),
      se_est = sd_est / sqrt(n()),  # Standard error
      lower_CI = mean_est - 1.96 * se_est,  # 95% CI lower bound
      upper_CI = mean_est + 1.96 * se_est   # 95% CI upper bound
    )
}

# Summarize for each metric
summary_DE <- summarize_results(resDE.df)
summary_SE0 <- summarize_results(resSE0.df)
summary_SE1 <- summarize_results(resSE1.df)


library(ggplot2)

# Function to create CI plots
plot_results <- function(summary_df, title) {
  ggplot(summary_df, aes(x = alpha, y = mean_est)) +
    geom_point(color = "blue", size = 2) +  # Mean estimates
    geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), width = 0.05, color = "red") +  # CI bars
    labs(title = title, x = "Alpha", y = "Estimate") +
    theme_minimal()
}

# Plot results
plot_DE <- plot_results(summary_DE, "Estimated Direct Effects")
plot_SE0 <- plot_results(summary_SE0, "Estimated Indirect Effects (Treatment = 0)")
plot_SE1 <- plot_results(summary_SE1, "Estimated Indirect Effects (Treatment = 1)")

# Display plots
print(plot_DE)
print(plot_SE0)
print(plot_SE1)


