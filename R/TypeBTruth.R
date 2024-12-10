install.packages("dplyr") 
library(dplyr)
library(truncnorm)
############################################
############################################
# Causal estimand of Type-B
# num_sub: Number of subsampling approximation

estimand_B<- function(N, alphas, num_sub,
                      OR_1=2, OR_2 = 3, OR_C=1.6,OR_WI=1,OR_WC=5,OR_H=1,em=0){
  estimands <- data.frame(alpha = alphas,
                          mu = 0, mu_1 = 0, mu_0 = 0)
  for(i in seq_len(length(alphas))){
    C<-runif(n=N, 1, 2);   vaxpos <- rtruncnorm(1, -1, 1, mean = 0, sd = 1)
    U1<-rbinom(n=N,size=1,prob=0.5) #affects both
    U2<-rbinom(n=N,size=1,prob=plogis(0.4+0.2*vaxpos)) #affects covid #shifted upwards for vaxpos blocks #it's like reckless behaviour specific to COVID
    # Step2. Treatment model
    ### Subsampling approximation ###
    V <- rbinom(n = num_sub*N, size = 1, prob = rep(alphas[i], num_sub))
    
    g.V = (rep(aggregate(V, list(rep(1:num_sub, each = N)), sum)[,-1], each = N) - V) / (N-1)
    #infection
    #Infection (with something) has some common risk factors U1 and C
    Infp <- plogis(0.25*C-3+1.5*U1)
    Infec<-rbinom(prob=rep(Infp, num_sub),size=1,n=num_sub*N)
    #Infected with COVID #vaccine more effective with high vaccination rates
    C <- rep(C, num_sub); vaxpos <- rep(vaxpos, num_sub); U1 <- rep(U1, num_sub); U2 <- rep(U2, num_sub)
    Infprob<-plogis(-3.5 + C - log(OR_1)*V - log(OR_2)* g.V - log(OR_C)*V*(g.V) +em*V*C + log(2)*U2-2*U1)
    #range(Infprob)
    # which(is.na(Infprob))
    Infec_COVID<-rbinom(prob = rep(Infprob, num_sub), size=1,n=num_sub*N) #0.018
    #symptoms based on infection
    W=W1=W2=rep(0,num_sub*N)
    W1[Infec==1]<-rbinom(prob=plogis(-2+0.5*C[Infec==1]-log(OR_WI)*V[Infec==1]-0.5*U1[Infec==1]),size=1, n=sum(Infec==1))
    W2[Infec_COVID==1]<-rbinom(prob=plogis(-1.5+1*C[Infec_COVID==1]-log(OR_WC)*V[Infec_COVID==1]-1*U1[Infec_COVID==1]+0.5*U2[Infec_COVID==1]*(1-V[Infec_COVID==1])),size=1, n=sum(Infec_COVID))
    W=(W1|W2)
    
    #hospitalization
    H=rep(0,num_sub*N)
    Hprob<-plogis(0.5*C[W==1]-log(OR_H)*V[W==1]-0.5*U1[W==1]+g.V[W==1])
    H[W==1]<-rbinom(prob=Hprob,size=1,n=sum(W==1))
    
    
    mu   = mean(H*Infec_COVID, na.rm = T)
    mu_1 = mean(H*Infec_COVID * V / rep(alphas[i], num_sub), na.rm = T)
    mu_0 = mean(H*Infec_COVID * (1-V) / rep(1-alphas[i], num_sub), na.rm = T)
    
    estimands[i, ] = c(alphas[i], mu, mu_1, mu_0)
    
  }
  return(estimands)
}


estimand_B(N = 10, alphas = c(0.3,0.5,0.7), num_sub = 10000)

############################################
############################################

### Libraries ###
library(dplyr)

### Load estimand computation function ###
source("/Users/congjiang/Documents/simulaINTtnd/Helpfunc_Cor.R")

time = proc.time()

D = 20                                                                         ## Number of simulations
nblocks = 1000                                                                     ## Number of clusters in one simulation
r = 100                                                                     ## Number of subsampling approximation
alphas = seq(0.2, 0.7, 0.05)                                                        ## alpha values of the Type-B policies

print("[Simulation setting]")
print(paste0("D: ", D))
print(paste0("nblocks: ", nblocks))
print(paste0("alphas: ", paste(signif(alphas, 4), collapse = ", ")))
print(paste0("num_sub: ", r))

###----------- Estimand computation  ---------------###

for(simul.id in 11:D){
  
  ### Cluster sizes ###
  N = sample(x = 400:500, size = nblocks, replace = T)
  
  ### Estimands computation ###
  estimands.list <- lapply(N, estimand_B, alphas = alphas, num_sub = r)
  
  #DE: Compute log(mu_1) - log(mu_0) and combine the results in one step
  combined_df <- do.call(rbind, lapply(estimands.list, function(df) {
    df$log_diff <- log(df$mu_1) - log(df$mu_0)
    return(df)
  }))
  
  # Calculate the average log difference for each alphas value
  DEavg_log_diff <- combined_df %>%
    group_by(alpha) %>%
    summarise(DE = exp(mean(log_diff, na.rm = TRUE)))
  
  #SE: Compute the log differences and combine the results
  SE1_log_diff_df <- do.call(rbind, lapply(estimands.list, function(df) {
    df$log_mu1 <- log(df$mu_1)
    standard_log_mu1 <- df$log_mu1[abs(df$alpha - 0.3) < 1e-6] # df$alpha == 0.30
    df$log_mu1_diff <- df$log_mu1 - standard_log_mu1
    return(df)
  }))
  
  # Calculate the average log difference for each alpha value
  SEavg_log_mu1_diff <- SE1_log_diff_df %>%
    group_by(alpha) %>%
    summarise(SE1 = exp(mean(log_mu1_diff, na.rm = TRUE)))
  
  SE0_log_diff_df <- do.call(rbind, lapply(estimands.list, function(df) {
    df$log_mu0 <- log(df$mu_0)
    standard_log_mu0 <- df$log_mu0[abs(df$alpha - 0.3) < 1e-6]
    df$log_mu0_diff <- df$log_mu0 - standard_log_mu0
    return(df)
  }))
  
  # Calculate the average log difference for each alpha value
  SEavg_log_mu0_diff <- SE0_log_diff_df %>%
    group_by(alpha) %>%
    summarise(SE0 = exp(mean(log_mu0_diff, na.rm = TRUE)))
  
  #OE: Compute the log differences and combine the results
  OE_log_diff_df <- do.call(rbind, lapply(estimands.list, function(df) {
    df$log_mu <- log(df$mu)
    standard_log_mu <- df$log_mu[abs(df$alpha - 0.3) < 1e-6]
    df$log_mu_diff <- df$log_mu - standard_log_mu
    return(df)
  }))
  
  # Calculate the average log difference for each alpha value
  OEavg_log_mu_diff <- OE_log_diff_df %>%
    group_by(alpha) %>%
    summarise(OE = exp(mean(log_mu_diff, na.rm = TRUE)))
  
  #TE: Compute the log differences and combine the results
  TE_log_diff_df <- do.call(rbind, lapply(estimands.list, function(df) {
    df$log_mu1 <- log(df$mu_1)
    df$log_mu0 <- log(df$mu_0)
    df$log_mu10_diff <- df$log_mu1 - df$log_mu0[abs(df$alpha - 0.3) < 1e-6]
    return(df)
  }))
  
  # Calculate the average log difference for each alpha value
  TEavg_log_mu10_diff <- TE_log_diff_df %>%
    group_by(alpha) %>%
    summarise(TE = exp(mean(log_mu10_diff, na.rm = TRUE)))
  
  # Combine results into one data frame
  estimands <- DEavg_log_diff %>%
    left_join(SEavg_log_mu1_diff, by = "alpha") %>%
    left_join(SEavg_log_mu0_diff, by = "alpha") %>%
    left_join(OEavg_log_mu_diff, by = "alpha")  %>%
    left_join(TEavg_log_mu10_diff, by = "alpha")
  
  ### Save output ###
  saveRDS(estimands, file = paste0("estimandB_id", simul.id,".rds"))
}

### Stop timer and report total run time ###
script.time = proc.time() - time
print(paste0("Total run time was ", script.time[3], " seconds"))


(readRDS("/n/home09/c55jiang/TNDintRes/estimandB_id10.rds"))

###------ Read estimand ------###
file.list <- list.files("/n/home09/c55jiang/TNDintRes", pattern = "estimandB.*rds")
D <- length(file.list)

# Define the initial data frame
estimandB <- data.frame(alpha = rep(0, length(alphas)),
                        DE = 0, SE_1 = 0, SE_0 = 0,
                        OE = 0, TE = 0)

# Specify the directory containing the .RDS files
directory <- "/n/home09/c55jiang/TNDintRes/"

# Loop through each file in the file list
for (file in file.list) {
  # Read and add the content of each RDS file to estimandB
  estimandB <- estimandB + readRDS(paste0(directory, file))
}

# Print the number of files loaded
print(paste0(length(file.list), " estimand Rdata files were loaded"))

# Calculate the average by dividing by the number of files
estimandB <- estimandB / length(file.list)

# Output the result
print(estimandB)



############################################
############################################
# Causal estimand of GLMM
# num_sub: Number of subsampling approximation


#' Function that calculates the random effect that corresponds to a specific
#' average probability of treatment.
#'
#' @param alpha Target average probability of treatment. Numberic in 0 - 1.
#' @param lin_pred Linear predictor for the probability of treatment among the
#' remaining units. Includes intercept and fixed effects.
#' @param alpha_re_bound The lower and upper end of the values for bi we will
#' look at. Defaults to 10, meaning we will look between - 10 and 10.
#'
FromAlphaToRE <- function(alpha, lin_pred, alpha_re_bound = 10) {
  
  alpha_re_bound <- abs(alpha_re_bound)
  
  r <- optimise(f = AlphaToBi, lower = - 1 * alpha_re_bound,
                upper = alpha_re_bound, alpha = alpha, lin_pred = lin_pred)
  r <- r$minimum
  if (alpha_re_bound - abs(r) < 0.1) {
    warning(paste0('bi = ', r, ', alpha_re_bound = ', alpha_re_bound))
  }
  return(r)
}

AlphaToBi <- function(b, alpha, lin_pred) {
  exp_lin_pred <- exp(lin_pred)
  r <- abs(mean(exp_lin_pred / (exp_lin_pred + exp(- b))) - alpha)
  return(r)
}

estimand_glmm<- function(N, betas, num_sub){
  estimands <- data.frame(beta = betas,
                          mu = 0, mu_1 = 0, mu_0 = 0)
  
  for(i in seq_len(length(betas))){
    C<-runif(n=N, 1, 2);   vaxpos <- rtruncnorm(1, -1, 1, mean = 0, sd = 1)
    U1<-rbinom(n=N,size=1,prob=0.5) #affects both
    U2<-rbinom(n=N,size=1,prob=plogis(0.4+0.2*vaxpos)) #affects covid #shifted upwards for vaxpos blocks #it's like reckless behaviour specific to COVID
    # Step2. Treatment model
    p_trt <- plogis(-0.25 + 0.3*C+0.5*vaxpos)
    
    beta <- betas[i]
    # calculates the random effect that corresponds to a specific
    #' average probability of treatment.
    re_beta <- FromAlphaToRE(alpha = beta, lin_pred = 0.3*C)
    
    pi.beta <- plogis(re_beta + 0.3*C)
    
    ### Subsampling approximation ###
    V <- rbinom(n = num_sub*N, size = 1, prob = rep(pi.beta, num_sub))
    
    g.V = (rep(aggregate(V, list(rep(1:num_sub, each = N)), sum)[,-1], each = N) - V) / (N-1)
    #infection
    #Infection (with something) has some common risk factors U1 and C
    Infp <- plogis(0.2+0.5*C-5+0.5*U1)
    Infec<-rbinom(prob=rep(Infp, num_sub),size=1,n=num_sub*N)
    #Infected with COVID #vaccine more effective with high vaccination rates
    C <- rep(C, num_sub); vaxpos <- rep(vaxpos, num_sub); U1 <- rep(U1, num_sub); U2 <- rep(U2, num_sub)
    Infprob=plogis(-3.5 + C -0.5*V -0.2* g.V - 1.2*log(OR_C)*V*(g.V) +em*V*C +log(2)*U2-2*U1)
    #range(Infprob)
    # which(is.na(Infprob))
    Infec_COVID<-rbinom(prob = rep(Infprob, num_sub), size=1,n=num_sub*N) #0.018
    #symptoms based on infection
    W=W1=W2=rep(0,num_sub*N)
    W1[Infec==1]<-rbinom(prob=plogis(-2+0.5*C[Infec==1]-log(OR_WI)*V[Infec==1]-0.5*U1[Infec==1]),size=1, n=sum(Infec==1))
    W2[Infec_COVID==1]<-rbinom(prob=plogis(-1.5+1*C[Infec_COVID==1]-log(OR_WC)*V[Infec_COVID==1]-1*U1[Infec_COVID==1]+0.5*U2[Infec_COVID==1]*(1-V[Infec_COVID==1])),size=1, n=sum(Infec_COVID))
    W=(W1|W2)+0
    
    #hospitalization
    H=rep(0,num_sub*N)
    Hprob<-plogis(0.5*C[W==1]-log(OR_H)*V[W==1]-0.5*U1[W==1]+g.V[W==1])
    H[W==1]<-rbinom(prob=Hprob,size=1,n=sum(W==1))
    
    
    mu   = mean(H*Infec_COVID, na.rm = T)
    mu_1 = mean(H*Infec_COVID * V / rep(pi.beta, num_sub), na.rm = T)
    mu_0 = mean(H*Infec_COVID * (1-V) / rep(1-pi.beta, num_sub), na.rm = T)
    
    estimands[i, ] = c(beta, mu, mu_1, mu_0)
    
  }
  return(estimands)
}


(estglm <- estimand_glmm(N = 10, betas = c(0.3,0.5,0.7), num_sub = 10000))

############################################
############################################

### Libraries ###
library(dplyr)

### Load estimand computation function ###
source("/Users/congjiang/Documents/simulaINTtnd/Helpfunc_Cor.R")

time = proc.time()

D = 10                                                                         ## Number of simulations
nblocks = 100                                                                      ## Number of clusters in one simulation
r = 100                                                                      ## Number of subsampling approximation
betas = c(0.3,0.5,0.7)                                                        ## beta values of the GLMM policies

print("[Simulation setting]")
print(paste0("D: ", D))
print(paste0("nblocks: ", nblocks))
print(paste0("betas: ", paste(signif(betas, 4), collapse = ", ")))
print(paste0("num_sub: ", r))

###----------- Estimand computation  ---------------###

for(simul.id in 1:D){
  
  ### Cluster sizes ###
  N = sample(x = 50:200, size = nblocks, replace = T)
  
  ### Estimands computation ###
  estimands.list <- lapply(N, estimand_glmm, betas = betas, num_sub = r)
  
  estimands <- data.frame(beta = rep(0, length(betas)),
                          mu = 0, mu_1 = 0, mu_0 = 0)
  
  for(i in 1:nblocks){
    estimands = estimands + estimands.list[[i]]
  }
  
  estimands = estimands / nblocks
  
  ### Causal effects computation using beta==0.5 as the standard ###
  standard = estimands %>% filter(beta == 0.5)
  
  estimands = estimands %>% mutate(de  =  mu_1/mu_0,
                                   se_1 =  mu_1/standard$mu_1,
                                   se_0 =  mu_0/standard$mu_0,
                                   oe  =  mu/standard$mu,
                                   te  =  mu_1/standard$mu_0)
  
  ### Save output ###
  saveRDS(estimands, file = paste0("estimandGLMM_id", simul.id,".rds"))
}

### Stop timer and report total run time ###
script.time = proc.time() - time
print(paste0("Total run time was ", script.time[3], " seconds"))


(estimandGLMM <- readRDS("estimandGLMM_id9.rds"))

###------ Read estimand ------###
file.list <- list.files("/Users/congjiang/Documents/simulaINTtnd", pattern = "estimandGLMM.*rds")
D <- length(file.list)

betas = c(0.3,0.5,0.7)
estimandGLMM <- data.frame(beta = rep(0, length(betas)),
                           mu = 0, mu_1 = 0, mu_0 = 0,
                           de = 0, se_1 = 0, se_0 = 0,
                           oe = 0, te = 0)

for(file in file.list){
  estimandGLMM <- estimandGLMM + readRDS(paste0("/Users/congjiang/Documents/simulaINTtnd/", file))
}
print(paste0(D, " estimand Rdata files were loaded"))
estimandGLMM = estimandGLMM / D
print(estimandGLMM)


print(estimandB)
