library(truncnorm)
datagen_int<-function(rangeN = 400:500, nblocks=1000,OR_1=3, OR_2 =10, OR_C=3.6,OR_WI=1,OR_WC=5,OR_H=1,em=0, return_full=FALSE){
  data.list <- list()
  N <- sample(x = rangeN, size = nblocks, replace = TRUE)
  vaxpos <- rtruncnorm(nblocks, -1, 1, mean = 0, sd = 1)
  for(i in 1:nblocks){
    # Step1. Data generation
    #generate data (WITH clustering for U2)
    C<-runif(n=N[i], 1, 2)
    U1<-rbinom(n=N[i],size=1,prob=0.5) #affects both
    U2<-rbinom(n=N[i],size=1,prob=plogis(0.4+0.2*vaxpos[i])) #affects covid #shifted upwards for vaxpos blocks #it's like reckless behaviour specific to COVID
    # Step2. Treatment model
    p_trt <- plogis(-0.25 + 0.3*C+0.5*vaxpos[i])
    V<-rbinom(prob= p_trt,size=1,n=N[i])
    g.V = (sum(V)-V)/(N[i]-1)
    # Step3. Outcome model
    #infection
    #Infection (with something) has some common risk factors U1 and C
    Infec<-rbinom(prob=plogis(0.25*C-3+1.5*U1),size=1,n=N[i]) #current prevalence around 0.007
    #Infected with COVID #vaccine more effective with high vaccination rates
    Infprob=plogis(-3.5 + C - log(OR_1)*V - log(OR_2)* g.V - log(OR_C)*V*(g.V) +em*V*C +log(2)*U2-2*U1)
    #range(Infprob)
    # which(is.na(Infprob))
    Infec_COVID<-rbinom(prob = Infprob, size=1,n=N[i])
    #symptoms based on infection
    W=W1=W2=rep(0,N[i])
    W1[Infec==1]<-rbinom(prob=plogis(-2+0.5*C[Infec==1]-log(OR_WI)*V[Infec==1]-0.5*U1[Infec==1]),size=1, n=sum(Infec==1))
    W2[Infec_COVID==1]<-rbinom(prob=plogis(-1.5+1*C[Infec_COVID==1]-log(OR_WC)*V[Infec_COVID==1]-1*U1[Infec_COVID==1]+0.5*U2[Infec_COVID==1]*(1-V[Infec_COVID==1])),size=1, n=sum(Infec_COVID))
    W<-(W1|W2)
    #hospitalization
    H=rep(0,N[i])
    Hprob<-plogis(0.5*C[W==1]-log(OR_H)*V[W==1]-0.5*U1[W==1]+g.V[W==1])
    H[W==1]<-rbinom(prob=Hprob,size=1,n=sum(W==1))
    #selection on outcome for testing (does not condition on infectious status, just being in the hospital)
    R <- as.vector(na.omit(sample((1:N[i])[H == 1], min(rangeN[1], sum(H == 1)))))
    
    if (is.logical(return_full) == FALSE) {
      data.list[[i]] <- as.data.frame(cbind(Y=Infec_COVID,V=V,C=C,block=i,f_m=g.V))[R,]
    } else {
      data.list[[i]] <- as.data.frame(cbind(Infec_COVID = Infec_COVID, Infec = Infec, H=H, W=W, V=V, C=C,block=i,f_m=g.V))
    }
  }
  return(dplyr::bind_rows(data.list))
}

(datTND<-datagen_int(nblocks=1000))

(datTND<-datagen_int(nblocks=1000, return_full=T))


sum(datTND$Y == 0)/nrow(datTND)
head(datTND)









############################################
############################################
# Causal estimand of CIPS: Incremental Propensity Score Interventions (Kennedy 2019) extension to Clustered Interference setting
# num_sub: Number of subsampling approximation
estimand_cips<- function(N, deltas, num_sub){
  ### Add `delta` = 1 if not in `deltas` ###
  if(!(1 %in% deltas)) deltas <- c(deltas, 1)
  estimands <- data.frame(delta = deltas,
                          mu = 0, mu_1 = 0, mu_0 = 0)
  
  for(i in seq_len(length(deltas))){
    C<-runif(n=N, 1, 2);   vaxpos <- rtruncnorm(1, -1, 1, mean = 0, sd = 1)
    U1<-rbinom(n=N,size=1,prob=0.5) #affects both
    U2<-rbinom(n=N,size=1,prob=plogis(0.4+0.2*vaxpos)) #affects covid #shifted upwards for vaxpos blocks #it's like reckless behaviour specific to COVID
    # Step2. Treatment model
    p_trt <- plogis(-0.25 + 0.3*C+0.5*vaxpos)
    
    delta <- deltas[i]
    pi.delta <- delta*p_trt / (delta*p_trt + 1 - p_trt)
    
    ### Subsampling approximation ###
    V <- rbinom(n = num_sub*N, size = 1, prob = rep(pi.delta, num_sub))
    
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
    W=(W1|W2)
    
    #hospitalization
    H=rep(0,num_sub*N)
    Hprob<-plogis(0.5*C[W==1]-log(OR_H)*V[W==1]-0.5*U1[W==1]+g.V[W==1])
    H[W==1]<-rbinom(prob=Hprob,size=1,n=sum(W==1))
    
    
    mu   = mean(H*Infec_COVID, na.rm = T)
    mu_1 = mean(H*Infec_COVID * V / rep(pi.delta, num_sub), na.rm = T)
    mu_0 = mean(H*Infec_COVID * (1-V) / rep(1-pi.delta, num_sub), na.rm = T)
    
    estimands[i, ] = c(delta, mu, mu_1, mu_0)
    
  }
  return(estimands)
}


estimand_int(N = 10, deltas = c(0.5,1,2), num_sub = 10000)

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
deltas <- c(0.5,1,2)                                                            ## delta values of the IPS policies

print("[Simulation setting]")
print(paste0("D: ", D))
print(paste0("nblocks: ", nblocks))
print(paste0("deltas: ", paste(signif(deltas, 4), collapse = ", ")))
print(paste0("num_sub: ", r))

###----------- Estimand computation  ---------------###

for(simul.id in 1:D){
  
  ### Cluster sizes ###
  N = sample(x = 50:200, size = nblocks, replace = T)
  
  ### Estimands computation ###
  estimands.list <- lapply(N, estimand_cips, deltas = deltas, num_sub = r)
  
  #DE: Compute log(mu_1) - log(mu_0) and combine the results in one step
  combined_df <- do.call(rbind, lapply(estimands.list, function(df) {
    df$log_diff <- log(df$mu_1) - log(df$mu_0)
    return(df)
  }))
  
  # Calculate the average log difference for each delta value
  DEavg_log_diff <- combined_df %>%
    group_by(delta) %>%
    summarise(DE = exp(mean(log_diff, na.rm = TRUE)))
  
  #SE: Compute the log differences and combine the results
  SE1_log_diff_df <- do.call(rbind, lapply(estimands.list, function(df) {
    df$log_mu1 <- log(df$mu_1)
    standard_log_mu1 <- df$log_mu1[df$delta == 1]
    df$log_mu1_diff <- df$log_mu1 - standard_log_mu1
    return(df)
  }))
  
  # Calculate the average log difference for each delta value
  SEavg_log_mu1_diff <- SE1_log_diff_df %>%
    group_by(delta) %>%
    summarise(SE1 = exp(mean(log_mu1_diff, na.rm = TRUE)))
  
  SE0_log_diff_df <- do.call(rbind, lapply(estimands.list, function(df) {
    df$log_mu0 <- log(df$mu_0)
    standard_log_mu0 <- df$log_mu0[df$delta == 1]
    df$log_mu0_diff <- df$log_mu0 - standard_log_mu0
    return(df)
  }))
  
  # Calculate the average log difference for each delta value
  SEavg_log_mu0_diff <- SE0_log_diff_df %>%
    group_by(delta) %>%
    summarise(SE0 = exp(mean(log_mu0_diff, na.rm = TRUE)))
  
  #OE: Compute the log differences and combine the results
  OE_log_diff_df <- do.call(rbind, lapply(estimands.list, function(df) {
    df$log_mu <- log(df$mu)
    standard_log_mu <- df$log_mu[df$delta == 1]
    df$log_mu_diff <- df$log_mu - standard_log_mu
    return(df)
  }))
  
  # Calculate the average log difference for each delta value
  OEavg_log_mu_diff <- OE_log_diff_df %>%
    group_by(delta) %>%
    summarise(OE = exp(mean(log_mu_diff, na.rm = TRUE)))
  
  #TE: Compute the log differences and combine the results
  TE_log_diff_df <- do.call(rbind, lapply(estimands.list, function(df) {
    df$log_mu1 <- log(df$mu_1)
    df$log_mu0 <- log(df$mu_0)
    df$log_mu10_diff <- df$log_mu1 - df$log_mu0[df$delta == 1]
    return(df)
  }))
  
  # Calculate the average log difference for each delta value
  TEavg_log_mu10_diff <- TE_log_diff_df %>%
    group_by(delta) %>%
    summarise(TE = exp(mean(log_mu10_diff, na.rm = TRUE)))
  
  # Combine results into one data frame
  estimands <- DEavg_log_diff %>%
    left_join(SEavg_log_mu1_diff, by = "delta") %>%
    left_join(SEavg_log_mu0_diff, by = "delta") %>%
    left_join(OEavg_log_mu_diff, by = "delta")  %>%
    left_join(TEavg_log_mu10_diff, by = "delta")
  
  ### Save output ###
  saveRDS(estimands, file = paste0("estimand_id", simul.id,".rds"))
}

### Stop timer and report total run time ###
script.time = proc.time() - time
print(paste0("Total run time was ", script.time[3], " seconds"))


(estimands <- readRDS("estimand_id10.rds"))


library(dplyr)
library(parallel)

### Parameters ###
time <- proc.time()

D <- 10                                # Number of simulations
nblocks <- 100                         # Number of clusters in one simulation
r <- 100                               # Number of subsampling approximation
deltas <- c(0.5, 1, 2)                 # Delta values of the IPS policies

print("[Simulation Setting]")
print(paste0("D: ", D))
print(paste0("nblocks: ", nblocks))
print(paste0("deltas: ", paste(signif(deltas, 4), collapse = ", ")))
print(paste0("num_sub: ", r))

### Function to compute log differences ###
compute_log_differences <- function(estimands_list, column, reference_delta = 1, group_name) {
  df <- do.call(rbind, lapply(estimands_list, function(df) {
    df[[paste0("log_", column)]] <- log(df[[column]])
    reference_log <- df[[paste0("log_", column)]][df$delta == reference_delta]
    df[[paste0("log_", column, "_diff")]] <- df[[paste0("log_", column)]] - reference_log
    return(df)
  }))
  avg_log_diff <- df %>%
    group_by(delta) %>%
    summarise(!!group_name := exp(mean(!!sym(paste0("log_", column, "_diff")), na.rm = TRUE)))
  return(avg_log_diff)
}

### Simulation Function ###
simulate_estimands <- function(simul_id) {
  # Generate cluster sizes
  N <- sample(x = 50:200, size = nblocks, replace = TRUE)
  
  # Estimands computation
  estimands_list <- lapply(N, estimand_cips, deltas = deltas, num_sub = r)
  
  # DE
  DEavg_log_diff <- compute_log_differences(estimands_list, column = "mu_1", group_name = "DE")
  
  # SE1 and SE0
  SEavg_log_mu1_diff <- compute_log_differences(estimands_list, column = "mu_1", group_name = "SE1")
  SEavg_log_mu0_diff <- compute_log_differences(estimands_list, column = "mu_0", group_name = "SE0")
  
  # OE
  OEavg_log_mu_diff <- compute_log_differences(estimands_list, column = "mu", group_name = "OE")
  
  # TE
  TE_log_diff_df <- do.call(rbind, lapply(estimands_list, function(df) {
    df$log_mu1 <- log(df$mu_1)
    df$log_mu0 <- log(df$mu_0)
    df$log_mu10_diff <- df$log_mu1 - df$log_mu0[df$delta == 1]
    return(df)
  }))
  TEavg_log_mu10_diff <- TE_log_diff_df %>%
    group_by(delta) %>%
    summarise(TE = exp(mean(log_mu10_diff, na.rm = TRUE)))
  
  # Combine results
  estimands <- DEavg_log_diff %>%
    left_join(SEavg_log_mu1_diff, by = "delta") %>%
    left_join(SEavg_log_mu0_diff, by = "delta") %>%
    left_join(OEavg_log_mu_diff, by = "delta") %>%
    left_join(TEavg_log_mu10_diff, by = "delta")
  
  # Save results
  saveRDS(estimands, file = paste0("estimand_id", simul_id, ".rds"))
}

### Execute Simulations in Parallel ###
num_cores <- detectCores() - 1
print(paste0("Using ", num_cores, " cores"))
mclapply(1:D, simulate_estimands, mc.cores = num_cores)

### Stop Timer ###
script_time <- proc.time() - time
print(paste0("Total run time was ", script_time[3], " seconds"))



# Check the current working directory
current_directory <- getwd()
print(current_directory)


###------ Read estimand ------###
file.list <- list.files("/Users/congjiang/Documents/simulaINTtnd", pattern = "estimand.*rds")
D <- length(file.list)

deltas <- c(0.5,1,2)
estimands <- data.frame(delta = rep(0, length(deltas)),
                        mu = 0, mu_1 = 0, mu_0 = 0,
                        de = 0, se_1 = 0, se_0 = 0,
                        oe = 0, te = 0)

for(file in file.list){
  estimands <- estimands + readRDS(paste0("/Users/congjiang/Documents/simulaINTtnd/", file))
}
print(paste0(D, " estimand Rdata files were loaded"))
estimands = estimands / D

print(estimands)








############################################
############################################
# Causal estimand of Type-B
# num_sub: Number of subsampling approximation

estimand_B<- function(N, alphas, num_sub){
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
nblocks = 1000                                                                      ## Number of clusters in one simulation
r = 100                                                                      ## Number of subsampling approximation
alphas = c(0.3,0.5,0.7)                                                        ## alpha values of the Type-B policies

print("[Simulation setting]")
print(paste0("D: ", D))
print(paste0("nblocks: ", nblocks))
print(paste0("alphas: ", paste(signif(alphas, 4), collapse = ", ")))
print(paste0("num_sub: ", r))

###----------- Estimand computation  ---------------###

for(simul.id in 1:D){
  
  ### Cluster sizes ###
  N = sample(x = 400:500, size = nblocks, replace = T)
  
  ### Estimands computation ###
  estimands.list <- lapply(N, estimand_B, alphas = alphas, num_sub = r)
  
  estimands <- data.frame(alpha = rep(0, length(alphas)),
                          mu = 0, mu_1 = 0, mu_0 = 0)
  
  for(i in 1:nblocks){
    estimands = estimands + estimands.list[[i]]
  }
  
  estimands = estimands / nblocks
  
  ### Causal effects computation using alpha==0.5 as the standard ###
  standard = estimands %>% filter(alpha == 0.5)
  
  estimands = estimands %>% mutate(de  = mu_1/mu_0,
                                   se_1 = mu_1/standard$mu_1,
                                   se_0 = mu_0/standard$mu_0,
                                   oe  = mu/standard$mu,
                                   te  = mu_1/standard$mu_0)
  
  ### Save output ###
  saveRDS(estimands, file = paste0("estimandB_id", simul.id,".rds"))
}

### Stop timer and report total run time ###
script.time = proc.time() - time
print(paste0("Total run time was ", script.time[3], " seconds"))


(readRDS("estimandB_id10.rds"))

###------ Read estimand ------###
file.list <- list.files("/Users/congjiang", pattern = "estimandB.*rds")
D <- length(file.list)

alphas = c(0.3,0.5,0.7)
estimandB <- data.frame(alpha = rep(0, length(alphas)),
                        mu = 0, mu_1 = 0, mu_0 = 0,
                        de = 0, se_1 = 0, se_0 = 0,
                        oe = 0, te = 0)

for(file in file.list){
  estimandB <- estimandB + readRDS(paste0("/Users/congjiang/", file))
}
print(paste0(D, " estimand Rdata files were loaded"))
estimandB = estimandB / D
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

