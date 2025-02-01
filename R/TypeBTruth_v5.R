# need settings on more repetiation on treatment assignments
# second, need to deal with extreme value (change mu_1, mu_0 =0 into very small number.)
# what is the extreme part appear on the estimator? 

# selection probability
# N = sample(x = 50:60, size = nblocks, replace = T)
# tune parameters
# N = sample(x = 20:30, size = nblocks, replace = T)
# num_sub=1000
# OR_1=1, OR_2 = 2
# OR_1=2, OR_2 = 3
# both work

#install.packages("dplyr") 
library(dplyr)
library(truncnorm)
############################################
############################################
# Causal estimand of Type-B
# The function is used for only simulate one cluster with size N
# num_sub: Number of repetitions for treatment realization

# mu(alpha), mu(1,alpha), mu(0,alpha)
estimand_B <- function(N, alphas, num_sub,
                       OR_1=2, OR_2 = 3, OR_C=1.6,OR_WI=1,OR_WC=5,OR_H=1,em=0){
  
  # N: the size of cluster 
  # one block repetition, get E(1/n_j sum_{i=1}^n_j Y_i(alpha))
  # E(1/n_j sum_{i=1}^n_j Y_i(1, alpha)), E(1/n_j sum_{i=1}^n_j Y_i(0, alpha))
  
  # rows for different coverage. 
  # mu_1 and mu_0
  estimands <- data.frame(alpha = alphas,
                          mu = 0, mu_1 = 0, mu_0 = 0)
  # simpler model, keep V, g.V 
  
  for(i in seq_len(length(alphas))){
    
    # process: vaccination -> if infected & having symption & hospitalized, then these people are included in TND sample
    # treated: in TND sample, who are vaccinated, in TND sample, who are not treated. 
    
    C<-runif(n=N, 1, 2);   vaxpos <- rtruncnorm(1, -1, 1, mean = 0, sd = 1)
    U1<-rbinom(n=N,size=1,prob=0.5) #affects both
    U2<-rbinom(n=N,size=1,prob=plogis(0.4+0.2*vaxpos)) #affects covid #shifted upwards for vaxpos blocks #it's like reckless behaviour specific to COVID
    
    # Step2. treatment model
    ### Subsampling approximation: treatment assignment realization ###
    # i.i.d. sample N treatment with probability alphas[i]
    # repeated num_sub times
    V <- rbinom(n = num_sub*N, size = 1, prob = alphas[i])
    # V <- rbinom(n = num_sub*N, size = 1, prob = rep(alphas[i], num_sub))
    
    # element i in g.V represents the proportion of treated neighbors of unit i.
    g.V = (rep(aggregate(V, list(rep(1:num_sub, each = N)), sum)[,-1], each = N) - V) / (N-1)
    
    # Infp: infected based on other disease than Covid-19 liked diseases
    # not influenced by neighbors' treatments.
    Infp <- plogis(0.25*C-3+1.5*U1)
    
    # infection with some diseases.
    Infec<-rbinom(prob=rep(Infp, num_sub),size=1,n=num_sub*N)
    #Infected with COVID #vaccine more effective with high vaccination rates
    # same parameter, across different blocks. 
    C <- rep(C, num_sub); vaxpos <- rep(vaxpos, num_sub); U1 <- rep(U1, num_sub); U2 <- rep(U2, num_sub)
    
    # it counts the effect that myself being vaccinated, my neighbor being vaccinated 
    # interaction effect. 
    
    # appropriate setting such that the probability not very extreme. 
    Infprob<-plogis(-3.5 + C - log(OR_1)*V - log(OR_2)* g.V - log(OR_C)*V*(g.V)+em*V*C + log(2)*U2-2*U1)
    #range(Infprob)
    # which(is.na(Infprob))
    #Infec_COVID<-rbinom(prob = rep(Infprob, num_sub), size=1,n=num_sub*N) #0.018
    ### ******* 
    Infec_COVID<-rbinom(n=num_sub*N, size=1, prob = Infprob)
    ### ******
    
    #symptoms based on infection
    W=W1=W2=rep(0,num_sub*N)
    # infection based on other disease, then less people have sympton
    W1[Infec==1]<-rbinom(prob=plogis(-2+0.5*C[Infec==1]-log(OR_WI)*V[Infec==1]-0.5*U1[Infec==1]),size=1, n=sum(Infec==1))
    # infected by covid, then more people having sympton
    W2[Infec_COVID==1]<-rbinom(prob=plogis(-1.5+1*C[Infec_COVID==1]-log(OR_WC)*V[Infec_COVID==1]-1*U1[Infec_COVID==1]+0.5*U2[Infec_COVID==1]*(1-V[Infec_COVID==1])),size=1, n=sum(Infec_COVID))
    # either infected in by other disease or infected by covid.
    W=(W1|W2)
    
    #hospitalization
    H=rep(0,num_sub*N)
    # (1) having sympton, then possibly hospitalized. 
    # (2) if having sympton but treated, then having smaller prob get hospitalized.
    # neighbor being treated and they themself get sympton, then easier to go hospital.
    Hprob<-plogis(0.5*C[W==1]-log(OR_H)*V[W==1]-0.5*U1[W==1]+g.V[W==1])
    H[W==1]<-rbinom(prob=Hprob,size=1,n=sum(W==1))
    
    # given one block with size N, repeated num_sub times
    # 1/num_sub sum_{m=1}^num_sub 1/n_i [ sum_{j=1}^{n_i} Y_j  ]
    # = 1/(num_sub*n_i) sum_{m=1}^num sum_{j=1}^{n_i} Y_j
    
    mu<- mu_1<- mu_0<-0
    
    for(j in 1:num_sub){
      
      # mu   = mean(H*Infec_COVID, na.rm = T)
      # mu_1 = mean(H*Infec_COVID * V / rep(alphas[i], num_sub), na.rm = T)
      # mu_0 = mean(H*Infec_COVID * (1-V) / rep(1-alphas[i], num_sub), na.rm = T)
      mu <- mu+ mean( H[ ((j-1)*N+1): (j*N) ]*Infec_COVID[((j-1)*N+1): (j*N)] , na.rm = T)
      mu_1<-mu_1+  mean( (H[ ((j-1)*N+1): (j*N) ]*Infec_COVID[((j-1)*N+1): (j*N)]* V[((j-1)*N+1): (j*N)])/alphas[i],na.rm = T )
      mu_0<-mu_0+  mean( (H[ ((j-1)*N+1): (j*N) ]*Infec_COVID[((j-1)*N+1): (j*N)]* (1-V[((j-1)*N+1): (j*N)]) )/(1-alphas[i]),na.rm = T )
      
    }
    
    mu<-mu/num_sub
    mu_1<-mu_1/num_sub
    mu_0<-mu_0/num_sub
    
    estimands[i, ] = c(alphas[i], mu, mu_1, mu_0)
    
  }
  return(estimands)
}

# infection-> sympton -> hospitalized 
# prob(Y_i=1)=P(H=1|W=1,Infec_COVID=1 or infec=1 )* P(W=1|Infec_COVID=1 or infec=1)*(P(Infec_COVID=1|parameter)+P(infec=1|parameter))
# integrate out all the situation, obtain the true values. 

# alpha increase
# OR1: V parameter, OR2: g.V parameter
estimand_B(N = 20, alphas = c(0.3,0.5,0.7), OR_1=2, OR_2 = 3, num_sub = 1000)

estimand_B(N = 20, alphas = c(0.3,0.5,0.7), OR_1=1, OR_2 = 2, num_sub = 1000)


############################################
############################################

### Libraries ###
library(dplyr)

### Load estimand computation function ###
#source("/Users/congjiang/Documents/simulaINTtnd/Helpfunc_Cor.R")

time = proc.time()

D = 10                                                                         ## Number of simulations
nblocks = 100                                                                     ## Number of clusters in one simulation
r = 1000                                                                     ## Number of subsampling approximation
alphas = seq(0.2, 0.7, 0.05)                                                        ## alpha values of the Type-B policies

print("[Simulation setting]")
print(paste0("D: ", D))
print(paste0("nblocks: ", nblocks))
print(paste0("alphas: ", paste(signif(alphas, 4), collapse = ", ")))
print(paste0("num_sub: ", r))

###----------- Estimand computation  ---------------###

for(simul.id in 1:D){
  
  ### Cluster sizes ###
  # each cluster sizes. 
  # total nblocks: having size from 50 to 100 
  N = sample(x = 20:30, size = nblocks, replace = T)
  
  ### Estimands computation ###
  # with each combination of (N,alphas,num_sub), apply the function estimand_B 
  # this is a list, each element in the list is mu(alpha), mu(1,alpha) and mu(0,alpha)
  # for given alpha, each cluster, sample that r times to obtain the mu(alpha), mu(1,alpha) and mu(0,alpha)
  estimands.list <- lapply(N, estimand_B, alphas = alphas, num_sub = 1000,OR_1=2, OR_2 = 3)
  
  #DE: Compute log(mu_1) - log(mu_0) and combine the results in one step
  # problem df$mu_1=0, log(df$mu_1)=-Inf
  # log(df$mu_1) - log(df$mu_0)=-Inf, mean(-Inf), exp(-Inf)=0, i.e., one value with extreme will drive 
  # the whole function to take 0, which you can hard to see the decreasing trend when alpha grows.
  # mu_1=0 means all Y_ij(1,Z_{-ij})=0 for all Z_{-ij}, which is rare, most likely is the approximation error
  
  # combined_df <- do.call(rbind, lapply(estimands.list, function(df) {
  #   # remove those cluster, which should not affect the estimand a lot if there is
  #   # a large number of clusters.
  #   df$log_diff <- log(df$mu_1) - log(df$mu_0)
  #   return(df)
  # }))
  
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
  combined_df <- do.call(rbind, lapply(estimands.list, function(df) {
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
  }))
  ############## *********
  
  
  # dim(combined_df)= 11000 and 5 
  # length(alphas)*length(N) (the number of blocks)
  
  #combined_df <- combined_df[order(combined_df$alpha), ]
  #combined_df[432,]
  # mu_1=0, mu_0=0.0348
  
  # Calculate the average log difference for each alphas value
  # exp{ mean[log(df$mu_1) - log(df$mu_0)] }
  # expectation taken over all clusters
  DEavg_log_diff <- combined_df %>%
    group_by(alpha) %>%
    summarise(DE = exp(mean(log_diff, na.rm = TRUE)))
  DEavg_log_diff
  
  # need to remove the one which has -Inf Inf, why does this happen?
  # exp(mean(combined_df$log_diff[combined_df$alpha==0.2],na.rm=T))
  # exp(mean(combined_df$log_diff[combined_df$alpha==0.2][-432]))
  # combined_df$log_diff[combined_df$alpha==0.2][432]
  
  SE1_log_diff_df <- do.call(rbind, lapply(estimands.list, function(df) {
    
    # if df$mu_1=0, then log(0)=-Inf, change to a fairly small number.
    # then when put back to exp, it will close to 0
    #df$log_mu1[which(df$log_mu1==-Inf)]<--
    # get the one which simulated as alpha super close to 0.3.
    # why set as 0.2 does not work???
    standard_log_mu1 <- df$log_mu1[abs(df$alpha - 0.3) < 1e-6] # df$alpha == 0.30
    df$log_mu1_diff <- df$log_mu1 - standard_log_mu1
    return(df)
  }))
  
  
  # #SE: Compute the log differences and combine the results
  # # comparing with policy 0.3
  # SE1_log_diff_df <- do.call(rbind, lapply(estimands.list, function(df) {
  #   
  #   # Identify rows where mu_1 or mu_0 is exactly 0
  #   # rare_case_rec <- which(abs(df$mu_1) == 0 )
  #   # 
  #   # # Remove those rows
  #   # if (length(rare_case_rec) > 0) {
  #   #   df <- df[-rare_case_rec, ]
  #   # }
  #   
  #   
  #   df$log_mu1 <- log(df$mu_1)
  #   
  #   # get the one which simulated as alpha super close to 0.3.
  #   # why set as 0.2 does not work???
  #   standard_log_mu1 <- df$log_mu1[abs(df$alpha - 0.40) < 1e-6] # df$alpha == 0.30
  #   df$log_mu1_diff <- df$log_mu1 - standard_log_mu1
  #   return(df)
  # }))
  
  # Calculate the average log difference for each alpha value
  SEavg_log_mu1_diff <- SE1_log_diff_df %>%
    group_by(alpha) %>%
    summarise(SE1 = exp(mean(log_mu1_diff, na.rm = TRUE)))
  SEavg_log_mu1_diff
  
  SE0_log_diff_df <- do.call(rbind, lapply(estimands.list, function(df) {
    rare_case_rec <- which(abs(df$mu_0) == 0)
    
    # Remove those rows, then not have that alpha's value
    # having dataset all rows have 0 exists.
    # if something equal to 0, then change to a pretty small number 
    if(length(rare_case_rec) > 0 & length(rare_case_rec)<length(df$mu_0) ) {
      df$mu_0[rare_case_rec]<-10^(-6)
      #df <- df[-rare_case_rec, ]
      df$log_mu0 <- log(df$mu_0)
    }else if(length(rare_case_rec)==length(df$mu_0)){
      
      # if mu_1 all =0 across alphas 
      df$mu_0<-10^(-6)
      df$log_mu0 <- log(df$mu_0)
      
    }else{
      df$log_mu0 <- log(df$mu_0)
    }
    
    #df$log_mu0 <- log(df$mu_0)
    standard_log_mu0 <- df$log_mu0[abs(df$alpha - 0.3) < 1e-6]
    df$log_mu0_diff <- df$log_mu0 - standard_log_mu0
    return(df)
  }))
  
  # Calculate the average log difference for each alpha value
  SEavg_log_mu0_diff <- SE0_log_diff_df %>%
    group_by(alpha) %>%
    summarise(SE0 = exp(mean(log_mu0_diff, na.rm = TRUE)))
  SEavg_log_mu0_diff
  
  #OE: Compute the log differences and combine the results
  OE_log_diff_df <- do.call(rbind, lapply(estimands.list, function(df) {
    df$log_mu <- log(df$mu)
    standard_log_mu <- df$log_mu[abs(df$alpha - 0.3) < 1e-6]
    df$log_mu_diff <- df$log_mu - standard_log_mu
    return(df)
  }))
  
  
  # Calculate the average log difference for each alpha value (?? Difference between OE and TE?)
  OEavg_log_mu_diff <- OE_log_diff_df %>%
    group_by(alpha) %>%
    summarise(OE = exp(mean(log_mu_diff, na.rm = TRUE)))
  OEavg_log_mu_diff
  
  #TE: Compute the log differences and combine the results (total effect)
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


#(readRDS("/n/home09/c55jiang/TNDintRes/estimandB_id10.rds"))

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
#' optimize 1/n sum p_i(gamma) - beta =0, optimize gamma 
#' @param alpha Target average probability of treatment. Numberic in 0 - 1.
#' @param lin_pred Linear predictor for the probability of treatment among the
#' remaining units. Includes intercept and fixed effects.
#' @param alpha_re_bound The lower and upper end of the values for bi we will
#' look at. Defaults to 10, meaning we will look between - 10 and 10.
#'
FromAlphaToRE <- function(alpha, lin_pred, alpha_re_bound = 10) {
  
  alpha_re_bound <- abs(alpha_re_bound)
  
  # the optimal b such that AlphaToBi reach the smallest?
  r <- optimise(f = AlphaToBi, lower = - 1 * alpha_re_bound,
                upper = alpha_re_bound, alpha = alpha, lin_pred = lin_pred)
  r <- r$minimum
  if (alpha_re_bound - abs(r) < 0.1) {
    warning(paste0('bi = ', r, ', alpha_re_bound = ', alpha_re_bound))
  }
  return(r)
}

# ?? what is the formula here, exp(0.3 C_i)/[exp(0.3 C_i)+exp(-b)] - alpha
# keep the average level as alpha by optimize b.
AlphaToBi <- function(b, alpha, lin_pred) {
  exp_lin_pred <- exp(lin_pred)
  r <- abs(mean(exp_lin_pred / (exp_lin_pred + exp(- b))) - alpha)
  return(r)
}

N = sample(x = 10:15, size = nblocks, replace = T)

estimand_glmm<- function(N, betas, num_sub,OR_1=2, OR_2 = 3, OR_C=1.6,OR_WI=1,OR_WC=5,OR_H=1,em=0){
  estimands <- data.frame(alpha = betas,
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
    #' 
    #'
    re_beta <- FromAlphaToRE(alpha = beta, lin_pred = 0.3*C)
    
    pi.beta <- plogis(re_beta + 0.3*C)
    #mean(pi.beta[1:11])
    #mean(pi.beta[12:31])
    # correct
    
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
#source("/Users/congjiang/Documents/simulaINTtnd/Helpfunc_Cor.R")

time = proc.time()

D = 10                                                                         ## Number of simulations
nblocks = 100                                                                      ## Number of clusters in one simulation
r = 1000                                                                      ## Number of subsampling approximation
betas = c(0.3,0.5,0.7)                                                        ## beta values of the GLMM policies

print("[Simulation setting]")
print(paste0("D: ", D))
print(paste0("nblocks: ", nblocks))
print(paste0("betas: ", paste(signif(betas, 4), collapse = ", ")))
print(paste0("num_sub: ", r))

###----------- Estimand computation  ---------------###

# increase num_sub, the results becomes normal 
for(simul.id in 1:D){
  
  ### Cluster sizes ###
  N = sample(x = 20:30, size = nblocks, replace = T)
  
  ### Estimands computation ###
  estimands.list <- lapply(N, estimand_glmm, betas = betas, num_sub = r)
  
  ##################################################
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
  combined_df <- do.call(rbind, lapply(estimands.list, function(df) {
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
  }))
  ############## *********
  
  
  # dim(combined_df)= 11000 and 5 
  # length(alphas)*length(N) (the number of blocks)
  
  #combined_df <- combined_df[order(combined_df$alpha), ]
  #combined_df[432,]
  # mu_1=0, mu_0=0.0348
  
  # Calculate the average log difference for each alphas value
  # exp{ mean[log(df$mu_1) - log(df$mu_0)] }
  # expectation taken over all clusters
  DEavg_log_diff <- combined_df %>%
    group_by(alpha) %>%
    summarise(DE = exp(mean(log_diff, na.rm = TRUE)))
  DEavg_log_diff
  
  # need to remove the one which has -Inf Inf, why does this happen?
  # exp(mean(combined_df$log_diff[combined_df$alpha==0.2],na.rm=T))
  # exp(mean(combined_df$log_diff[combined_df$alpha==0.2][-432]))
  # combined_df$log_diff[combined_df$alpha==0.2][432]
  
  SE1_log_diff_df <- do.call(rbind, lapply(estimands.list, function(df) {
    
    # if df$mu_1=0, then log(0)=-Inf, change to a fairly small number.
    # then when put back to exp, it will close to 0
    #df$log_mu1[which(df$log_mu1==-Inf)]<--
    # get the one which simulated as alpha super close to 0.3.
    # why set as 0.2 does not work???
    standard_log_mu1 <- df$log_mu1[abs(df$alpha - 0.3) < 1e-6] # df$alpha == 0.30
    df$log_mu1_diff <- df$log_mu1 - standard_log_mu1
    return(df)
  }))
  
  
  # #SE: Compute the log differences and combine the results
  # # comparing with policy 0.3
  # SE1_log_diff_df <- do.call(rbind, lapply(estimands.list, function(df) {
  #   
  #   # Identify rows where mu_1 or mu_0 is exactly 0
  #   # rare_case_rec <- which(abs(df$mu_1) == 0 )
  #   # 
  #   # # Remove those rows
  #   # if (length(rare_case_rec) > 0) {
  #   #   df <- df[-rare_case_rec, ]
  #   # }
  #   
  #   
  #   df$log_mu1 <- log(df$mu_1)
  #   
  #   # get the one which simulated as alpha super close to 0.3.
  #   # why set as 0.2 does not work???
  #   standard_log_mu1 <- df$log_mu1[abs(df$alpha - 0.40) < 1e-6] # df$alpha == 0.30
  #   df$log_mu1_diff <- df$log_mu1 - standard_log_mu1
  #   return(df)
  # }))
  
  # Calculate the average log difference for each alpha value
  SEavg_log_mu1_diff <- SE1_log_diff_df %>%
    group_by(alpha) %>%
    summarise(SE1 = exp(mean(log_mu1_diff, na.rm = TRUE)))
  SEavg_log_mu1_diff
  
  SE0_log_diff_df <- do.call(rbind, lapply(estimands.list, function(df) {
    rare_case_rec <- which(abs(df$mu_0) == 0)
    
    # Remove those rows, then not have that alpha's value
    # having dataset all rows have 0 exists.
    # if something equal to 0, then change to a pretty small number 
    if(length(rare_case_rec) > 0 & length(rare_case_rec)<length(df$mu_0) ) {
      df$mu_0[rare_case_rec]<-10^(-6)
      #df <- df[-rare_case_rec, ]
      df$log_mu0 <- log(df$mu_0)
    }else if(length(rare_case_rec)==length(df$mu_0)){
      
      # if mu_1 all =0 across alphas 
      df$mu_0<-10^(-6)
      df$log_mu0 <- log(df$mu_0)
      
    }else{
      df$log_mu0 <- log(df$mu_0)
    }
    
    #df$log_mu0 <- log(df$mu_0)
    standard_log_mu0 <- df$log_mu0[abs(df$alpha - 0.3) < 1e-6]
    df$log_mu0_diff <- df$log_mu0 - standard_log_mu0
    return(df)
  }))
  
  # Calculate the average log difference for each alpha value
  SEavg_log_mu0_diff <- SE0_log_diff_df %>%
    group_by(alpha) %>%
    summarise(SE0 = exp(mean(log_mu0_diff, na.rm = TRUE)))
  SEavg_log_mu0_diff
  
  #OE: Compute the log differences and combine the results
  OE_log_diff_df <- do.call(rbind, lapply(estimands.list, function(df) {
    df$log_mu <- log(df$mu)
    standard_log_mu <- df$log_mu[abs(df$alpha - 0.3) < 1e-6]
    df$log_mu_diff <- df$log_mu - standard_log_mu
    return(df)
  }))
  
  
  # Calculate the average log difference for each alpha value (?? Difference between OE and TE?)
  OEavg_log_mu_diff <- OE_log_diff_df %>%
    group_by(alpha) %>%
    summarise(OE = exp(mean(log_mu_diff, na.rm = TRUE)))
  OEavg_log_mu_diff
  
  #TE: Compute the log differences and combine the results (total effect)
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
  
  ##################################################
  
  # estimands <- data.frame(beta = rep(0, length(betas)),
  #                         mu = 0, mu_1 = 0, mu_0 = 0)
  # 
  # for(i in 1:nblocks){
  #   estimands = estimands + estimands.list[[i]]
  # }
  # 
  # estimands = estimands / nblocks
  # 
  # ### Causal effects computation using beta==0.5 as the standard ###
  # standard = estimands %>% filter(beta == 0.5)
  # 
  # estimands = estimands %>% mutate(de  =  mu_1/mu_0,
  #                                  se_1 =  mu_1/standard$mu_1,
  #                                  se_0 =  mu_0/standard$mu_0,
  #                                  oe  =  mu/standard$mu,
  #                                  te  =  mu_1/standard$mu_0)
  # 
  # ### Save output ###
  # saveRDS(estimands, file = paste0("estimandGLMM_id", simul.id,".rds"))
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
