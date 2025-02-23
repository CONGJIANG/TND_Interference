datagen_int<-function(rangeN = 400:500, nblocks=1000,
                      OR_1=2, OR_2 =3, OR_C=1.6,OR_WI=1,OR_WC=5,OR_H=1,em=0, return_full=F
                      ,treat_type="GLMM",parameter_vec){
  #' treat_type="GLMM"
  #' parameter_vec=c(alpha,beta_0,beta_1,include_RE=F, sigma=sigma)
  #' alpha: cluster level proportion of treated 
  #' include_RE=F 
  #' beta_0: NA, need to solve, specify beta_1 
  #' include_RE=T
  #' specify beta_0, beta_1, alpha=NA, sigma (random effect's standard error)
  #' rangeN: size of each cluster, length = the number of clusters
  #' nblocks: the number of blocks
  
  if (return_full == FALSE) {
    data.list <- list()
  } else if(return_full == T){
    data.list <- list()
  }else if(return_full == "BOTH"){
    data.list_select <- list()
    data.list_full <- list()
  }
  
  N <- sample(x = rangeN, size = nblocks, replace = TRUE)
  #vaxpos <- rtruncnorm(nblocks, -1, 1, mean = 0, sd = 1)
  # ?? should we have it in our case or not.
  # rangeN each 
  vaxpos <- rtruncnorm(nblocks, -1, 1, mean = 0, sd = parameter_vec["sigma"] )
  # in total nblocks blocks, each block has N[i] people. 
  for(i in 1:nblocks){
    # Step1. Data generation
    #generate data (WITH clustering for U2)
    C<-runif(n=N[i], 1, 2)
    U1<-rbinom(n=N[i],size=1,prob=0.5) #affects both
    U2<-rbinom(n=N[i],size=1,prob=plogis(0.4+0.2*vaxpos[i])) #affects covid #shifted upwards for vaxpos blocks #it's like reckless behaviour specific to COVID
    
    # Step2. Treatment model
    if(treat_type=="GLMM"){
      
      if(parameter_vec["include_RE"]==F){
        # need psi to make the average to be alpha within each cluster. 
        # If having random effect, do the optimization, still get the re_beta, which 
        re_beta <- FromAlphaToRE(alpha = parameter_vec["alpha"], lin_pred = parameter_vec["beta_1"]*C)
        p_trt <- plogis(re_beta + parameter_vec["beta_1"]*C)
        # if having the random effect, the optimal solution is not be the one for 
        # average alpha 
        #p_trt <- plogis(re_beta + parameter_vec["beta_1"]*C + parameter_vec["beta_2"]*vaxpos[i])
      }else{
        ## 
        # vaxpos <- rtruncnorm(nblocks, -1, 1, mean = 0, sd = parameter_vec["sigma"] )
        p_trt <- plogis(parameter_vec["beta_0"] + parameter_vec["beta_1"]*C+0.5*vaxpos[i])
      }
      
    }else if(treat_type=="Type-B"){
      p_trt<-rep(parameter_vec,N[i])
    }
    
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
    # **** realized: set to log(1)*g.V[W==1], this will be more realistic.
    Hprob<-plogis(0.5*C[W==1]-log(OR_H)*V[W==1]-0.5*U1[W==1]+g.V[W==1])
    H[W==1]<-rbinom(prob=Hprob,size=1,n=sum(W==1))
    #selection on outcome for testing (does not condition on infectious status, just being in the hospital)
    R <- as.vector(na.omit(sample((1:N[i])[H == 1], min(rangeN[1], sum(H == 1)))))
    # rho^0_i = Hprob*plogis(-1.5+1*C[Infec_COVID==1]-log(OR_WC)*V[Infec_COVID==1]-1*U1[Infec_COVID==1]+0.5*U2[Infec_COVID==1]*(1-V[Infec_COVID==1]))*Infprob
    
    if (return_full == FALSE) {
      data.list[[i]] <- as.data.frame(cbind(Y=Infec_COVID,V=V,C=C,block=i,f_m=g.V))[R,]
    } else if(return_full == T){
      data.list[[i]] <- as.data.frame(cbind(Y = Infec_COVID, Infec = Infec, H=H, W=W, V=V, C=C,block=i,f_m=g.V))
    }else if(return_full == "BOTH"){
      data.list_select[[i]] <- as.data.frame(cbind(Y=Infec_COVID,V=V,C=C,block=i,f_m=g.V))[R,]
      data.list_full[[i]] <- as.data.frame(cbind(Y = Infec_COVID, Infec = Infec, H=H, W=W, V=V, C=C,block=i,f_m=g.V, p_trt=p_trt))
    }
  }
  
  if (return_full == FALSE) {
    data<-dplyr::bind_rows(data.list)
  } else if(return_full == T){
    data<-dplyr::bind_rows(data.list)
  }else if(return_full == "BOTH"){
    data_select<-dplyr::bind_rows(data.list_select)
    data_full<-dplyr::bind_rows(data.list_full)
    data<-list(data_select=data_select,data_full=data_full)
  }
  
  return(data)
}



datagen_int_old1<-function(rangeN = 400:500, nblocks=1000,
                      OR_1=2, OR_2 =3, OR_C=1.6,OR_WI=1,OR_WC=5,OR_H=1,em=0, return_full=F
                      ,treat_type="GLMM",parameter_vec){
  
  if (return_full == FALSE) {
    data.list <- list()
  } else if(return_full == T){
    data.list <- list()
  }else if(return_full == "BOTH"){
    data.list_select <- list()
    data.list_full <- list()
  }
  
  N <- sample(x = rangeN, size = nblocks, replace = TRUE)
  #vaxpos <- rtruncnorm(nblocks, -1, 1, mean = 0, sd = 1)
  # ?? should we have it in our case or not.
  vaxpos <- rtruncnorm(nblocks, -1, 1, mean = 0, sd = 0.1)
  # in total nblocks blocks, each block has N[i] people. 
  for(i in 1:nblocks){
    # Step1. Data generation
    #generate data (WITH clustering for U2)
    C<-runif(n=N[i], 1, 2)
    U1<-rbinom(n=N[i],size=1,prob=0.5) #affects both
    U2<-rbinom(n=N[i],size=1,prob=plogis(0.4+0.2*vaxpos[i])) #affects covid #shifted upwards for vaxpos blocks #it's like reckless behaviour specific to COVID
    
    # Step2. Treatment model
    if(treat_type=="GLMM"){
      
      # need psi to make the average to be alpha within each cluster. 
      # If having random effect, do the optimization, still get the re_beta, which 
      # will make the average level as alpha?
      re_beta <- FromAlphaToRE(alpha = parameter_vec["alpha"], lin_pred = parameter_vec["beta_1"]*C)
      p_trt <- plogis(re_beta + parameter_vec["beta_1"]*C)
      # if having the random effect, the optimal solution is not be the one for 
      # average alpha 
      #p_trt <- plogis(re_beta + parameter_vec["beta_1"]*C + parameter_vec["beta_2"]*vaxpos[i])
      
      
    }else if(treat_type=="Type-B"){
      p_trt<-rep(parameter_vec,N[i])
    }
    
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
    # rho^0_i = Hprob*plogis(-1.5+1*C[Infec_COVID==1]-log(OR_WC)*V[Infec_COVID==1]-1*U1[Infec_COVID==1]+0.5*U2[Infec_COVID==1]*(1-V[Infec_COVID==1]))*Infprob
    
    if (return_full == FALSE) {
      data.list[[i]] <- as.data.frame(cbind(Y=Infec_COVID,V=V,C=C,block=i,f_m=g.V))[R,]
    } else if(return_full == T){
      data.list[[i]] <- as.data.frame(cbind(Infec_COVID = Infec_COVID, Infec = Infec, H=H, W=W, V=V, C=C,block=i,f_m=g.V))
    }else if(return_full == "BOTH"){
      data.list_select[[i]] <- as.data.frame(cbind(Y=Infec_COVID,V=V,C=C,block=i,f_m=g.V))[R,]
      data.list_full[[i]] <- as.data.frame(cbind(Infec_COVID = Infec_COVID, Infec = Infec, H=H, W=W, V=V, C=C,block=i,f_m=g.V, p_trt=p_trt))
    }
  }
  
  if (return_full == FALSE) {
    data<-dplyr::bind_rows(data.list)
  } else if(return_full == T){
    data<-dplyr::bind_rows(data.list)
  }else if(return_full == "BOTH"){
    data_select<-dplyr::bind_rows(data.list_select)
    data_full<-dplyr::bind_rows(data.list_full)
    data<-list(data_select=data_select,data_full=data_full)
  }
  
  return(data)
}


estimand_glmm_one_iter_old1<- function(N, betas,datFull){
  estimands <- data.frame(alpha = betas,
                          mu = 0, mu_1 = 0, mu_0 = 0)
  
  for(i in seq_len(length(betas))){
    
    # the datFull$p_trt will be correct p_trt? 
    mu   = mean(datFull$H*datFull$Infec_COVID, na.rm = T)
    mu_1 = mean(datFull$H*datFull$Infec_COVID * datFull$V / (datFull$p_trt), na.rm = T)
    mu_0 = mean(datFull$H*datFull$Infec_COVID * (1-datFull$V) / (1-datFull$p_trt), na.rm = T)
    
    estimands[i, ] = c(betas, mu, mu_1, mu_0)
    
  }
  return(estimands)
}

# data generation
datagen_int_old1<-function(rangeN = 400:500, nblocks=1000,
                      OR_1=2, OR_2 =3, OR_C=1.6,OR_WI=1,OR_WC=5,OR_H=1,em=0, return_full=F
                      ,treat_type="GLMM",parameter_vec){
  
  if (return_full == FALSE) {
    data.list <- list()
  } else if(return_full == T){
    data.list <- list()
  }else if(return_full == "BOTH"){
    data.list_select <- list()
    data.list_full <- list()
  }
  
  N <- sample(x = rangeN, size = nblocks, replace = TRUE)
  #vaxpos <- rtruncnorm(nblocks, -1, 1, mean = 0, sd = 1)
  vaxpos <- rtruncnorm(nblocks, -1, 1, mean = 0, sd = 0.1)
  # in total nblocks blocks, each block has N[i] people. 
  for(i in 1:nblocks){
    # Step1. Data generation
    #generate data (WITH clustering for U2)
    C<-runif(n=N[i], 1, 2)
    U1<-rbinom(n=N[i],size=1,prob=0.5) #affects both
    U2<-rbinom(n=N[i],size=1,prob=plogis(0.4+0.2*vaxpos[i])) #affects covid #shifted upwards for vaxpos blocks #it's like reckless behaviour specific to COVID
    
    # Step2. Treatment model
    if(treat_type=="GLMM"){
      
      # need psi to make the average to be alpha within each cluster. 
      # If having random effect, do the optimization, still get the re_beta, which 
      # will make the average level as alpha?
      re_beta <- FromAlphaToRE(alpha = parameter_vec["alpha"], lin_pred = parameter_vec["beta_1"]*C)
      p_trt <- plogis(re_beta + parameter_vec["beta_1"]*C)
      # if having the random effect, the optimal solution is not be the one for 
      # average alpha 
      #p_trt <- plogis(re_beta + parameter_vec["beta_1"]*C + parameter_vec["beta_2"]*vaxpos[i])
      
    }else if(treat_type=="Type-B"){
      p_trt<-rep(parameter_vec,N[i])
    }
    
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
    
    if (return_full == FALSE) {
      data.list[[i]] <- as.data.frame(cbind(Y=Infec_COVID,V=V,C=C,block=i,f_m=g.V))[R,]
    } else if(return_full == T){
      data.list[[i]] <- as.data.frame(cbind(Infec_COVID = Infec_COVID, Infec = Infec, H=H, W=W, V=V, C=C,block=i,f_m=g.V))
    }else if(return_full == "BOTH"){
      data.list_select[[i]] <- as.data.frame(cbind(Y=Infec_COVID,V=V,C=C,block=i,f_m=g.V))[R,]
      data.list_full[[i]] <- as.data.frame(cbind(Infec_COVID = Infec_COVID, Infec = Infec, H=H, W=W, V=V, C=C,block=i,f_m=g.V))
    }
  }
  
  if (return_full == FALSE) {
    data<-dplyr::bind_rows(data.list)
  } else if(return_full == T){
    data<-dplyr::bind_rows(data.list)
  }else if(return_full == "BOTH"){
    data_select<-dplyr::bind_rows(data.list_select)
    data_full<-dplyr::bind_rows(data.list_full)
    data<-list(data_select=data_select,data_full=data_full)
  }
  
  return(data)
}


# data generation
datagen_int_old1<-function(rangeN = 400:500, nblocks=1000,
                      OR_1=2, OR_2 =3, OR_C=1.6,OR_WI=1,OR_WC=5,OR_H=1,em=0, return_full=F
                      ,treat_type="GLMM",parameter_vec){
  
  if (return_full == FALSE) {
    data.list <- list()
  } else if(return_full == T){
    data.list <- list()
  }else if(return_full == "BOTH"){
    data.list_select <- list()
    data.list_full <- list()
  }

  N <- sample(x = rangeN, size = nblocks, replace = TRUE)
  vaxpos <- rtruncnorm(nblocks, -1, 1, mean = 0, sd = 1)
  # in total nblocks blocks, each block has N[i] people. 
  for(i in 1:nblocks){
    # Step1. Data generation
    #generate data (WITH clustering for U2)
    C<-runif(n=N[i], 1, 2)
    U1<-rbinom(n=N[i],size=1,prob=0.5) #affects both
    U2<-rbinom(n=N[i],size=1,prob=plogis(0.4+0.2*vaxpos[i])) #affects covid #shifted upwards for vaxpos blocks #it's like reckless behaviour specific to COVID
    
    # Step2. Treatment model
    if(treat_type=="GLMM"){
      p_trt <- plogis(parameter_vec[1] + parameter_vec[2]*C + parameter_vec[3]*vaxpos[i])
    }else if(treat_type=="Type-B"){
      p_trt<-rep(parameter_vec,N[i])
    }
    
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
    
    if (return_full == FALSE) {
      data.list[[i]] <- as.data.frame(cbind(Y=Infec_COVID,V=V,C=C,block=i,f_m=g.V))[R,]
    } else if(return_full == T){
      data.list[[i]] <- as.data.frame(cbind(Infec_COVID = Infec_COVID, Infec = Infec, H=H, W=W, V=V, C=C,block=i,f_m=g.V))
    }else if(return_full == "BOTH"){
      data.list_select[[i]] <- as.data.frame(cbind(Y=Infec_COVID,V=V,C=C,block=i,f_m=g.V))[R,]
      data.list_full[[i]] <- as.data.frame(cbind(Infec_COVID = Infec_COVID, Infec = Infec, H=H, W=W, V=V, C=C,block=i,f_m=g.V))
    }
  }
  
  if (return_full == FALSE) {
    data<-dplyr::bind_rows(data.list)
  } else if(return_full == T){
    data<-dplyr::bind_rows(data.list)
  }else if(return_full == "BOTH"){
    data_select<-dplyr::bind_rows(data.list_select)
    data_full<-dplyr::bind_rows(data.list_full)
    data<-list(data_select=data_select,data_full=data_full)
  }
  
  return(data)
}


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


estimand_glmm_one_iter<-function(N, betas,datFull){
  
  nblocks<-length(unique(datFull$block))
  
  estimands <- data.frame(mu = rep(0,nblocks ), mu_1 = rep(0, nblocks), mu_0 = rep(0,nblocks) )
  
  # for each block, compute mu,mu_1,mu_0
  # obtain 
  # block 1, mu, mu_1, mu_0
  # block 2, mu, mu_1, mu_0
  
  for(i in 1:nblocks){
    
    block_ind<-which(datFull$block==i)
    # the datFull$p_trt will be correct p_trt? 
    estimands$mu[i] = mean(datFull$H[block_ind]*datFull$Infec_COVID[block_ind], na.rm = T)
    estimands$mu_1[i] = mean(datFull$H[block_ind]*datFull$Infec_COVID[block_ind] * datFull$V[block_ind] / (datFull$p_trt[block_ind]), na.rm = T)
    estimands$mu_0[i] = mean(datFull$H[block_ind]*datFull$Infec_COVID[block_ind] * (1-datFull$V[block_ind]) / (1-datFull$p_trt[block_ind]), na.rm = T)
  }
  return(estimands)
}


