#' TND Partial Interference DATA generation
#' 
#'
#'
#'
library(truncnorm)
library(Matrix)
library(lme4)
#*** key for lme4 works.
#utils::install.packages("lme4", type = "source")

source("TNDintIPWestimator.R")

datagen_int<-function(rangeN = 400:500, nblocks=1000,
                      OR_1=2, OR_2 =3, OR_C=1.6,OR_WI=1,OR_WC=5,OR_H=1,em=0, return_full=F){
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
    
    if (return_full == FALSE) {
      data.list[[i]] <- as.data.frame(cbind(Y=Infec_COVID,V=V,C=C,block=i,f_m=g.V))[R,]
    } else {
      data.list[[i]] <- as.data.frame(cbind(Infec_COVID = Infec_COVID, Infec = Infec, H=H, W=W, V=V, C=C,block=i,f_m=g.V))
    }
  }
  return(dplyr::bind_rows(data.list))
}
set.seed(2025)
(datTND<-datagen_int(OR_1=10, OR_2 =10,nblocks=1000))

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
  
  GM_IE_sumv1 <- list(); GM_IE_sumv2 <- list()
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
    GM_IE_sumv1[[alpha2]] <- temp
  }
  
  GM_IE_sumv1 <- bind_rows(GM_IE_sumv1);  GM_IE_sumv0 <- bind_rows(GM_IE_sumv0)
  
  return(list(GM_DE = GM_DE_summary, GM_IEv1 = GM_IE_sumv1, GM_IEv0 = GM_IE_sumv0))
}

# Run MC replications and store results
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


