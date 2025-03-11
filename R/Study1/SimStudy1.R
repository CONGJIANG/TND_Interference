library(data.table)

source("TNDintIPWestimator.R")
source("function_class.R")

datagen_int <- function(rangeN = 400:500, nblocks = 1000,
                        OR_0 = 3, OR_1 = 2, OR_2 = 3, OR_C = 1.6, OR_WI = 1, OR_WC = 5, OR_H = 1, OR_H_g_V = 1,
                        em = 0, Infec_a = 0.25, Infec_b = 3, Infec_c = 1.5, 
                        treat_type = "GLMM", parameter_vec) {
  #' treat_type="GLMM"
  #' parameter_vec=c(alpha,beta_0,beta_1,include_RE=F, sigma=sigma)
  #' alpha: cluster level proportion of treated 
  #' include_RE=F 
  #' beta_0: NA, need to solve, specify beta_1 
  #' include_RE=T
  #' specify beta_0, beta_1, alpha=NA, sigma (random effect's standard error)
  #' rangeN: size of each cluster, length = the number of clusters
  #' nblocks: the number of blocks
  # Preallocate lists for three datasets
  data.list_TND <- vector("list", nblocks)
  data.list_full <- vector("list", nblocks)
  data.list_SRS <- vector("list", nblocks)
  
  N <- sample(x = rangeN, size = nblocks, replace = TRUE)
  vaxpos <- rnorm(nblocks, mean = 0, sd = parameter_vec["sigma"])
  
  for (i in 1:nblocks) {
    C1 <- runif(n = N[i], 1, 2)
    C2 <- rnorm(n = N[i], 0.6, 0.5)
    U1 <- rbinom(n = N[i], size = 1, prob = 0.5)
    
    if (treat_type == "GLMM") {
      if (!parameter_vec["include_RE"]) {
        re_beta <- FromAlphaToRE(alpha = parameter_vec["alpha"],  lin_pred = parameter_vec["beta_1"] * C1 + parameter_vec["beta_2"] * C2)
        p_trt <- plogis(re_beta + parameter_vec["beta_1"] * C1 + parameter_vec["beta_2"] * C2)
      } else {
        p_trt <- plogis(parameter_vec["beta_0"] + parameter_vec["beta_1"] * C1 + parameter_vec["beta_2"] * C2 + vaxpos[i])
      }
    } else if (treat_type == "Type-B") {
      p_trt <- rep(parameter_vec, N[i])
    }
    
    V <- rbinom(prob = p_trt, size = 1, n = N[i])
    g.V <- (sum(V) - V) / (N[i] - 1)
    Infec <- rbinom(prob = plogis(Infec_a * C1 - Infec_b + Infec_c * U1), size = 1, n = N[i])
    which(Infec == 1)
    logit_infprob <- -OR_0 + C1 - log(OR_1) * V - log(OR_2) * g.V - log(OR_C) * V * g.V + em * V * C1 + log(2) * C2 - 2 * U1
    Infprob <- plogis(logit_infprob)
    Infprob[which(Infec == 1)] <- 0
    Infec_COVID <- rbinom(prob = Infprob, size = 1, n = N[i])
    
    W1 <- W2 <- W <- rep(0, N[i])
    W1[Infec == 1] <- rbinom(prob = plogis(-2 + 0.5 * C1[Infec == 1] - log(OR_WI) * V[Infec == 1] - 0.5 * U1[Infec == 1]), size = 1, n = sum(Infec == 1))
    pw2 <- plogis(-1.5 + C1[Infec_COVID == 1] - log(OR_WC) * V[Infec_COVID == 1] 
                  - U1[Infec_COVID == 1] + 0.5 * C2[Infec_COVID == 1] * (1 - V[Infec_COVID == 1]))
    W2[Infec_COVID == 1] <- rbinom(prob = pw2, size = 1, n = sum(Infec_COVID))
    W <- (W1 | W2)
    
    H <- rep(0, N[i])
    Hprob <- plogis(0.5 * C1[W == 1] - log(OR_H) * V[W == 1] - 0.5 * U1[W == 1] + log(OR_H_g_V) * g.V[W == 1])
    H[W == 1] <- rbinom(prob = Hprob, size = 1, n = sum(W == 1))
    
    R <- if (sum(H == 1) > 0) sample((1:N[i])[H == 1], min(rangeN[1], sum(H == 1))) else integer(0)
    R_SRS <- sample(seq_len(N[i]), size = floor(0.03 * rangeN[1]), replace = FALSE)
    
    
    data.list_TND[[i]] <- as.data.frame(cbind(Y=Infec_COVID,V=V,C1=C1,C2=C2,block=i,f_m=g.V))[R,]
    data.list_full[[i]] <- as.data.frame(cbind(Y = Infec_COVID, Infec = Infec, H=H, W=W, V=V, C1=C1,C2=C2,block=i,f_m=g.V, p_trt=p_trt))
    data.list_SRS[[i]] <- data.list_full[[i]][R_SRS, ]
  }   
  
  data_TND<-dplyr::bind_rows(data.list_TND)
  data_full<-dplyr::bind_rows(data.list_full)
  data_SRS<-dplyr::bind_rows(data.list_SRS)
  data<-list(data_TND=data_TND,data_full=data_full, data_SRS =data_SRS)
  
  # Return all three datasets
  return(data)
}





MC_sim <- function(r, parameter_vec) {
  output_file <- paste0("/n/home09/c55jiang/TND_Interference/Mar5_7", ".txt")
  
  for (i in 1:r) {
    data <- datagen_int(rangeN<-5000:8000, nblocks = 2000,parameter_vec = parameter_vec)
    
    datTND<-data$data_TND
    datfull <-data$data_full
    datSRS <-data$data_SRS
    
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
    phi_hat_tnd <- PS_model(datTND, glm_form = glm_form, method = 'TND_IPW')
    # ***tnd data, ori_IPW
    phi_hat_org <- PS_model(datTND, glm_form = glm_form, method = 'ORG_IPW')
    # ***srs data, ori_IPW
    phi_hat_srs <- PS_model(datSRS, glm_form = glm_form, method = 'ORG_IPW')
    # ***full data, ori_IPW
    phi_hatfull <- PS_model(datfull, glm_form = glm_form, method = 'ORG_IPW')
    
    # Save results
    results <- c(i, c(phi_hat_tnd$coefs, phi_hat_tnd$re_var), 
                 c(phi_hat_srs$coefs, phi_hat_srs$re_var), 
                 c(phi_hat_org$coefs, phi_hat_org$re_var), 
                 c(phi_hatfull$coefs, phi_hatfull$re_var) )
    
    write(results, file = output_file, ncolumns = length(results), append = TRUE)
  }
}


# Example usage
r <- 1000   # Number of replicates
parameter_vec <- c(alpha=NA,beta_0=-0.6,beta_1=0.3,beta_2=0.5,include_RE=T, sigma=0.5)
sim.res <- MC_sim(r = r, parameter_vec = parameter_vec)


