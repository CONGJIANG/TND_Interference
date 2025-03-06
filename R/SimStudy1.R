MC_sim <- function(r, parameter_vec = c(alpha=NA,beta_0=-0.25,beta_1=0.3,include_RE=T, sigma=0.3)) {
  type <- match.arg(type)
  output_file <- paste0("/n/home09/c55jiang/mediation/FEB19_B5", type, "_n_", n_obs, "rate", censor_rate, ".txt")
  
  for (i in 1:r) {
    
    datTND<-datagen_int(rangeN<-100:120, nblocks = 1000,parameter_vec = parameter_vec)
    datfull <-datagen_int(rangeN<-100:120, nblocks = 1000, return_full=T, parameter_vec=parameter_vec)
    
    
    cov_names <- c('C')
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
    # ***full data, ori_IPW
    phi_hat_org <- PS_model(datTND, glm_form = glm_form, method = 'ORG_IPW')
    # ***full data, ori_IPW
    phi_hatfull <- PS_model(datfull, glm_form = glm_form, method = 'ORG_IPW')
    
    # Save results
    results <- c(i, c(phi_hat_tnd$coefs, phi_hat_tnd$re_var), 
                 c(phi_hat_org$coefs, phi_hat_org$re_var), 
                 c(phi_hatfull$coefs, phi_hatfull$re_var) )
    
    write(results, file = output_file, ncolumns = length(results), append = TRUE)
  }
  
  return(beta_m)  # Returning the final beta_m if needed
}


# Example usage
r <- 1000   # Number of replicates
sim.res <- MC_sim(r = r)