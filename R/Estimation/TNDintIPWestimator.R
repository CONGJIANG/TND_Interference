# Clear the entire environment
rm(list = ls())
#' Propensity score function to get phi_hat
#'
#' @param datTND TND Data frame including treatment, outcome, and covariates.
#' @param glm_form PS model formula 
#' @param method Character, either 'ORG_IPW' or 'TND_IPW'. If 'ORG_IPW' is specified, then the
#' PS is estimated by using the whole TND data (i.e., original partial interference methods). If the method
#' is set equal to 'TND_IPW', use the proposed TND IPW estimator for PS, using only controls. 
#' 
#' @return phi_hat A list with two elements. The first one is a vector of
#' coefficients of the PS, and the second one is the random effect variance.
#' If the second element is 0, the propensity score excludes random effects.
#'
#' @export

PS_model <- function(datTND, glm_form, method = c('ORG_IPW', 'TND_IPW')){
  if (method == 'TND_IPW') {
    glmer_fit_col <- lme4::glmer(
      data = datTND[datTND$Y==0, ],
      formula = as.formula(glm_form),
      family = stats::binomial,
      nAGQ = 2
    )
    phi_hat <- list(coefs = summary(glmer_fit_col)$coef[, 1],
                    re_var = as.numeric(summary(glmer_fit_col)$varcor))
  }
  
  if (method == 'ORG_IPW') {
    glmer_fit <- lme4::glmer(
      data = datTND,
      formula = as.formula(glm_form),
      family = stats::binomial,
      nAGQ = 2
    )
    phi_hat <- list(coefs = summary(glmer_fit)$coef[, 1],
                    re_var = as.numeric(summary(glmer_fit)$varcor))
  }
  return(phi_hat)
}


#' Propensity score function to get phi_hat
#'
#' @param datTND TND Data frame including treatment, outcome, and covariates.
#' @param gamma_form PS model formula 
#' @param method Character, either 'ORG_IPW' or 'TND_IPW'. If 'ORG_IPW' is specified, then the
#' PS is estimated by using the whole TND data (i.e., original partial interference methods). If the method
#' is set equal to 'TND_IPW', use the proposed TND IPW estimator for PS, using only controls. 
#' 
#' @return gamma_numer The coefficients of the ps model in the numerator.
#'
#' @export
Policy_coef <- function(datTND, gamma_form, method = c('ORG_IPW', 'TND_IPW')){
  if (method == 'TND_IPW') {
    gamma_glmod <- lme4::glmer(as.formula(gamma_form), data = datTND[datTND$Y==0, ], family = 'binomial',
                               control = glmerControl(optimizer = "bobyqa",
                                                      optCtrl = list(maxfun = 2e5)))
    gamma_numer <- summary(gamma_glmod)$coef[, 1]
    rm(list = c('gamma_form', 'gamma_glmod'))
  }
  
  if (method == 'ORG_IPW') {
    gamma_glmod <- lme4::glmer(as.formula(gamma_form), data = datTND, family = 'binomial',
                               control = glmerControl(optimizer = "bobyqa",
                                                      optCtrl = list(maxfun = 2e5)))
    gamma_numer <- summary(gamma_glmod)$coef[, 1]
    rm(list = c('gamma_form', 'gamma_glmod'))
  }
  return(gamma_numer)
}


###########################################
### GroupIPW needs
### (1)FromAlphaToRE and AlphaToBi
### (2)CalcNumerator
### (3)Denominator
###########################################
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

CalcNumerator <- function(Ai_j, Xi_j, gamma_numer, alpha, re_alpha,
                          include_alpha = TRUE) {
  
  gamma_numer <- matrix(gamma_numer, nrow = length(gamma_numer), ncol = 1)
  
  lin_pred <- cbind(1, as.matrix(Xi_j)) %*% gamma_numer
  lin_pred <- lin_pred + re_alpha
  probs <- expit(lin_pred)
  
  if (include_alpha) {
    r <- (probs / alpha) ^ Ai_j * ((1 - probs) / (1 - alpha)) ^ (1 - Ai_j)
  } else {
    r <- probs ^ Ai_j * (1 - probs) ^ (1 - Ai_j)
  }
  
  return(list(prob = prod(r), re_alpha = re_alpha))
}

expit <- function(x) {
  return(exp(x) / (1 + exp(x)))
}




#' Calculating the denominator of the group estimator.
#' 
#' @param A The treatment vector for units in the group.
#' @param X The covariates for units in the group.
#' @param phi_hat A list with two elements. The first one is a vector of
#' coefficients of the ps, and the second one is the random effect variance.
#' If the second element is 0, the propensity score excludes random effects.
#' @param alpha This value of alpha is used to stabilize calculations. If
#' include_alpha is set to TRUE, alpha needs to coincide with the alpha used
#' for calculating the numerator.
#' @param integral_bound If the propensity score includes a random effect, the
#' integral in the denominator is calculated over the normal distribution of
#' the random effects from - integral_bound to integral_bound.
#' @param include_alpha If include_alpha is set to true, the probabilities in
#' the denominator are divided by the specified value of alpha to stabilize the
#' integral calculation.
#' 
#' @export
Denominator <- function(A, X, phi_hat, alpha = NULL, integral_bound = 10,
                        include_alpha = TRUE) {
  
  integral_bound <- abs(integral_bound)
  if (is.null(alpha) & include_alpha) {
    stop('No alpha provided.')
  }
  
  X <- as.matrix(cbind(1, X))
  re_sd <- sqrt(phi_hat[[2]])
  phi_hat[[1]] <- matrix(phi_hat[[1]], nrow = length(phi_hat[[1]]), ncol = 1)
  
  # Creating the function that we will integrate over.
  f_int <- function(b) {
    r <- 1
    lin_pred <- X %*% phi_hat[[1]]  # Includes intercept.
    for (ii in 1:length(A)) {
      prob_trt <- expit(lin_pred[ii] + b)
      if (include_alpha) {
        success_weight <- prob_trt / alpha
        failure_weight <- (1 - prob_trt) / (1 - alpha)
      } else {
        success_weight <- prob_trt
        failure_weight <- 1 - prob_trt
      }
      r <- r * success_weight ^ A[ii] * failure_weight ^ (1 - A[ii])
    }
    if (re_sd > 0) { # If re_sd = 0, there is no random effect in the ps model.
      r <- r * dnorm(b, mean = 0, sd = re_sd)
    }
    return(r)
  }
  
  if (re_sd > 0) {
    ans <- integrate(f_int, lower = - integral_bound * re_sd,
                     upper = integral_bound * re_sd)
  } else {
    ans <- list(value = f_int(0))
  }
  
  return(ans)
}


#' Estimating the group average potential outcome
#' 
#' IPW estimator of the group average potential outcome.
#'
#' @param dta Data frame including treatment, outcome and covariates.
#' @param cov_cols The indices including the covariates of the ps model.
#' @param phi_hat A list with two elements. The first one is a vector of
#' coefficients of the ps, and the second one is the random effect variance.
#' If the second element is 0, the propensity score excludes random effects.
#' @param gamma_numer The coefficients of the ps model in the numerator.
#' If left NULL and estimand is 1, the coefficients in phi_hat will be used
#' instead.
#' @param alpha The values of alpha for which we want to estimate the group
#' average potential outcome.
#' @param neigh_ind List. i^{th} element is a vector with the row indices of
#' dta that are in cluster i. Can be left NULL.
#' @param trt_col If the treatment is not named 'A' in dta, specify the
#' treatment column index.
#' @param out_col If the outcome is not named 'Y', specify the outcome column
#' index.
#' @param alpha_re_bound The lower and upper end of the values for bi we will
#' look at. Defaults to 10, meaning we will look between - 10 and 10.
#' @param integral_bound The number of standard deviations of the random effect
#' that will be used as the lower and upper limit.
#' @param keep_re_alpha Logical. If set to TRUE the "random" effect that makes
#' the average probability of treatment equal to alpha will be returned along
#' with the estimated group average potential outcome.
#' @param estimand Character, either 'GLMM' or 'TypeB'. If 'GLMM' is specified, then the
#' estimand with numerator depending on covariates is estimated. If estimand
#' is set equal to 'TypeB', the numerator considered is the product of Bernoulli.
#' @param verbose Whether printing of progress is wanted. Defaults to TRUE.
#' 
#' @export
GroupIPW <- function(dta, cov_cols, phi_hat, gamma_numer = NULL, alpha,
                     neigh_ind = NULL, trt_col = NULL, out_col = NULL, 
                     alpha_re_bound = 10, integral_bound = 10,
                     keep_re_alpha = FALSE, estimand = c('GLMM', 'TypeB'),
                     verbose = TRUE) {
  
  #estimand <- match.arg(estimand)
  estimand <- match.arg(estimand, choices=c("GLMM",'TypeB'))
  integral_bound <- abs(integral_bound)
  alpha_re_bound <- abs(alpha_re_bound)
  phi_hat[[1]] <- matrix(phi_hat[[1]], ncol = 1)
  dta <- as.data.frame(na.omit(dta))
  
  # We only return the ksi's if we are estimating estimand 'GLMM'.
  keep_re_alpha <- keep_re_alpha & (estimand == 'GLMM')

  # Specifyling neigh_ind will avoid re-running the following lines.
  if (is.null(neigh_ind)) {
    neigh_ind <- sapply(1 : max(dta$block), function(x) which(dta$block == x))
  }
  # removing there is no units in the block (removing because all 0 or 1 are removed)
  neigh_ind<-neigh_ind[-which(unlist(lapply(neigh_ind, function(x){length(x)==0})))]
  n_neigh <- length(neigh_ind)
  
  yhat_group <- array(NA, dim = c(n_neigh, 2, length(alpha)))
  dimnames(yhat_group) <- list(neigh = 1:n_neigh, trt = c(0, 1), alpha = alpha)
  
  # Names of treatment and outcome column.
  if (!is.null(trt_col)) {
    names(dta)[trt_col] <- 'V'
  }
  if (!is.null(out_col)) {
    names(dta)[out_col] <- 'Y'
  }
  if (is.null(gamma_numer)) {
    gamma_numer <- matrix(phi_hat[[1]], ncol = 1)
  }
  
  # If we want to return the ksis that make average propensity alpha.
  if (keep_re_alpha) {
    re_alphas <- matrix(NA, nrow = n_neigh, ncol = length(alpha))
    dimnames(re_alphas) <- list(neigh = 1 : n_neigh, alpha = alpha)
  }
  
  for (aa in 1 : length(alpha)) {
    if (verbose) {
      print(paste('alpha =', round(alpha[aa], 3)))
    }
    curr_alpha <- alpha[[aa]]
    
    for (nn in 1 : n_neigh) {
      
      # For estimand 'GLMM', we need to calculate numerator depending on covariates.
      if (estimand == 'GLMM') {
        
        # Calculating the random effect that gives alpha.
        Xi <- dta[neigh_ind[[nn]], cov_cols]
        lin_pred <- cbind(1, as.matrix(Xi)) %*% gamma_numer
        re_alpha <- FromAlphaToRE(alpha = curr_alpha, lin_pred = lin_pred,
                                  alpha_re_bound = alpha_re_bound)
        
        # Keeping the intercept that makes cluster average propensity alpha.
        if (keep_re_alpha) {
          re_alphas[nn, aa] <- re_alpha
        }
        
      }
      
      for (curr_it in c(0, 1)) {
        bern_prob <- curr_alpha ^ curr_it * (1 - curr_alpha) ^ (1 - curr_it)
        prob_ind <- list(prob = 1)  # For estimand 'TypeB'.
        y_curr <- 0
        
        for (ind in neigh_ind[[nn]]) {
          if (dta$V[ind] == curr_it) {
            
            if (estimand == 'GLMM') {
              wh_others <- setdiff(neigh_ind[[nn]], ind)
              Ai_j <- dta$V[wh_others]
              Xi_j <- dta[wh_others, cov_cols]
              prob_ind <- CalcNumerator(Ai_j = Ai_j, Xi_j = Xi_j,
                                        gamma_numer = gamma_numer,
                                        alpha = curr_alpha, re_alpha = re_alpha)
            }
            
            y_curr <- y_curr + dta$Y[ind] * prob_ind$prob
          }
        }
        
        denom <- Denominator(A = dta$V[neigh_ind[[nn]]],
                             X = dta[neigh_ind[[nn]], cov_cols],
                             phi_hat = phi_hat, alpha = curr_alpha,
                             integral_bound = integral_bound)
        denom <- length(neigh_ind[[nn]]) * denom$value * bern_prob
        
        yhat_group[nn, curr_it + 1, aa] <- pmax(y_curr / denom, 0.0000001)
      }
    }
  }
  if (keep_re_alpha) {
    return(list(yhat_group = yhat_group, re_alpha = re_alphas))
  }
  return(list(yhat_group = yhat_group))
}


###########################################
### Estimates and asymptotic variance of the population average potential
### (1)CalcScore
### (2)Ypop
###########################################

#' Matrix of score for the propensity score.
#' 
#' @param dta The data set including (at least) the treatment and covaratiates.
#' @param neigh_ind A list including one element for each neighborhood. That
#' element includes the indices of the observations in dta that belong to each
#' neighborhood. Can be left NULL if dta includes a column neigh.
#' @param phi_hat A list with two elements. The first one is a vector of
#' coefficients of the propensity score, and the second one is the random
#' effect variance.
#' @param cov_cols The indeces including the covariates of the propensity score
#' model.
#' @param trt_name The name of the treatment column. If it is 'A', you can
#' leave NULL.
#' 
#' @export
CalcScore <- function(dta, neigh_ind = NULL, phi_hat, cov_cols,
                      trt_name = NULL, integral_bound = 10) {
  
  dta <- as.data.frame(dta)
  if (is.null(neigh_ind)) {
    n_neigh <- max(dta$neigh)
    neigh_ind <- sapply(1 : n_neigh, function(x) which(dta$neigh == x))
  }
  n_neigh <- length(neigh_ind)
  
  re_var <- phi_hat[[2]]
  num_gamma <- length(phi_hat$coefs) + (re_var > 0)
  phi_hat$coefs <- matrix(phi_hat$coefs, ncol = 1)
  
  if (is.null(trt_name)) {
    trt_name <- 'V'
  }
  trt_col <- which(names(dta) == trt_name)
  
  scores <- matrix(NA, nrow = num_gamma, ncol = n_neigh)
  
  for (nn in 1 : n_neigh) {
    
    Ai <- dta[neigh_ind[[nn]], trt_col]
    Xi <- dta[neigh_ind[[nn]], cov_cols]
    
    hess_function <- function(gamma) {
      
      phi_hat_hess <- NULL
      if (re_var == 0) {
        phi_hat_hess$coefs <- gamma
        phi_hat_hess$re_var <- 0
      } else {
        phi_hat_hess$coefs <- gamma[- num_gamma]
        phi_hat_hess$re_var <- gamma[num_gamma]
      }
      
      likelihood <- Denominator(A = Ai, X = Xi, phi_hat = phi_hat_hess,
                                include_alpha = FALSE,
                                integral_bound = integral_bound)
      return(log(likelihood$value))
    }
    
    hess_x <- as.numeric(phi_hat$coefs)
    if (re_var > 0) {
      hess_x <- c(hess_x, re_var)
    }
    scores[, nn] <- numDeriv::grad(hess_function, x = hess_x)
    
  }
  return(scores)  
}





#' Geometric Direct effect estimates and asymptotic variance.
#' 
#' @param ypop A matrix with rows corresponding to the potential outcome under
#' control and treatment, and columns corresponding to the cluster-average
#' propensity of treatment.
#' @param ypop_var An 3-dimensional array, where the first two dimensions are
#' equal to 2 and include the variance covariance matrix of the population
#' average potential outcome for each alpha. Dimension 3 is alpha.
#' @param boots The results of BootVar() function including estimates of the
#' potential outcomes from the bootstrap samples.
#' @param alpha The values of alpha we consider. If ypop has column names,
#' alpha can be left null.
#' @param alpha_level Numeric. The alpha level of the confidence intervals
#' based on the quantiles of the bootstrap estimates.
#' 
#' @return A matrix with rows including the estimate and variance of the direct
#' effect and columns corresponding to alpha.
#' 
#' @export
#' 
GM_DE <-function(ygroup, boots = NULL, alpha = NULL,alpha_level = 0.05, scores = NULL,dta = NULL, use = 'everything'){
  quants <- c(0, 1) + c(1, - 1) * alpha_level / 2
  norm_quant <- - qnorm(alpha_level / 2)
  
  if (is.null(alpha)) {
    if (is.null(colnames(ypop))) {
      stop('Specify alpha.')
    }
    alpha <- as.numeric(colnames(ypop))
  }
  
  dim_names <- c('est', 'var', 'low_int', 'high_int')
  if (!is.null(boots)) {
    dim_names <- c(dim_names, 'boot_var', 'boot_var_LB', 'boot_var_UB',
                   'boot_low_quant', 'boot_high_quant')
  }
  
  de <- array(NA, dim = c(length(dim_names), length(alpha)))
  dimnames(de) <- list(stat = dim_names, alpha = alpha)
  
  DE <- exp(apply( (log(ygroup[, 2, ]) - log(ygroup[, 1, ])), 2, mean))
  var_log <- DEvar(ygroup,scores = scores,dta = dta, use = use)
  de[1, ] <- DE
  de[2, ] <- DE^2 * var_log$DE_var
  de[3, ] <- de[1, ] - norm_quant * sqrt(de[2, ])
  de[4, ] <- de[1, ] + norm_quant * sqrt(de[2, ])
  
  if (!is.null(boots)) {
    de[5, ] <- apply(boots[2, , ] - boots[1, , ], 1, var)
    de[6, ] <- de[1, ] - norm_quant * sqrt(de[5, ])
    de[7, ] <- de[1, ] + norm_quant * sqrt(de[5, ])
    de[8 : 9, ] <- apply(boots[2, , ] - boots[1, , ], 1, quantile,
                         probs = quants)
  }
  return(de)
}



DEvar <- function(ygroup,scores = NULL,dta = NULL, use = 'everything') {
# for the equations of cluster-level average
  n_neigh <- dim(ygroup)[1]
  alpha <- as.numeric(dimnames(ygroup)[[3]])
  
  logdiff <- log(ygroup[, 2, ]) - log(ygroup[, 1, ])
  ypop_var <- apply(logdiff, 2, var)
  # In order to get 1 / N, instead of 1 / (N - 1) in the variance estimates.
  ypop_var <- ypop_var * (n_neigh - 1) / n_neigh
  
  #For the equations of PS
  
  neigh_ind <- sapply(1 : n_neigh, function(x) which(dta$neigh == x))
  
  if (is.null(scores)) {
    stop('scores needs to be specified.')
  }
  
  num_gamma <- dim(scores)[1]
  
  var_est_ps <- rep(NA, length(alpha))
  names(var_est_ps) <- alpha
  
  # --- Calculating B11, the information matrix of the cluster ps.
  B11 <- matrix(0, nrow = num_gamma, ncol = num_gamma)
  for (nn in 1 : n_neigh) {
    scores_nn <- scores[, nn, drop = FALSE]
    B11 <- B11 + scores_nn %*% t(scores_nn)
  }
  B11 <- B11 / n_neigh
  B11_inv <- chol2inv(chol(B11))
  
  # ---- Calculating B21, F21 for a single treatment level.
  for (aa in 1:length(alpha)) {
    B21 <- numeric(num_gamma)
    F21 <- numeric(num_gamma)

    for (nn in 1:n_neigh) {
      logdiff_mean <- mean(logdiff[, aa])
      diff_vector <- logdiff[nn, aa] - logdiff_mean
      F21 <- F21 + scores[, nn] * diff_vector  # Ensure this multiplication is size-appropriate
    }
    
    # Normalize B21 and F21 by the number of neighborhoods
    B21 <- B21 / n_neigh
    F21 <- F21 / n_neigh
    
    # Calculate variance component for the adjusted variance estimation
    chol_B11_inv <- chol(B11_inv)
    mat1 <- (B21 - F21) %*% chol_B11_inv %*% B21
    mat <- B21 %*% B11_inv %*% F21
    
    # Store or utilize mat1 and mat as needed for further variance calculations
    var_est_ps[aa] <- mat1 - mat + ypop_var[aa]
    var_est_ps[aa] <- var_est_ps[aa] / n_neigh
  }
  return(list(DE_var = var_est_ps))
  
}


#' Geometric Indirect effect estimates and asymptotic variance.
#' 
#' @param ygroupM An matrix including the group average potential outcome
#' estimates where rows correspond to group, and columns to values of alpha.
#' #' @param boots The results of BootVar() function including estimates of the
#' potential outcomes from the bootstrap samples.
#' @param ps String. Can take values 'true', or 'estimated' for known or
#' estimated propensity score. Defaults to 'true'.
#' @param scores A matrix with rows corresponding to the parameters of the
#' propensity score model and columns for groups. Includes the score of the
#' propensity score evaluated for the variables of each group. Can be left
#' NULL for ps = 'true'.
#' @param alpha_level Numeric. The alpha level of the confidence intervals
#' based on the quantiles of the bootstrap estimates.
#' 
#' @export
GM_IE <- function(ygroupM, boots = NULL, ps = c('true', 'estimated'),
               scores = NULL, alpha_level = 0.05) {
  
  alpha <- as.numeric(dimnames(ygroupM)[[2]])
  quants <- c(0, 1) + c(1, - 1) * alpha_level / 2
  norm_quant <- - qnorm(alpha_level / 2)
  
  dim_names <- c('est', 'var', 'LB', 'UB')
  if (!is.null(boots)) {
    dim_names <- c(dim_names, 'boot_var', 'boot_var_LB', 'boot_var_UB',
                   'boot_low_quant', 'boot_high_quant')
  }
  
  ie <- array(NA, dim = c(length(dim_names), length(alpha), length(alpha)))
  dimnames(ie) <- list(stat = dim_names, alpha1 = alpha, alpha2 = alpha)
  
  for (a1 in 1 : length(alpha)) {
    for (a2 in 1 : length(alpha)) {
      ie[1, a1, a2] <- exp(mean(log(ygroupM[, a2]) - log(ygroupM[, a1]))) # estimate of SE
      ie[2, a1, a2] <- IEvar(ygroupV = ygroupM, a1, a2, scores = scores)   # variance of SE
      ie_sd <- sqrt(ie[2, a1, a2])
      ie[c(3, 4), a1, a2] <- ie[1, a1, a2] + norm_quant * c(- 1, 1) * ie_sd
    }
  }
  
  if (!is.null(boots)) {
    ie_var_boots <- array(NA, dim = c(length(alpha), length(alpha)))
    for (a1 in 1 : length(alpha)) {
      for (a2 in 1 : length(alpha)) {
        ie[5, a1, a2] <- var(boots[1, a1, ] - boots[1, a2, ])
        ie_sd <- sqrt(ie[5, a1, a2])
        ie[c(6, 7), a1, a2] <- ie[1, a1, a2] + norm_quant * c(- 1, 1) * ie_sd
        ie[c(8, 9), a1, a2] <- quantile(boots[1, a2, ] - boots[1, a1, ],
                                        probs = quants)
      }
    }
  }
  
  return(ie)
}


IEvar <- function(ygroupV, a1, a2, scores){
  n_neigh <- dim(ygroupV)[1]
  alpha <- as.numeric(dimnames(ygroupV)[[2]])
  
  ie_var <- cov(log(ygroupV))
  ie_var <- ie_var * (n_neigh - 1) / (n_neigh ^ 2)
  
  # Check scores input
  if (is.null(scores)) {
    stop('Provide score matrix.')
  }
  
  var_est_ps <- array(0, dim = dim(ie_var))
  num_gamma <- dim(scores)[1]
  
  # --- B11: Information matrix of the cluster ps
  B11 <- matrix(0, nrow = num_gamma, ncol = num_gamma)
  for (nn in 1:n_neigh) {
    scores_nn <- scores[, nn, drop = FALSE]
    B11 <- B11 + scores_nn %*% t(scores_nn)
  }
  B11 <- B11 / n_neigh
  B11_inv <- chol2inv(chol(B11))
  
  # Initialize C21 and D12 as vectors
  C21 <- rep(0, num_gamma)
  D12 <- rep(0, num_gamma)
  
  se_mean <- numeric(n_neigh)
  
  for (nn in 1:n_neigh) {
    log_y_nn <- log(ygroupV[nn, ])  # Extract row for clarity
    se_mean[nn] <- log_y_nn[a2] - log_y_nn[a1]
  }
  
  mean_se_mean <- mean(se_mean)
  
  for (nn in 1:n_neigh) {
    log_y_nn <- log(ygroupV[nn, ])
    logdiff <- log_y_nn[a2] - log_y_nn[a1]  # Scalar difference
    
    # Update D12 as a vector
    D12 <- D12 + scores[, nn, drop = FALSE] * (logdiff - mean_se_mean)
  }
  
  C21 <- C21 / n_neigh
  D12 <- D12 / n_neigh
  
  chol_B11_inv <- chol(B11_inv)
  mat1 <- (C21 %*% t(chol_B11_inv)) %*% t(C21 %*% t(chol_B11_inv))
  mat <- C21 %*% B11_inv %*% D12
  
  var_est_ps <- mat1 + mat + t(mat)
  var_est_ps <- ie_var[a1, a2] + var_est_ps / n_neigh
  
  return(var_est_ps) 
}




