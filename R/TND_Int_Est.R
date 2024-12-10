# Using control data
glmer_fit_col <- lme4::glmer(
  data = datTND[datTND$Y==0, ],
  formula = V ~ C + (1 | block),
  family = stats::binomial,
  nAGQ = 2
)


phi_hat <- list(coefs = summary(glmer_fit_col)$coef[, 1], re_var = re_var)

glmer_fit <- lme4::glmer(
  data = datTND,
  formula = V ~ C + (1 | block),
  family = stats::binomial,
  nAGQ = 2
)

phi_hat <- list(coefs = summary(glmer_fit)$coef[, 1], re_var = re_var)

resB <- GroupIPW(dta = datTND, cov_cols = 'C', phi_hat = phi_hat, gamma_numer = NULL, alpha = c(0.3, 0.5),
                 neigh_ind = NULL, trt_col = NULL, out_col = NULL, 
                 alpha_re_bound = 10, integral_bound = 10,
                 keep_re_alpha = FALSE, estimand = '1',
                 verbose = TRUE)
ygroup = resB$yhat_group

resGLMM <- GroupIPW(dta = datTND, cov_cols = 'C', phi_hat = phi_hat, gamma_numer = NULL, alpha = c(0.3, 0.5),
                    neigh_ind = NULL, trt_col = NULL, out_col = NULL, 
                    alpha_re_bound = 10, integral_bound = 10,
                    keep_re_alpha = FALSE, estimand = '2',
                    verbose = TRUE)
ygroup = resGLMM$yhat_group
head(resGLMM$yhat_group)

Score.est <- CalcScore(dta = datTND, neigh_ind = datTND$block, phi_hat = phi_hat, cov_cols = 'C',
                       trt_name = NULL, integral_bound = 10)
Ypop.est <- Ypop(ygroup, ps = 'estimated', scores = Score.est,
                 dta = NULL, use = 'everything')
DE(ypop = Ypop.est$ypop, ypop_var = Ypop.est$ypop_var, boots = NULL, alpha = NULL,
   alpha_level = 0.05)

#. indices corresponding to treatment = 0,  and [ygroup[, 2, ] is for treatment = 0]
head(ygroup[, 1 , ])
IE(ygroupV = ygroup[, 1, ], scores = Score.est)




PS_model <- function(TND, method = c('ORG_IPW', 'TND_IPW')){
  if (method == 'TND_IPW') {
    glmer_fit_col <- lme4::glmer(
      data = datTND[datTND$Y==0, ],
      formula = V ~ C + (1 | block),
      family = stats::binomial,
      nAGQ = 2
    )
    phi_hat <- list(coefs = summary(glmer_fit_col)$coef[, 1], re_var = re_var)
  }
  
  if (method == 'ORG_IPW') {
    glmer_fit <- lme4::glmer(
      data = datTND,
      formula = V ~ C + (1 | block),
      family = stats::binomial,
      nAGQ = 2
    )
    phi_hat <- list(coefs = summary(glmer_fit)$coef[, 1], re_var = re_var)
  }
  return(phi_hat)
}

r <- 50
# Initialize a list to store results
resDE.org <- vector("list", length = r)
resSE.org <- vector("list", length = r)

# Run the function r times and store the results
for (i in 1:r) {
  datTND <- datagen_int(nblocks=1000)
  phi_hat <- PS_model(datTND, method = 'ORG_IPW')
  resB <- GroupIPW(dta = datTND, cov_cols = 'C', phi_hat = phi_hat, gamma_numer = NULL, alpha = c(0.3, 0.5),
                   neigh_ind = NULL, trt_col = NULL, out_col = NULL, 
                   alpha_re_bound = 10, integral_bound = 10,
                   keep_re_alpha = FALSE, estimand = '1',
                   verbose = TRUE)
  ygroup = resB$yhat_group
  
  Score.est <- CalcScore(dta = datTND, neigh_ind = datTND$block, phi_hat = phi_hat, cov_cols = 'C',
                         trt_name = NULL, integral_bound = 10)
  Ypop.est <- Ypop(ygroup, ps = 'estimated', scores = Score.est,
                   dta = NULL, use = 'everything')
  resDE.org[[i]] <- DE(ypop = Ypop.est$ypop, ypop_var = Ypop.est$ypop_var, boots = NULL, alpha = NULL,
     alpha_level = 0.05)
  
  #. indices corresponding to treatment = 0,  and [ygroup[, 2, ] is for treatment = 0]
  #resSE[[i]] <- IE(ygroupV = ygroup[, 1, ], scores = Score.est)
}


# Initialize a list to store results
resDE.tnd <- vector("list", length = r)
resSE.tnd <- vector("list", length = r)

# Run the function r times and store the results
for (i in 1:r) {
  datTND<-datagen_int(nblocks=1000)
  phi_hat <- PS_model(datTND, method = 'TND_IPW')
  resB <- GroupIPW(dta = datTND, cov_cols = 'C', phi_hat = phi_hat, gamma_numer = NULL, alpha = c(0.3, 0.5),
                   neigh_ind = NULL, trt_col = NULL, out_col = NULL, 
                   alpha_re_bound = 10, integral_bound = 10,
                   keep_re_alpha = FALSE, estimand = '1',
                   verbose = TRUE)
  ygroup = resB$yhat_group
  
  Score.est <- CalcScore(dta = datTND, neigh_ind = datTND$block, phi_hat = phi_hat, cov_cols = 'C',
                         trt_name = NULL, integral_bound = 10)
  Ypop.est <- Ypop(ygroup, ps = 'estimated', scores = Score.est,
                   dta = NULL, use = 'everything')
  resDE.tnd[[i]] <- DE(ypop = Ypop.est$ypop, ypop_var = Ypop.est$ypop_var, boots = NULL, alpha = NULL,
                       alpha_level = 0.05)
  
  #. indices corresponding to treatment = 0,  and [ygroup[, 2, ] is for treatment = 0]
  resSE[[i]] <- IE(ygroupV = ygroup[, 1, ], scores = Score.est)

}

# Calculate the average among the first row of the output for the r replicates
apply(sapply(resDE.org, function(x) x[1, ]), 1, mean)
apply(sapply(resDE.tnd, function(x) x[1, ]), 1, mean)
# Calculate the average over the r replicates of the output
sapply(resSE, function(arr) apply(arr[1,,], 1, mean))




datTND<-datagen_int(nblocks=1000)
phi_hat <- PS_model(datTND, method = 'ORG_IPW')
resB <- GroupIPW(dta = datTND, cov_cols = 'C', phi_hat = phi_hat, gamma_numer = NULL, alpha = c(0.3, 0.5),
                 neigh_ind = NULL, trt_col = NULL, out_col = NULL, 
                 alpha_re_bound = 10, integral_bound = 10,
                 keep_re_alpha = FALSE, estimand = '1',
                 verbose = TRUE)
ygroup = resB$yhat_group

Score.est <- CalcScore(dta = datTND, neigh_ind = datTND$block, phi_hat = phi_hat, cov_cols = 'C',
                       trt_name = NULL, integral_bound = 10)
Ypop.est <- Ypop(ygroup, ps = 'estimated', scores = Score.est,
                 dta = NULL, use = 'everything')
de <- DE(ypop = Ypop.est$ypop, ypop_var = Ypop.est$ypop_var, boots = NULL, alpha = NULL,
   alpha_level = 0.05)

#. indices corresponding to treatment = 0,  and [ygroup[, 2, ] is for treatment = 0]
ie <- IE(ygroupV = ygroup[, 1, ], scores = Score.est)

alpha_range <- quantile(datTND$V, probs = c(0.2, 0.8))
alpha <- seq(alpha_range[1], alpha_range[2], length.out = num_alphas)
alpha <- sort(c(alpha, 0.1, 0.4))



plot_boot <- 3
index_low <- ifelse(plot_boot == 1, 3, ifelse(plot_boot == 2, 6, 8))
index_high <- index_low + 1

de_plot <- data.frame(alpha = alpha, de = de[1, ],
                      low = de[index_low, ], high = de[index_high, ])
a1 <- which(alpha == 0.1)
ie1_plot <- data.frame(ie = ie[1, a1, ], low = ie[index_low, a1, ],
                       high = ie[index_high, a1, ])

a1 <- which(alpha == 0.4)
ie2_plot <- data.frame(ie = ie[1, a1, ], low = ie[index_low, a1, ],
                       high = ie[index_high, a1, ])

res_array <- array(NA, dim = c(length(alpha), 3, 3))
dimnames(res_array) <- list(alpha = alpha, quant = c('DE', 'IE1', 'IE2'),
                            stat = c('est', 'LB', 'HB'))
res_array[, 1, ] <- as.matrix(de_plot[, - 1])
res_array[, 2, ] <- as.matrix(ie1_plot)
res_array[, 3, ] <- as.matrix(ie2_plot)

res_df <- plyr::adply(res_array[, , 1], 1 : 2)
res_df$LB <- c(de_plot$low, ie1_plot$low, ie2_plot$low)
res_df$UB <- c(de_plot$high, ie1_plot$high, ie2_plot$high)

f_names <- list('DE' = expression(DE(alpha)),
                'IE1' = expression(IE(0.1,alpha)),
                'IE2' = expression(IE(0.4,alpha)))
f_labeller <- function(variable, value){
  return(f_names[value])
}
res_df$alpha <- as.numeric(levels(res_df$alpha))[res_df$alpha]


ggplot(data = res_df, aes(x = alpha, y = V1, group = quant)) +  geom_line() +
  facet_wrap(~ quant, nrow = 1, labeller = f_labeller) +
  geom_ribbon(data = res_df, aes(ymin = LB, ymax = UB, group = quant), alpha = 0.3) +
  xlab(expression(alpha)) + ylab('') +
  theme(axis.title = element_text(size = 12),
        strip.text = element_text(size = 13),
        axis.text = element_text(size = 10)) +
  scale_x_continuous(breaks = seq(0.1, 0.4, by = 0.1)) +
  geom_hline(yintercept = 0, linetype = 2)


res <- list(yhat_group = yhat_group, scores = scores, yhat_pop = yhat_pop,
            yhat_pop_var = yhat_pop_var, boots = boots, de = de, ie = ie)
res$specs <- list(seed = 1234, B = B, num_alphas = num_alphas,
                  date = Sys.Date(), ps_with_re = ps_with_re,
                  estimand = estimand,
                  numerator_with_re = numerator_with_re)

# save(res, file = '~/Documents/Research/Interference/Revisions/Application/ps_wo_re/numerator_with_re/results.dat')

library(plot3D)
par(mar = c(1, 1, 1, 3))
persp3D(alpha, alpha, ie[1, , ], theta=50, phi=50, axes=TRUE,
        nticks=5, ticktype = "detailed", xlab='α1', ylab='α2', zlab='',
        colkey = list(length = 0.5, shift = -0.1))
