(datTND<-datagen_int(nblocks=1000, return_full=T))

  
prevafun <- function(nn, OR_1=3,OR_2=0.8,OR_C=3,OR_WC=5,em=-0.5) {
  stats = matrix(NA, nrow = nn, ncol = 18)
  for (j in 1:nn) {
    da = datagen_int(return_full=T, em= em, OR_1=OR_1, OR_2 = OR_2, OR_C = OR_C, OR_WC = OR_WC)
    tnd = sample(which(da$H==1))
    stats[j, ] = c(mean(da$V), mean(da$Infec), 
                   mean(da$Infec_COVID), mean(da$Infec==1 & da$Infec_COVID==1), 
                   mean(da$W==1), 
                   mean(da$W[da$Infec_COVID==1]), 
                   mean(da$W[da$Infec==1|da$Infec_COVID==1]),
                   mean(da$W[da$Infec_COVID==1& da$V==1]), mean(da$W[da$Infec_COVID==1& da$V==0]), 
                   mean(da$H), mean(da$H[da$Infec_COVID==1]), mean(da$H[da$W==1]), 
                   mean(da$H[da$Infec_COVID==1& da$V==1]), mean(da$H[da$Infec_COVID==1& da$V==0]), 
                   
                   mean(da[tnd,]$Infec_COVID), mean(da[tnd,]$Infec),
                   sum(rowSums(da[tnd, c("Infec","Infec_COVID")])==2)/length(tnd), mean(da[tnd,]$V))
  }
  return(stats= stats)
}

prevalences_em1 = prevafun(nn=3, OR_1=3,OR_2=5)
prevalences_em2 = prevafun(nn=3, OR_1=3,OR_2=8)
prevalences_em3 = prevafun(nn=3, OR_1=3,OR_2=10)


prev_c = round(data.frame(rbind(colMeans(prevalences_em1), colMeans(prevalences_em2), colMeans(prevalences_em3))), digits =4)
colnames(prev_c)= c("p_V","p_I1","p_I2", "p_[I1=I2=1]", "p_W", "p_W[I2=1]", "p_W[I1=1|I2=1]",  
                    "p_W[V=1,I2=1]","p_W[V=0,I2=1]",
                    "p_H","p_H[I2=1]","p_H[W=1]","p_H[V=1,I2=1]","p_H[V=0,I2=1]",
                    "p_TND_I2", "p_TND_I1", "p_TND_[I1=I2=1]", "p_TND_V")
rownames(prev_c)= c("co_inf_para = 0.00001","co_inf_para = 0.001", "co_inf_para = 0.1")
t(prev_c)

