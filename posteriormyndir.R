##figure 3 rating curve logscale with conf interval

data=data.frame(W=RC$w, Q=RC$y)
data$l_m=l_m
seq=seq(2000,20000,5)
quantypo2=ypo2[,seq]
quantypo3=ypo3[,seq]
quantypo4=ypo4[,seq]
quantmatrix=t(cbind(quantypo1,quantypo2,quantypo3,quantypo4))
prctile=t(apply(quantmatrix, 2, quantile, probs = c(0.025,0.5, 0.975),  na.rm = TRUE))
data$fit=prctilep[,2]
data$lower=prctile[,1]
data$upper=prctile[,3]
rclog=ggplot(data)+geom_line(mapping=aes(l_m,fit))+theme_bw()+geom_point(mapping=aes(l_m,Q))+geom_line(mapping=aes(l_m,lower),linetype="dashed")+
  geom_line(mapping=aes(l_m,upper),linetype="dashed")

##figure 6 rating curve real scale with conf interval
rcreal=ggplot(data)+theme_bw()+geom_point(aes(exp(Q),W))+geom_line(aes(exp(fit),W))+geom_line(aes(exp(lower),W),linetype="dashed")+
  geom_line(aes(exp(upper),W),linetype="dashed")