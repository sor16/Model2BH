Densevalm22 <- function(th,RC){
  phi_b=th[3]
  sig_b2=th[2]
  zeta=th[1]
  lambda=th[4:9]
  
  f=lambda[1:5]-lambda[6]
  l=c(log(RC$w_tild+exp(th[1])))
  
  varr=c(exp(RC$B%*%lambda))
  Sig_eps=diag(c(varr,0))
  R_Beta=(1+sqrt(5)*RC$dist/exp(phi_b)+5*RC$dist^2/(3*exp(phi_b)^2))*exp(-sqrt(5)*RC$dist/exp(phi_b))+diag(RC$n)*RC$nugget#36.7 micro
  Sig_x=rbind(cbind(RC$Sig_ab,RC$m1),cbind(RC$m2,exp(sig_b2)*R_Beta))#23 micro
  
  X=Matrix(rbind(cbind(1,l,Matrix(diag(l))%*%RC$A),RC$Z)) #Matrix 2.76 milli/matrix 48.2 micro
  L=t(chol(as.matrix(X%*%Sig_x%*%t(X)+Sig_eps)))#Matrix 3.96 milli/matrix 66.3 micro
  w=solve(L,RC$y-X%*%RC$mu_x)#Matrix 3.13 milli/matrix 66.3 micro
  p=as.numeric(-0.5%*%t(w)%*%w-sum(log(diag(L)))-
    (RC$v+5-1)/2*log(RC$v*RC$s+f%*%RC$P%*%f)+
    sig_b2-exp(sig_b2)/RC$mu_sb+zeta-exp(zeta)/RC$mu_c-0.5/RC$tau_pb2*(phi_b-RC$mu_pb)^2) #Matrix 2.64 milli/matrix 2.47 milli
  
  W=solve(L,X%*%Sig_x)#Matrix 1.56 milli/matrix 2.47 milli
  x_u=RC$mu_x+t(chol(Sig_x))%*%rnorm(RC$n+2)
  sss=(X%*%x_u)-RC$y+rbind(sqrt(varr)*as.matrix(rnorm(RC$N)),0)#Matrix 3.48 milli/matrix 23.6085 micro
  x=as.matrix(x_u-t(W)%*%solve(L,sss))
  yp=X %*% x #M 17.46 micro/m 1.5 
  yp=yp[1:RC$N,]
  ypo=yp+as.matrix(rnorm(RC$n))*sqrt(varr)
  
  D=-2*sum(log(dlnorm(exp(RC$y[1:RC$N,]),yp,sqrt(varr))))
  
  return(list("p"=p,"x"=x,"yp"=yp,"ypo"=ypo,"D"=D,"varr"=varr))
  
}