library(Matrix)
Densevalm22 <- function(th,RC){
  
  phi_b=th[3,]
  sig_b=th[2,]
  zeta=th[1,]
  lambda=th[4:9,]
  
  
  f=(lambda[1:5]-lambda[6])
  l=log(RC$w_tild)+exp(t[1,])
  
  varr=exp(RC$B*lambda)
  Sig_eps=diag(rbind(varr,0))
  R_Beta=(1+sqrt(5)*RC$dist/exp(phi_b)+5*RC$dist^2/(3*exp(phi_b)^2))*exp(-sqrt(5)*RC$dist/exp(phi_b))+diag(RC$n)*R$nugget
  Sig_x=rbind(cbind(RC$Sig_ab,matrix(0,nrow=2,ncol=RC$n)),cbind(matrix(0,nrow=25,ncol=RC$n),exp(sig_b2)%*%R_Beta))
  
  X=Matrix(rbind(rep(1,length(l)),l,rbind(Matrix(diag(l))*RC$A,RC$Z)[kornecker(1:length(RC$Z),c(0,length(RC$Z)))]))
  L=chol(as.matrix(X)%*%Sig_x%*%t(as.matrix(X))+Sig_eps)
  w=solve(L,RC$y-X%*%RC$mu_x)
  p=p=-0.5%*%t(w)%*%w-sum(log(diag(L)))-
    (RC$v+5-1)/2*log(RC$v*RC$s+t(f)%*%RC$P%*%f)+
    sig_b2-exp(sig_b2)/RC$mu_sb+zeta-exp(zeta)/RC.mu_c-0.5/RC.tau_pb2*(phi_b-RC$mu_pb)^2
  
  W=solve(L,X%*%Sig_x)
  x_u=RC$mu_x+t(chol(Sig_x))%*%as.matrix(rnorm(RC$n+2)
  sss=(X%*%x_u)-RC$y+rbind(sqrt(varr)*as.matrix(rnorm(RC$N)),0)
  x=x_u-t(W)%*%solve(L,sss)
  yp=X %*% x
  yp=yp[1:RC$N,]
  ypo=yp+as.matrix(rnorm(RC$n))*sqrt(varr)
  
  D=-2*sum(log(dlnorm(exp(RC$y[1:RC$N,]),yp,sqrt(varr))))
  
  return(list("pmin"=-p,"p"=p,"x"=x,"yp"=yp,"ypo"=ypo,"D"=D,"varr"=varr))
  
  
}