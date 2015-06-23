Densp2 <- function(th,RC){
  th=as.matrix(th)
  phi_b=th[3,]
  sig_b2=th[2,]
  zeta=th[1,]
  lambda=th[4:9, ,drop=FALSE]
  
  
  f=lambda[1:5,,drop=FALSE]-lambda[6,]
  l=as.vector(log(RC$w_tild+exp(th[1,])))
  
  varr=as.vector(exp(RC$B%*%lambda))
  Sig_eps=diag(c(varr,0))
  R_Beta=(1+sqrt(5)*RC$dist/exp(phi_b)+5*RC$dist^2/(3*exp(phi_b)^2))*exp(-sqrt(5)*RC$dist/exp(phi_b))+diag(RC$n)*RC$nugget
  Sig_x=rbind(cbind(RC$Sig_ab,matrix(0,nrow=2,ncol=RC$n)),cbind(matrix(0,nrow=RC$n,ncol=2),exp(sig_b2)*R_Beta))
  
  X=Matrix(rbind(cbind(1,l,Matrix(diag(l))%*%RC$A),RC$Z))
  L=t(chol(as.matrix(X%*%Sig_x%*%t(X)+Sig_eps)))
  w=solve(L,RC$y-X%*%RC$mu_x)
  p=as.numeric(-0.5%*%t(w)%*%w-sum(log(diag(L)))-
                 (RC$v+5-1)/2*log(RC$v*RC$s+t(f)%*%RC$P%*%f)+
                 sig_b2-exp(sig_b2)/RC$mu_sb+zeta-exp(zeta)/RC$mu_c-0.5/RC$tau_pb2*(phi_b-RC$mu_pb)^2)
  
  return(p)
  
}