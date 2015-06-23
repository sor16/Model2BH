#t_new=t_old+solve(t(LH),as.matrix(rnorm(2,0,1)))
library(MCMCpack)
# Dens <- function(th){Densevalm22(th,RC)}
Densmin=optim(par=c(0,0),Densp2,RC=RC,hessian=TRUE)
t_m=as.matrix(Densmin$par)
#post.samp <- MCMCmetrop1R(Densp,theta.init=t_m,RC=RC,mcmc=20000)
ptm <- proc.time()
for(i in 1:4){
  post.samp <- MCMCmetrop1R(Dens,theta.init=t_m,RC=RC,mcmc=20000)
  ypo=matrix(0,nrow=nrow(wq),ncol=Nit)
  for(i in 1:nrow(post.samp)){
    ypo[,i]=Densevalm11(post.samp[i,],RC)$ypo
    
  }
  if(j==1){
    ypo1=ypo
  } else if(j==2){
    ypo2=ypo
  } else if(j==3){
    ypo3=ypo
  } else if(j==4){
    ypo4=ypo
  }
}

proc.time() - ptm
#post.samp <- MCMCmetrop1R(Densp,theta.init=t_m,RC=RC,mcmc=20000)

ptm <- proc.time()
ypo=list()
for(i in 1:4){
post.samp <- MCMCmetrop1R(Densp2fast,theta.init=t_m,RC=RC,mcmc=20000)
ypo=apply(post.samp,1,function(x) Densevalm22fast(x,RC)$ypo) 
}
proc.time() - ptm