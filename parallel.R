# ypo=list()
ptm <- proc.time()
cl <- makeCluster(getOption('cl.cores', detectCores()))
clusterEvalQ(cl, { source("Adist.R");source("spline_functions.R");source("B_splines.R");source("Densp2fast.R");source("Densevalm22fast.R");source("model2.R");library(MCMCpack)})
m=matrix(NA,nrow=4,ncol=1)
# for(i in 1:4){
  post.samp <- rbind(parRapply(cl=cl,m,function(x) {MCMCmetrop1R(Densp2fast,theta.init=t_m,RC=RC,mcmc=50000)}))
  ypo=parRapply(cl,post.samp,function(x) Densevalm22fast(x,RC)$ypo) 
  #ypo=apply(post.samp,1,function(x) Densevalm22fast(x,RC)$ypo) 
  
proc.time()-ptm
# }
library(foreach)
post.temp=list()
post.samp <- foreach(i=1:4) %dopar% MCMCmetrop1R(Densp2fast,theta.init=t_m,RC=RC,mcmc=50000)