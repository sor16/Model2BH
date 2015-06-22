#Model2BH
library(stats)
library(Matrix)
library(ggplot2)

#axel/begin/08.06.15


 Nit=50000;
 dataset=15;

RC=list()

RC$mu_a=3.20
RC$mu_b=2.29
RC$sig_a=sqrt(1.21)
RC$sig_b=sqrt(0.48)
RC$p_ab=-0.61
RC$mu_c=1.9000
RC$nugget=10^-8


RC$mu_sb=0.5
RC$mu_pb=0.5
RC$tau_pb2=0.25^2
RC$s=3
RC$v=5



# %import data from text file that has water level measurements in cm in left
# %column and corresponding discharge measurements in m^3/s in right column
# wq=importdata([num2str(dataset) '.txt']);
# %wq=importdata(['Jokdal.txt']);

 wq = as.matrix(read.table('15.txt'))

RC$y=rbind(as.matrix(log(wq[,2])),0)
RC$w=0.01*wq[,1]
RC$w_tild=RC$w-min(RC$w)
# 
 H=RC$w
 Q=wq[,2]
 dat=data.frame(H,Q)
 
 
 ggplot(dat,aes(x=H,y=Q))+geom_point(shape=1)+theme_bw()

#axel/end/08.06.15

Adist1 <- Adist(RC$w)
RC$A=Adist1$A
RC$dist=Adist1$dist
RC$n=Adist1$n
RC$N=Adist1$N
 
 
RC$A=Matrix(RC$A,sparse=TRUE)
RC$P=diag(nrow=5,ncol=5,6)-matrix(nrow=5,ncol=5,1)
 
RC$Sig_ab= rbind(c(RC$sig_a^2, RC$p_ab*RC$sig_a*RC$sig_b), c(RC$p_ab*RC$sig_a*RC$sig_b, RC$sig_b^2))

RC$mu_x=as.matrix(c(RC$mu_a,RC$mu_b, rep(0,RC$n))) #Setja i RC
 
RC$B=B_splines(t(RC$w_tild)/RC$w_tild[length(RC$w_tild)])
RC$Z=cbind(t(rep(0,2)),t(rep(1,RC$n)))

Dens = function(th) {-Densevalm22(th,RC)$p}

Densmin=optim(par=as.matrix(rep(0,9)),Dens,method="BFGS",hessian=TRUE)

t_m =Densmin$par
H=Densmin$hessian



phi_b=t_m[3]
sig_b2=t_m[2]
zeta=t_m[1]
lambda=t_m[4:9]

l=log(RC$w_tild+exp(t_m[1])) #as.matrix

varr_m=exp(RC$B%*%lambda)
Sig_eps=diag(as.numeric(rbind(varr_m,0)))
R_Beta=(1+sqrt(5)*RC$dist/exp(phi_b)+5*RC$dist^2/(3*exp(phi_b)^2))*exp(-sqrt(5)*RC$dist/exp(phi_b))+diag(1,RC$n,RC$n)*RC$nugget
Sig_x=rbind(cbind(RC$Sig_ab,matrix(0,nrow=2,ncol=RC$n)),cbind(matrix(0,nrow=RC$n,ncol=2),exp(sig_b2)*R_Beta))

X=Matrix(rbind(cbind(matrix(1,dim(l)),l,Matrix(diag(as.numeric(l)),sparse=TRUE)%*%RC$A),RC$Z),sparse=TRUE)

#check

L=t(chol(as.matrix(X%*%Sig_x%*%t(X)+Sig_eps)))

w=solve(L,(-RC$y+X%*%RC$mu_x))
mu=RC$mu_x-Sig_x%*%(t(X)%*%(solve(t(L),w)))

#
                      
ymu=X%*%mu
ymu=ymu[1:RC$N]
plot(RC$w,exp(ymu))

W=solve(L,(X%*%Sig_x))
vartem=diag(X%*%(Sig_x-t(W)%*%W)%*%t(X))

vartem=vartem[1:RC$N]
varaprr=+vartem+varr_m


#%emp bayes
ymu #%mat
oryggisbil=cbind(ymu+qnorm(0.025,0,sqrt(varaprr)), ymu+qnorm(0.975,0,sqrt(varaprr))) #oryggisbil a log

#[norminv(0.025,0,sqrt(varaprr)) norminv(0.975,0,sqrt(varaprr))]





#LH=chol(H)'/0.42
LH=t(chol(H))/0.8
#LH=chol(H)';


#accept=0;



t1=matrix(0,9,Nit)
t2=matrix(0,9,Nit)
t3=matrix(0,9,Nit)
t4=matrix(0,9,Nit)


xsiz=max(dim(mu));
x1=matrix(0,xsiz,Nit)
x2=matrix(0,xsiz,Nit)
x3=matrix(0,xsiz,Nit)
x4=matrix(0,xsiz,Nit)


for(j in 1:4){
  t_old=t_m
  t=matrix(0,nrow=9,ncol=Nit)
  x=matrix(0,nrow=xsiz,ncol=Nit)
  yp=matrix(0,nrow=RC$N,ncol=Nit)
  ypo=matrix(0,nrow=RC$N,ncol=Nit)
  varr=matrix(0,nrow=RC$N,ncol=Nit)
  D=matrix(0,nrow=1,ncol=Nit)
  
  
  
  Dens<-Densevalm22(t_old,RC)
  p_old=Dens$p
  x_old=Dens$x
  yp_old=Dens$yp
  ypo_old=Dens$ypo
  D_old=Dens$D
  varr_old=Dens$varr
  
  for(i in 1:Nit){
    t_new=t_old+solve(t(LH),as.matrix(rnorm(9,0,1)))
    
    Densnew<-Densevalm22(t_new,RC)
    p_new=Densnew$p
    x_new=Densnew$x
    yp_new=Densnew$yp
    ypo_new=Densnew$ypo
    D_new=Densnew$D
    varr_new=Densnew$varr
    
    logR=p_new-p_old
    
    if (logR>log(runif(1))){
      t_old=t_new
      x_old=x_new
      p_old=p_new
      yp_old=yp_new
      ypo_old=ypo_new
      D_old=D_new
      varr_old=varr_new
    }
    
    t[,i]=t_old
    yp[,i]=yp_old
    ypo[,i]=ypo_old
   
    D[1,i]=D_old
    varr[,i]=varr_old
  }
  
  if(j==1){
    t1=t
    yp1=yp
    ypo1=ypo
    D1=D
    varr1=varr
  } else if(j==2){
    t2=t
    yp2=yp
    ypo2=ypo
    D2=D
    varr2=varr
  } else if(j==3){
    t3=t
    yp3=yp
    ypo3=ypo
    D3=D
    varr3=varr
  } else if(j==4){
    t4=t
    yp4=yp
    ypo4=ypo
    D4=D
    varr4=varr
  }
}


Dhat=(-2)*sum(log(dlnorm(exp(RC$y[1:RC$N]),ymu,sqrt(varr_m))))
Davg=mean(c(D1[seq(2000,20000,5)],D2[seq(2000,20000,5)],D3[seq(2000,20000,5)],D4[seq(2000,20000,5)]))
pd=Davg-Dhat
DIC=Dhat+2*pd
B=1/(mean(0.5*c(D1[seq(2000,20000,5)],D2[seq(2000,20000,5)],D3[seq(2000,20000,5)],D4[seq(2000,20000,5)])))

c(Dhat, Davg, DIC, pd, B) #afhverju thessi vigur?



#              % 
#              % parfor j=1:4
#              %     t_old=t_m;
#              %     t=zeros(9,20000);
#              %     for i=1:5000
#              %     t_new=t_old+(LH'\normrnd(0,1,[9,1]));
# %     logR=Dens(t_old)-Dens(t_new);
# %         if logR>log(rand(1))
# %         t_old=t_new;        
# %         end
# %     t(:,i)=t_old;
# %     end
# % end
# % toc
# % 
# % for i=1:5000
# % d(i)=(10^-3+(x(4:8,i)-x(9,i))'*P*(x(4:8,i)-x(9,i)))/chi2rnd(10^-3+5);
#         % end
#         % 
#         % v=10^-3+5;
#         % s2=(10^-3+(t_m(4:8)-t_m(9))'*P*(t_m(4:8)-t_m(9)))/v;
# % 
# % v*s2/(v+2);
# % 
