#Model2BH
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

Adist <- Adist(RC$w)
RC$A=Adist$A
RC$dist=Adist$dist
RC$n=Adist$n
RC$N=Adist$N
 
 
RC$A=Matrix(RC$A,sparse=TRUE)
RC$P=diag(nrow=5,nrow=5,6)-matrix(nrow=5,ncol=5,1)
 
 RC$Sig_ab= rbind((cbind(t(as.matrix(RC$sig_a^2)), t(as.matrix(RC$p_ab%*%RC$sig_a%*%RC$sig_b))) ), cbind)
# 
# 
# RC$Sig_ab=[RC.sig_a^2 RC.p_ab*RC.sig_a*RC.sig_b;RC.p_ab*RC.sig_a*RC.sig_b RC.sig_b^2]; %Setja Ã� RC
# RC$mu_x=[RC.mu_a RC.mu_b zeros(1,RC.n)]'; %Setja Ã� RC
# 
# RC.B=B_splines(RC.w_tild'/RC.w_tild(length(RC.w_tild)));
# RC.Z=[zeros(2,1);ones(RC.n,1)]';
# 
# Dens =@(t)-DensEvalm22(t,RC);
# 
# [t_m,~,~,~,~,H]=fminunc(Dens,zeros(9,1));
# 
# phi_b=t_m(3);
# sig_b2=t_m(2);
# zeta=t_m(1);
# lambda=t_m(4:9);
# 
# l=log(RC.w_tild+exp(t_m(1)));
# 
# varr_m=exp(RC.B*lambda);
# Sig_eps=diag([varr_m;0]);
# R_Beta=(1+sqrt(5)*RC.dist/exp(phi_b)+5*RC.dist.^2/(3*exp(phi_b)^2)).*exp(-sqrt(5)*RC.dist/exp(phi_b))+eye(RC.n)*RC.nugget;
# Sig_x=[RC.Sig_ab,zeros(2,RC.n);zeros(RC.n,2),exp(sig_b2)*R_Beta];
# 
# X=sparse([ones(size(l)),l,sparse(diag(l))*RC.A;RC.Z]);
# L=chol(full(X*Sig_x*X'+Sig_eps))';
# 
# w=L\(-RC.y+X*RC.mu_x);
# mu=RC.mu_x-Sig_x*(X'*(L'\w));
#                       
#                       ymu=X*mu;
#                       ymu=ymu(1:RC.N);
#                       hold on;plot(RC.w,exp(ymu))
#                       
#                       W=L\(X*Sig_x);
#                       vartem=diag(X*(Sig_x-W'*W)*X');
# vartem=vartem(1:RC.N)
# varaprr=+vartem+varr_m;
# 
# 
# %emp bayes
# ymu %mat
# [ymu+norminv(0.025,0,sqrt(varaprr)) ymu+norminv(0.975,0,sqrt(varaprr))] %Ã¶ryggisbil Ã¡ log
# 
# %[norminv(0.025,0,sqrt(varaprr)) norminv(0.975,0,sqrt(varaprr))]
# 
# 
# 
# 
# 
# %LH=chol(H)'/0.42;
# LH=chol(H)'/0.8;
# %LH=chol(H)';
# 
# 
# accept=0;
# 
# 
# 
# t1=zeros(9,Nit);
# t2=zeros(9,Nit);
# t3=zeros(9,Nit);
# t4=zeros(9,Nit);
# 
# 
# xsiz=max(size(mu));
# x1=zeros(xsiz,Nit);
# x2=zeros(xsiz,Nit);
# x3=zeros(xsiz,Nit);
# x4=zeros(xsiz,Nit);
# 
# tic
# for j=1:4
# t_old=t_m;
# t=zeros(9,Nit);
# x=zeros(xsiz,Nit);
# yp=zeros(RC.N,Nit);
# ypo=zeros(RC.N,Nit);
# varr=zeros(RC.N,Nit);
# D=zeros(1,Nit);
# [p_old,x_old,yp_old,ypo_old,D_old,varr_old]=DensEvalm22(t_old,RC);
# for i=1:Nit
# t_new=t_old+(LH'\normrnd(0,1,[9,1]));
#              [p_new,x_new,yp_new,ypo_new,D_new,varr_new]=DensEvalm22(t_new,RC);
#              logR=p_new-p_old;
#              if logR>log(rand(1))
#              t_old=t_new;
#              x_old=x_new;
#              p_old=p_new;
#              yp_old=yp_new;
#              ypo_old=ypo_new;
#              D_old=D_new;
#              varr_old=varr_new;
#              end
#              t(:,i)=t_old;
#              yp(:,i)=yp_old;
#              ypo(:,i)=ypo_old;
#              D(1,i)=D_old;
#              x(:,i)=x_old;
#              varr(:,i)=varr_old;
#              end
#              
#              
#              switch j
#              case 1
#              t1=t;
#              yp1=yp;
#              ypo1=ypo;
#              D1=D;
#              x1=x;
#              varr1=varr;
#              case 2
#              t2=t;
#              yp2=yp;
#              ypo2=ypo;
#              D2=D;
#              x2=x;
#              varr2=varr;
#              case 3
#              t3=t;
#              yp3=yp;
#              ypo3=ypo;
#              D3=D;
#              x3=x;
#              varr3=varr;
#              case 4
#              t4=t;
#              ypo4=ypo;
#              yp4=yp;
#              D4=D;
#              x4=x;
#              varr4=varr;
#              end
#              end
#              toc
#              
#              Dhat=-2*sum(log(lognpdf(exp(RC.y(1:RC.N)),ymu,sqrt(varr_m))));
#              Davg=mean([D1(2000:5:20000) D2(2000:5:20000) D3(2000:5:20000) D4(2000:5:20000)]);
#              pd=Davg-Dhat;
#              DIC=Dhat+2*pd;
#              B=1/mean(exp(0.5*[D1(2000:5:20000) D2(2000:5:20000) D3(2000:5:20000) D4(2000:5:20000)]));
#              [Dhat Davg DIC pd B]
#              
#              
#              
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
