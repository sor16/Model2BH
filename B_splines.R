B_splines <- function(ZZ){
  
#A script to test the B-splines.
  
#The number of equally spaced interior knots.
kx=2

#Delta x and Delta y.
dx=1/(kx+1)

#The order of the splines.
M = 4

#Determine the number of functions.
Nx = kx + M
# 
#The epsilon-knots
# epsilon_x = dx*[0:(kx+1)];
epsilon_x = dx*seq(0,kx+1,by=1)


#the tau-knots.
# tau_x = zeros(1,kx+2*M);
# tau_x(1:M) = epsilon_x(1)*ones(1,M);
# tau_x(M+1:kx+M) = epsilon_x(2:kx+1);
# tau_x(kx+M+1:kx+2*M) = epsilon_x(kx+2)*ones(1,M);
tau_x = matrix(0,nrow=1,ncol=(kx+2*M))
tau_x[1:M] = epsilon_x[1]*matrix(1,nrow=1,ncol=M)
tau_x[(M+1):(kx+M)]=epsilon_x[2:(kx+1)]
tau_x[(kx+M+1):(kx+2*M)]=epsilon_x[kx+2]*matrix(1,nrow=1,ncol=M)

#Vector with values of x and y.
lx = length(ZZ)
 
#Compute the x-splines and the y-splines.
#[XX] = spline_functions(ZZ,tau_x,dx,kx,M);
XX <- spline_functions(ZZ,tau_x,dx,kx,M)

  spline_functions <- function(ZZZ,tau,dx,k,M){
# function [XX] = spline_functions(ZZZ,tau,dx,k,M)
# 
# XX = zeros(k+M,length(ZZZ));
  XX = matrix(0,nrow=(k+M),ncol=length(ZZZ))
# 
# % i = 1
XX[1,] = (1/dx^3)*(tau[M+1]-ZZZ)*(tau[M+1]-ZZZ)*(tau[M+1]-ZZZ)*(tau[M]<=ZZZ)*(ZZZ<tau[M+1]);

# % i = 2
 XX[2,] = (1/dx^3)*(ZZZ-tau[2])*(tau[M+1]-ZZZ)*(tau[M+1]-ZZZ)*(tau[M]<=ZZZ)*(ZZZ<tau[M+1])+
  (1/2/dx^3)*(tau[M+2]-ZZZ)*(ZZZ-tau[3])*(tau[M+1]-ZZZ)*(tau[M]<=ZZZ)*(ZZZ<tau[M+1])+
  (1/4/dx^3)*(tau[M+2]-ZZZ)*(tau[M+2]-ZZZ)*(ZZZ-tau[M])*(tau[M]<=ZZZ)*(ZZZ<tau[M+1])+
  (1/4/dx^3)*(tau[M+2]-ZZZ)*(tau[M+2]-ZZZ)*(tau[M+2]-ZZZ)*(tau[M+1]<=ZZZ)*(ZZZ<tau[M+2])
 
# % i = 3
 XX[3,] = (1/2/dx^3)*(ZZZ-tau[3])*(ZZZ-tau[3])*(tau[M+1]-ZZZ)*(tau[M]<=ZZZ)*(ZZZ<tau[M+1])+
  (1/4/dx^3)*(ZZZ-tau[3])*(tau[M+2]-ZZZ)*(ZZZ-tau[M])*(tau[M]<=ZZZ)*(ZZZ<tau[M+1])+ 
  (1/4/dx^3)*(ZZZ-tau[3])*(tau[M+2]-ZZZ)*(tau[M+2]-ZZZ)*(tau[M+1]<=ZZZ)*(ZZZ<tau[M+2])+ 
  (1/6/dx^3)*(tau[M+3]-ZZZ)*(ZZZ-tau[M])*(ZZZ-tau[M])*(tau[M]<=ZZZ)*(ZZZ<tau[M+1])+
  (1/6/dx^3)*(tau[M+3]-ZZZ)*(ZZZ-tau[M])*(tau[M+2]-ZZZ)*(tau[M+1]<=ZZZ)*(ZZZ<tau[M+2])+ 
  (1/6/dx^3)*(tau[M+3]-ZZZ)*(tau[M+3]-ZZZ)*(ZZZ-tau[M+1])*(tau[M+1]<=ZZZ)*(ZZZ<tau[M+2])+ 
  (1/6/dx^3)*(tau[M+3]-ZZZ)*(tau[M+3]-ZZZ)*(tau[M+3]-ZZZ)*(tau[M+2]<=ZZZ)*(ZZZ<tau[M+3])
 
# % i = 4,...,k + 1
 for (kk in M:(k + 1)){
 XX[kk,] = (1/6/dx^3)*(ZZZ-tau[kk])*(ZZZ-tau[kk])*(ZZZ-tau[kk])*(tau[kk]<=ZZZ)*(ZZZ<tau[kk+1])+
  (1/6/dx^3)*(ZZZ-tau[kk])*(ZZZ-tau[kk])*(tau[kk+2]-ZZZ)*(tau[kk+1]<=ZZZ)*(ZZZ<tau[kk+2])+
  (1/6/dx^3)*(ZZZ-tau[kk])*(tau[kk+3]-ZZZ)*(ZZZ-tau[kk+1])*(tau[kk+1]<=ZZZ)*(ZZZ<tau[kk+2])+ 
  (1/6/dx^3)*(ZZZ-tau[kk])*(tau[kk+3]-ZZZ)*(tau[kk+3]-ZZZ)*(tau[kk+2]<=ZZZ)*(ZZZ<tau[kk+3])+
  (1/6/dx^3)*(tau[kk+4]-ZZZ)*(ZZZ-tau[kk+1])*(ZZZ-tau[kk+1])*(tau[kk+1]<=ZZZ)*(ZZZ<tau[kk+2])+ 
  (1/6/dx^3)*(tau[kk+4]-ZZZ)*(ZZZ-tau[kk+1])*(tau[kk+3]-ZZZ)*(tau[kk+2]<=ZZZ)*(ZZZ<tau[kk+3])+
  (1/6/dx^3)*(tau[kk+4]-ZZZ)*(tau[kk+4]-ZZZ)*(ZZZ-tau[kk+2])*(tau[kk+2]<=ZZZ)*(ZZZ<tau[kk+3])+
  (1/6/dx^3)*(tau[kk+4]-ZZZ)*(tau[kk+4]-ZZZ)*(tau[kk+4]-ZZZ)*(tau[kk+3]<=ZZZ)*(ZZZ<tau[kk+4]) 
 }
 
# % i = k + 2
XX[k+2,] =  -(1/6/dx^3)*(tau[k+2]-ZZZ)*(tau[k+2]-ZZZ)*(tau[k+2]-ZZZ)*(tau[k+2]<=ZZZ)*(ZZZ<tau[k+3])-
 (1/6/dx^3)*(tau[k+2]-ZZZ)*(tau[k+2]-ZZZ)*(ZZZ-tau[k+4])*(tau[k+3]<=ZZZ)*(ZZZ<tau[k+4])-
 (1/6/dx^3)*(tau[k+2]-ZZZ)*(ZZZ-tau[k+5])*(tau[k+3]-ZZZ)*(tau[k+3]<=ZZZ)*(ZZZ<tau[k+4])-
 (1/6/dx^3)*(tau[k+2]-ZZZ)*(ZZZ-tau[k+5])*(ZZZ-tau[k+5])*(tau[k+4]<=ZZZ)*(ZZZ<tau[k+5])-       
 (1/4/dx^3)*(ZZZ-tau[k+6])*(tau[k+3]-ZZZ)*(tau[k+3]-ZZZ)*(tau[k+3]<=ZZZ)*(ZZZ<tau[k+4])-
 (1/4/dx^3)*(ZZZ-tau[k+6])*(tau[k+3]-ZZZ)*(ZZZ-tau[k+5])*(tau[k+4]<=ZZZ)*(ZZZ<tau[k+5])-
 (1/2/dx^3)*(ZZZ-tau[k+6])*(ZZZ-tau[k+6])*(tau[k+4]-ZZZ)*(tau[k+4]<=ZZZ)*(ZZZ<tau[k+5])

# % i = k + 3
XX[k+3,] = - (1/4/dx^3)*(tau[k+3]-ZZZ)*(tau[k+3]-ZZZ)*(tau[k+3]-ZZZ)*(tau[k+3]<=ZZZ)*(ZZZ<tau[k+4])-
 (1/4/dx^3)*(tau[k+3]-ZZZ)*(tau[k+3]-ZZZ)*(ZZZ-tau[k+5])*(tau[k+4]<=ZZZ)*(ZZZ<tau[k+5])-
 (1/2/dx^3)*(tau[k+3]-ZZZ)*(ZZZ-tau[k+6])*(tau[k+4]-ZZZ)*(tau[k+4]<=ZZZ)*(ZZZ<tau[k+5])-
 (1/dx^3)*(ZZZ-tau[k+7])*(tau[k+4]-ZZZ)*(tau[k+4]-ZZZ)*(tau[k+4]<=ZZZ)*(ZZZ<tau[k+5])

# % i = k + 4
XX[k+4,] = -(1/dx^3)*(tau[k+4]-ZZZ)*(tau[k+4]-ZZZ)*(tau[k+4]-ZZZ)*(tau[k+4]<=ZZZ)*(ZZZ<=tau[k+5])

XX = t(XX)

  }
}
