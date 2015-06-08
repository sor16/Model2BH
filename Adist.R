
Adist <- function(w){

library(Matrix)
  # %Zdist creates the 
  # 
  # %Input w: water level measurements
  # 
  # %Outputs
  # %Z:     Matrix \mathbold{Z} linking unique water level measurements (mathbold{w}') to actual
  # %       water level measurements (w) such that \mathbold{w}=mathbold{Zw}'
  # %dist:  Matrix of distances between unique water level measurements
  # %       dist_{ij}=|w_{i}'-w_{j}'|
  # %n:     Number of unique measurements 
  # %N:     Number of measurements
#w=RC$w
#w=t(as.matrix(w))  
N=length(w)
O=t(as.matrix(w[1]))
A=matrix()
A[1,1]=1
e=1
for(ee in 2:N){
  if( w[ee]==w[ee-1]){
    
#A[ee,e]=1

  }else{
    e=e+1
#A[ee,e]=1
O[e]=w[ee]
  }
  }
O=t(O)
A=diag(1,nrow=N,ncol=N)
A=Matrix(A,sparse=TRUE)
w2=O
n=e

W=O

for(ee in 2:e){
  W=cbind(W,O)
}
  dist=abs(W-t(W))
return(dist,A,n,N)
}
