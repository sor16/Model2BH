
library(Matrix)
Adist <- function(w){
  
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
  
N=length(w)
A=1
O=w[1]
e=1
for(ee in 2:N){
  if( w[ee]==w[ee-1]){
    
A[ee,e]=1
  }else{
    e=e+1
A[ee,e]=1
O[e]=w[ee]
  }
  }
O=t(O)

A=Matrix(A,sparse=TRUE)
w2=O
n=e

for(ee in 1:e){
  W=cbind(t(as.matrix(W)),O)
}
  dist=abs(W-t(W))
return(dist,A,n,N)
}


# W=[];
# for ee=1:e
#     W=[W O];
# end
# dist=abs(W-W');
# 
# end