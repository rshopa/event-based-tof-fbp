############################
# author: Roman Shopa
# Roman.Shopa[at]ncbj.gov.pl
############################
# These functions have general math usage for multiple issues

# normalise intensity
normalize <- function(X, Type="ranges"){
  if(Type=="ranges") return((X-min(X))/(max(X)-min(X))) 
  else if(Type=="max") return(X/max(abs(X)))
  else return(X/sum(X))
}
# pixel ID and vectorised version for voxel - for negatives and positives
findPixelID  <- function(X, Delta) as.integer(round(X/Delta))
# non-vectorised version using the scale, i.e. vector of coordinates
findPixelIDByScale <- function(X, Scale.1D) which.min(abs(Scale.1D - X))
# ID of R vector element - starting from 1, for positive only!
findVectorID  <- function(X, Delta) as.integer(floor(X/Delta))+1L
# row of an unrolled 2D array
detectRow <- function(i,j,imax) (j-1L)*imax + i
# row of an unrolled 3D array
detect3DRow <- function(i,j,k,imax,jmax) (k-1L)*imax*jmax + (j-1L)*imax + i
# check if value is inside span []
inInterval <- function(x, interval) x <= interval[2] & x >= interval[1]
# distance between points
distanceBetweenPoints <- function(P1,P2) sqrt(sum((P1-P2)^2))
# quadratic equation
solveQuadraticEq <- function(a,b,c){
  D <- b^2-4*a*c
  if(D<0) return(NA) else if(D==0) return(-b/(2*a)) else return((-b+c(-1,1)*sqrt(D))/(2*a))
}
