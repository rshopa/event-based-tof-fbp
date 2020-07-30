############################
# author: Roman Shopa
# Roman.Shopa[at]ncbj.gov.pl
############################

# since no closed forms of high-pass filtering functions can be evaluated in image space,
# dummy profiles (Ram-Lak, inverse Gaussian) are created and then interpolated

# Arbitrary FBP filter (as for STIR)
# WARNING! Not shifted!
filterFBPfourier <- function(W.s, W.c=.5, Alpha=1.)
  if(W.c<=0) W.s*0 else abs(W.s)*ifelse(abs(W.s) > W.c, 0, Alpha+(1-Alpha)*cos(pi*W.s/W.c))
# no TOF regularisation
filterFBPfourierTOF <- function(W.s, Tau, W.c=.5, Alpha=1.){
  x <- (pi*Tau*W.s)^2
  ramp <- pi*exp(x)/integrate(function(fi) exp(-x*cos(fi)), 0, pi)[["value"]]
  if(abs(W.s)>W.c) 0 else ramp*(Alpha + (1-Alpha)*cos(pi*W.s/W.c))
}
# create FBP dummy for 3D kernel 
# IMPORTANT! Here Nyquist frequency is assumed to be w_Nq = 0.5 (in cycles), 
# hence W.s is in range (0,0.5]!
createFBPdummy <- function(Half.span = 25L, W.cut = .5, Alpha = 1., TOF=TRUE, Tau.reg = 10){
  if(TOF) fourier.FBP <- sapply(seq(-.5,.5,length.out = 2*Half.span+1),
                                filterFBPfourierTOF, Tau=Tau.reg, W.c=W.cut, Alpha=Alpha)
  else fourier.FBP <- filterFBPfourier(seq(-.5,.5,length.out = 2*Half.span+1), W.cut, Alpha)
  # shift to proper Fourier form (Nyquist at the centre)
  fourier.FBP <- c(tail(fourier.FBP,Half.span+1), head(fourier.FBP,Half.span))
  w.function <- Re(fft(fourier.FBP,inverse = T))
  # shift so that zero is at the centre + normalise!
  return(c(tail(w.function, Half.span), head(w.function, Half.span+1))/max(w.function))
}
# create dummy for inverse Gaussian
createInverseGaussDummy <- function(Sigma, Half.N, Cut.factor=NULL, Alpha=1L){ # Half.N - number N for points [-N:N]
  inv.DFT <- 1/abs(fft(dnorm(Sigma*(-Half.N:Half.N), 0, Sigma)))
  w.s <- .5/Half.N * c(0:Half.N, Half.N:1) # frequency
  if(!is.null(Cut.factor)) inv.DFT[inv.DFT>=(max(inv.DFT)/Cut.factor)]<-0
  w.c <- max(w.s[which.max(inv.DFT)])
  # inverse Fourier transform IFT, noise as in Ram-Lak (STIR)
  IFT <- Re(fft(inv.DFT*(Alpha + (1-Alpha)*cos(pi*w.s/w.c)), inverse = T)) # inverse 'anti-Gaussian', real part
  IFT <- c(tail(IFT,Half.N), head(IFT,Half.N+1L)) # shift so that zero is at the centre
  return(IFT/max(IFT)) # return normalized
}