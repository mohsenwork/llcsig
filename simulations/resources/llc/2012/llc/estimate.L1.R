estimate.L1<-function(eqs, lambda=0.01, prev=NA) {
  #L1 regularized solution of a linear equations system.
  #INPUT:
  # eqs - System of equations eqs$K%*%b =eqs$k.
  # lambda - Regularization penalizing the number of parents.
  # prev - 
  #OUTPUT:
  # b - solution vector or B directs effects matrix if 
  #     indexing structure eqs$P given


  gamma<-2

  #function to be minimized
  f<-function(b) {
    sum( (eqs$K%*%b-eqs$k)^2 ) + lambda*(1/gamma)*sum( log(cosh(gamma*b)) )
  }

  #gradient
  g<-function(b) {
    as.vector(2*t(eqs$K)%*%(eqs$K%*%b-eqs$k) + lambda*tanh(gamma*b)) #notice that the 2 disappears
  }

  #determining the starting value
  #this is a fairly good estimate better to use some regularization already here
  if ( !all(is.na(prev)) ) {
    b0<-prev[t(eqs$P)]
  } else {
    b0<-as.vector( mpinv(t(eqs$K)%*%eqs$K+lambda*diag(ncol(eqs$K)))%*%t(eqs$K)%*%eqs$k )
  }

  if ( is.na(f(b0)) || is.infinite(f(b0)) ) {
    b0<-rep(0,ncol(eqs$K))
  }

  #this is the optimization
  opt<-optim(b0,f,gr=g,method ='BFGS',control=list(maxit=10000,reltol=1e-15))
  b<-opt$par
  #print(c(sum(g(b)^2),max(abs(g(b)))))

  if ( all( !is.na( eqs$P ) ) ) {
    b.to.B(b, eqs$P, eqs$n )
  } else {
    b
  }
} 