estimate.L2<-function(eqs, lambda=0.01) {
  #L2 regularized solution of a linear equations system.
  #INPUT:
  # eqs - System of equations eqs$K%*%b =eqs$k.
  # lambda - Regularization penalizing the number of parents.
  #OUTPUT:
  # b - solution vector or B directs effects matrix if 
  #     indexing structure eqs$P given

  b<-as.vector( mpinv(t(eqs$K)%*%eqs$K+lambda*diag(ncol(eqs$K)))%*%t(eqs$K)%*%eqs$k )

  if ( all( !is.na( eqs$P ) ) ) {
    b.to.B(b, eqs$P, eqs$n )
  } else {
    b
  }
} 

