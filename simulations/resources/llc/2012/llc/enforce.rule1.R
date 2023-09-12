enforce.rule1<-function( G, intervened, sign=NULL, sepset=NULL, verbose=0 ) {
  # Enforces rule 1, aka skeleton rules.
  #INPUT:
  # G          - Skeleton from the pc_skeleton function.
  # intervened - Vector of variable indexes where intervened.
  # sign       - Significances for printing, from pc_skeleton call.
  # sepset     - Separating sets for printing, from pc_skeleton call.
  # verbose    - Amount of printing, 0 for none.
  #OUTPUT:
  # An list with
  # Z - Which elements of B should be zero, integer matrix.
  # E - Which elements of Ce should be zero

  Z <- 1*(G == 0)
  diag(Z)<-0 #central is not set to zero
  Z[intervened,]<-0 #arcs to intervened are not marked as zeros

  E<-1*(G == 0)
  diag(E)<-0 #central 
  E[intervened,]<-0 #any covariances where an intervened variable is, are not set to zero
  E[,intervened]<-0 

  if ( verbose > 0 ) { #this is just for printing
    for ( i in 1:nrow(Z) ) {
      for ( j in 1:ncol(Z) ) {
        if ( Z[i,j] > 0 ) {
          cat('B[',i,',',j,']=0 by R1.')
          if ( !is.null(sign) &&  !is.null(sepset) ) {
            cat('', i,'_|_',j,'|(',sepset[[i]][[j]],') by ',sign[[i]][[j]])
          }
          cat('\n')
       }
      }
    }
  }

  list(Z=Z,E=E)
}
