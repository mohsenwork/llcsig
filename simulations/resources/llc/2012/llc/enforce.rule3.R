enforce.rule3<-function(Cmut, N, intervened, observed, pind, pdep, verbose=0 ) {
  # Enforces rule 3, aka Experimental covariance rule.
  #INPUT:
  # Cmut           - Covariance matrix from the experiment
  # N              - The number of samples used for creating Cmut. For statistical tests.
  # intervened     - Vector of variable indexes where intervened.
  # observed       - Vector of variable indexes where not intervened but observed.
  # pind           - If the test p-value is above this, variables are independent.
  # pdep           - If the test p-value is below this, variables are dependent.
  # verbose        - Amount of printing, 0 for none.
  #OUTPUT:
  # An list with
  # Z - Which elements of B should be zero, integer matrix.
  # E - Which elements of Ce should be zero
  #experimental covariance rule

  Z<-array(0,dim(Cmut) )
  E<-array(0,dim(Cmut) )

  #preparing the data for use with a test

  for ( x in intervened ) {
    for ( y in observed ) {
      for ( z in observed ) {
        if ( y != z ) {
          #x is in J, y,z in U
          #if t(x->y||J)!=0 and
          #   x _|_ z || y 
          # then b(z->y)=0 and the confounding covariance is zero

          p1<-isIndependent(Cmut, x, y, given=NULL, N=N )
          p2<-isIndependent(Cmut, x, z, given=y, N=N )

          if ( p1 < pdep  &&  p2 > pind  ) {
            Z[y,z]<-Z[y,z]+1
            E[y,z]<-E[y,z]+1
            E[z,y]<-E[z,y]+1

            if ( verbose > 0 ) {
              cat('B[',y,',',z,']=0 by R3.')
              cat(' t(',x, '->', y,')!=0 by', p1, '&')
              cat('', x, '_|_', z,'|',y,'by', p2, '\n')
            }#if verbose
          }#if dep & ind
        } #if y!= z
      }# for z
    }# for y
  }# for x
  list(Z=Z,E=E)
}