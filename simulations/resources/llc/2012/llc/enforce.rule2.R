enforce.rule2<-function( Cmut, N, intervened, observed, pind, pdep, verbose=0 ) {
  # Enforces rule 2, aka P-rule.
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
  # E - Which elements of Ce should be zero (Non are marked in this function.)

  Z<-array(0,dim(Cmut) )

  for ( i in intervened ) {
    for ( o1 in observed ) {
      for ( o2 in observed ) {
        if ( o1 != o2 ) {
          p1<-isIndependent(Cmut, i, o1, given=NULL, N=N ) # i and o1 are independent
          p2<-isIndependent(Cmut, i, o2,  given=NULL, N=N ) # i and o2 are dependent 

          if (  p1 > pind  &&  p2 < pdep  ) {
            #then there is no arc from o2 to o1
            Z[o1,o2]<-Z[o1,o2]+1
            if ( verbose > 0 ) {
              cat('B[',o1,',',o2,']=0 by R2.')
              cat(' t(',i, '->', o1,')=0 by', p1, '&')
              cat(' t(', i, '->', o2,')!=0 by', p2, '\n')
            }
          } #if ind && dep
        } #o1 != o2
      } #o2
    } #o1
  } #i

  list(Z=Z,E=array(0,dim(Z)))
}