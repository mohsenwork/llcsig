perturb<-function( B, i, j, delta ) {
  #Perturbs the direct effects matrix B according to pair (i,j) with delta.
  #INPUT:
  # B - direct effects matrix
  # i,j - pair with respect to which to do the perturbation.
  # delta - delta parameter
  #OUTPUT:
  # B - perturbed direct effects B-matrix.

  n<-nrow(B)

  hB<-B
  hB[j,i]<- B[j,i] + delta

  k<-(-c(i,j))
  hB[j,k]<-B[j,k]-delta*(B[i,k]+B[i,j]*B[j,k])/(1-B[i,j]*B[j,i])

  hB
}