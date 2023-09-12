b.to.B<-function(b, P, n) {
  #Function for getting from vector b to matrix B.
  #INPUT:
  # n    - The number of nodes in the network.
  # P    - Indexing structure for indexing b -vector and matrix B. Used
  #        by lin2pair and pair2lin.
  #OUTPUT:
  # B      - Vector b rearranged as matrix B.

  if ( any(is.na(b) ) ) return(NA)

  B <- matrix(0,n,n)

  for (k in 1:dim(P)[2]) {
    ij <- lin2pair(P,k)
    B[ij[1],ij[2]] <- b[k] #just putting into place
  }

  B
}