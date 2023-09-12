pair2lin <- function(P,i,j) {
  #indexing function for switching btw. B and b.
  #INPUT:
  # P - indexing structure (often eqs$P)
  # i,j - index (row,col) of the matrix B
  #OUTPUT:
  # ind index of b
  if (i == j) stop('i and j should not be equal!')
  ind <- which((P[1,]==i) & (P[2,]==j))
  ind
}
