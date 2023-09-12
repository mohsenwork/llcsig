cyclic<-function(G) {
  #Returns true iff the graph G is cyclic.
  A<-accessible(G)
  any( diag(A) != 0 )
}