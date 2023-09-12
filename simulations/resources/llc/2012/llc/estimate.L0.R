estimate.L0<-function(eqs, maxpa=3, lambda=0.01) {
  #L0 regularized solution of a linear equations system.
  #INPUT:
  # eqs - System of equations eqs$K%*%b =eqs$k.
  # maxpa - Maximum number of parents considered.
  # lambda - Regularization penalizing the number of parents.
  #OUTPUT:
  # b - solution vector or B directs effects matrix if 
  #     indexing structure eqs$P given

  n<-eqs$n

  #insert the solution with maxpa=0
  best.SSE<-sum(eqs$k^2)
  b<-rep(0,ncol(eqs$K))

  for (npa  in 1:maxpa ) {
    pasets<-combinations(ncol(eqs$K),maxpa) 

    for ( ipa in 1:nrow(pasets) ) {
      #sk<-eqs$k[index.k]
      sK<-eqs$K[,pasets[ipa,],drop=FALSE]
      estimate<-mpinv(t(sK)%*%sK)%*%t(sK)%*%eqs$k

      #evaluation criterion
      SSE<-sum( ( sK%*%estimate - eqs$k )^2 ) + lambda*npa


      if ( SSE < best.SSE ) {
	print(SSE)

	#print(paset)
	best.SSE<-SSE
	b<-rep(0,ncol(eqs$K))
	b[pasets[ipa,]]<-estimate
      }
    }
  }
  if ( all( !is.na( eqs$P ) ) ) {
    b.to.B(b, eqs$P, eqs$n )
  } else {
    b
  }
}
# dream.estimate.L0<-function(eqs, maxpa=5, lambda=0.01) {
#   n<-eqs$n
# 
#   b<-rep(0,ncol(eqs$K))
# 
#   pasets<-array(NA,c(2^n,n))# not like this for the 100 node data!
#   for ( i in 1:2^n) {
#     pasets[i,]<-dec.to.bin(i-1,n)
#   }
# 
#   for ( i in 1:n ) {   #iterate through the variables
# 
#     #take indexes of b corresponding to b( -> i )
#     index.b <- eqs$P[1,] == i
# 
#     index.k <- rowSums( abs(eqs$K[,index.b]) > 1e-10 )!=0 #the equations that affect i
# 
#     sk<-eqs$k[index.k]
#     sK<-eqs$K[index.k,index.b]
# 
#     best.SSE<-Inf
#     estimate<-NA
# 
#     for ( j in 1:nrow(pasets) ) {
#       paset<-pasets[j,]
#       if ( sum(paset) > maxpa) next #no more that the given number o
#       if ( paset[i] == 1 ) next #not considering sets where the node
# 
#       #print(paset)
#       if ( sum(paset) == 0 ) {
# 	#ssK<-sK[,which(paset[-i]==1)]
# 	estimate<-c()
# 	SSE<- sum( sk^2 ) + lambda*sum(paset)
#       } else {
# 	ssK<-sK[,which(paset[-i]==1),drop=FALSE]
# 	estimate<-mpinv(t(ssK)%*%ssK)%*%t(ssK)%*%sk
# 	SSE<- sum( ( ssK%*%estimate - sk )^2 ) + lambda*sum(paset)
#       }
#       if ( SSE < best.SSE ) {
# 	#print(paset)
# 	best.SSE<-SSE
# 	apu<-rep(0,n-1)
# 	apu[which(paset[-i]==1)]<-estimate
# 	b[index.b]<-apu
#       }
#     }               
#   }
#   if ( all( !is.na( eqs$P ) ) ) {
#     b.to.B(b, eqs$P, eqs$n )
#   } else {
#     b
#   }
# }