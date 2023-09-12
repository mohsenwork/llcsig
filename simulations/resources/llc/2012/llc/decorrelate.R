decorrelate<-function(Cx,e) {
#Decorrelates the intervened variables. Only works if the intervened variables are
#randomized independently of the error terms of the observed variables. (Lemma 5)
#INPUT:
# Cx - Covariance matrix
# e  - Experiment vector e[i]=1 if x_i intervened, otherwise 0
  if( all(e==0) ) {
    return(Cx)
  }
  J<-which(e==1)
  U<-which(e==0)
 # browser()
  Tx<-Cx[,J,drop=FALSE]%*%mpinv( Cx[J,J,drop=FALSE] )


  Cx-Tx%*%Cx[J,J]%*%t(Tx)+Tx%*%t(Tx)
}