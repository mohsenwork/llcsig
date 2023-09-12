remove.selfcycles<-function(M) {
  #Removes all self cycles from a Model.
  #INPUT
  # M - model
  #OUTPUT
  # M2 - Model with self cycles removed
 
  M2<-M
  n<-nrow(M$B)

  A<-(diag(n)-M$B)
  iA<-mpinv(A)

  #first take out the self loops one by one.
  for ( i in 1:n ) {
    Ui<-array(0,c(n,n))
    Ui[i,i]<-1
    M2$B<-M$B-M$B[i,i]/(1-M$B[i,i])*Ui%*%A
  }

  #then fixing the covariance matrix
  A2<-diag(n)-M2$B
  M2$Ce<- A2%*%iA%*%M$Ce%*%t(iA)%*%t(A2) 

  M2
}