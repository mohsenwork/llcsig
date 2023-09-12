createDataSet<-function(M,type='single',total_samples=10000) {
  #Creates a dataset.
  #INPUT:
  # M - model from randomModel().
  # type - type of experiments 
  #     'single' - single intervention experiments + passive observational data
  #     'random' - random experiments (not including a passive observational data)
  #                the number of experiments is drawn uniformly from 1 to n
  #
  #total_samples - Total number of samples divided among the experiments. Can be Inf.
  #OUTPUT:
  # list structure of the data. Each list item is an item from createData().
  n<-nrow(M$B)

  if ( type == 'single') {
    nexp<-n+1
    E<-rbind(diag(n),rep(0,n))
  } else if ( type == 'random' ) {
    nexp<-sample(n,1)
    E<-array(NA,c(nexp,n))
    for ( i in 1:nexp) {
      E[i,]<-dec.to.bin(sample(2^n-1,1),n )
    }
  } else {
    cat('Bad type of experiments given!\n')
    return(NA)
  }
  
  if ( is.infinite(total_samples) ) {
    samples<-rep(Inf,nexp)
  } else {
    #dividing samples evenly among the experiments
    samples<-divideSamples(total_samples,nexp)
  }
  D<-list()
  for ( i in 1:nexp) {
    D[[i]]<-createData(M,E[i,],samples[i])
  }

  D
}