resampleData<-function(D) {
  #Resamples the data in D.
  #INPUT:
  # either a data set from createDataSet() or data from createData().
  #OUTPUT:
  # Data as from createData() or Dataset as from createDataSet()

  if ( is.list(D[[1]])) {
    for ( i in 1:length(D) ) {
      D[[i]]<-resampleData(D[[i]])
    }
  } else {
    #creating a resampling version of D
    n<-nrow(D$Cx)
    if ( is.infinite(D$N) ) {
      #cat('Warning trying to resample an infinite dataset\n')
      return(D)
    }
    D$Cx<-cov(D$data[sample(D$N,D$N,replace=TRUE),])

    #making the data unusable, resampling resampled data is nonsense
    D$data<-NA
  }

  D
}
