index<-function(from,to) {
  #matlab type of indexing
  if ( from > to ) {
    R<-c()
  } else {
    R<-from:to
  }
  R
}