loadfile<-function(file) {
  #Loads matrix from a file.
  as.matrix(read.table(file=file))
}