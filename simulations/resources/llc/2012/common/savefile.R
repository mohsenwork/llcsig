savefile<-function(M,file,append=FALSE) {
  #Saves a matdix M to a file.
  write.table(M,file=file,append=append,col.names=FALSE,row.names=FALSE)
}