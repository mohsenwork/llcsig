Null<-function (M,tol=1e-10) {
  S <- svd(M,nv=ncol(M))
  values<-S$d
  #values[(min(dim(M))+1):length(values)]<-0
  non_zero<-  abs(values) > tol*max(abs(values) ) #max(dim(M))*
  sum_non_zero<-sum( non_zero)
  #cat('Singular vals:', values, '\n' )
  #cat('Non-zero singularvals:', sum_non_zero, 'min_nonzero=', min(c(abs(values[non_zero]),Inf)) ,
   #    'maxzero=',max(c(abs(values[!non_zero]),-Inf)) ,'\n' )
  #browser()
  N <- S$v[,index((sum_non_zero+1),ncol(M)), drop = FALSE]

  N
}