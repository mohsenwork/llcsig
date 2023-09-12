test.stability<-function( M,  maxabs=0.9, verbose=1 ) { #eigen=TRUE,

 #Tests asymptotic stability of the model.
 #INPUT:
 # M - canonical model
 # maxabs - thresold for maximum absolute eigenvalue
 # verbose - 0 no printing, 1 nonstable printing, 2 stable printing

 #this just tests the invertibility
  #invstable<-TRUE
  eigenstable<-TRUE

  n<-nrow(M$B)
  #eye<-diag(n)

  for ( i in 0:(2^n - 1 ) ) {
    Ei<-dec.to.bin(i,n)

# old code for testing invstability
#    E<-experiment.canonical(M,Ei) #using here this function, faster option is 
                                  #to use condition number somehow straight up
                                  #without actually inverting the matrix
#     if ( any(is.na(E)) ) {
#       invstable<-FALSE
#       if ( verbose ) { cat('NOT INVSTABLE WITH INT. SET={',which(Ei==1),'}\n')
#       #browser()
#     } else if ( printstable ) {
#       cat('Invstable with int. set={',which(Ei==1),'}\n')
#     }
    U<-diag(1-Ei)
    UB<-U%*%M$B

    e<-eigen(UB,only.values=TRUE)
    #here could only get the abs max value by power method or something

    if ( any(abs(e$values) >= maxabs ) ) {
      eigenstable<-FALSE
      if ( verbose >= 1 ) {
        cat('NOT ASYMPTOTICALLY STABLE WITH INTERVENTION SET={',which(Ei==1),'} ',e$values,'\n')
      }
    } else {
      if ( verbose >= 2 ) {
        cat('Asymptotically stable with intervention set={',which(Ei==1),'} ',e$values,'\n')
      }
    }
  } #for i

  eigenstable
}

