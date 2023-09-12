addEquations<-function( e, Cmut, eqs=NULL ) {
  #Function adding regular equations to the equation structure.
  #INPUT:
  # e       - Vector of 1 where intervened. 
  # Cmut    - Covariance matrix of the mutilated network.
  #         - Can be replaced by Emut, total effects matrix.
  # eqs     - List structure of equations. If NULL it is created.

  #OUTPUT:
  # list with
  #  $n    - The number of nodes in the network.
  #  $nexp - Experiments added in total.
  #  $P    - Indexing structure for indexing b -vector and matrix B. Used
  #          by lin2pair and pair2lin.
  #  $K    - K matrix of the paper.
  #  $k    - k vector of the paper.
  #  $sat  - A matrix where sat[i,j] is TRUE if pair condition (j,i) is satisfied.

  if ( is.null(eqs) ) {
    eqs<-list()

    eqs$n <- length(e) # Total number of observed variables
    eqs$nexp <- 0 # Number of experiments added
  
    eqs$P <- matrix(0,2,0) #this is an indexing structure btw. b vector and B matrix
    for (i in 1:eqs$n) {
      for (j in 1:eqs$n) {
        if (i==j) next; #skipping the diagonal element
        eqs$P <- cbind(eqs$P,c(i,j))
      }
    }
  
    # Initialize (no rows yet, will get them with cbind)
    eqs$K <- matrix(0,0,eqs$n*(eqs$n-1))
    eqs$k <- matrix(0,0,1)

    eqs$sat<-array(FALSE,c(eqs$n,eqs$n)) # which pair conditions are satisfied
  }

  # Go through all the experiments
  eqs$nexp <- eqs$nexp + 1 

  # Get the experimental indices and the data
  evec <- (e > 0) #logical array of which elements were intervened on
  eveci <- which(evec) #vector of indexes where intervened
  evecni <- which(!evec) #vector of indexes where observed

  for (j in eveci) {
    for (u in evecni) {

      #right side of the equation
      b <- Cmut[u,j] #this is t(j->u||J)

      a <- rep(0,eqs$n*(eqs$n-1))

      #setting the first
      a[pair2lin(eqs$P,u,j)] <- 1 #because in the equation b(j->u) has coefficient 1

      for (k in evecni) { #for nodes that where observed
          if ( k != u ) {
            a[pair2lin(eqs$P,u,k)] <- Cmut[k,j] # t(j->k||J)*b(k->u) 
          }
      }

      eqs$K <- rbind(eqs$K,a)
      eqs$k <- rbind(eqs$k,b)

      eqs$sat[u,j]<-TRUE #marking the pair (j,u) satisfied
    }
  }

  eqs
}