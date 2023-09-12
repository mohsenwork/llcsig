llc_plain<-function(DATA, maxsetsize, rules, penalty='none', lambda=0, null=FALSE, alpha=0.05) {
  #Runs the plain LLC algorithm, close to what is given in the article.
  #INPUT:
  # DATA          - data in the from from createDataSet.
  # maxsetsize    - The maximum size of a set that is used to condition on when using faithfulness rules.
  #                 -1 (or NA) for no restrictions.
  # rules         - which of the faithfulness rules (1,2,3) should be used.
  # penalty       - 'none', 'L0', 'L1' or 'L2' - type of penalization/regularization used
  #                 in the solving of the linear equation systems
  # lambda        - regularization parameter
  # null          - specifies whether to request identifiability information about the null-space. Use null=True to request such information.
  # alpha         - specifies the significance level for the conditional independence tests

  #OUTPUT:
  # list with
  #     $B                - Estimate for the direct effects.
  #     $Ce               - Estimate of the covariance matrix. May contain NA if cov-condition is not satisfied for all pairs.
  #     $PC               - Which pair conditions were satisfied by the experiments?
  #     $COV              - Which cov conditions were satisfied by the experiments? (symmetric matrix)
  #     $Bcon             - Conservative estimate of which elements of B are identified. See article for details.
  #     $Bnull            - Norms of the parameters projected to the null-space. (only generated when null=True)
  #                          Values close to 0 mean that the element of B is identified. See article for details.
  #     $Bvar             - Posterior variances of the parameter estimates. More robust way of determining
  #                         the identified parameters. (only generated when null=True)
  #                         Values close to 0 mean that the element of B is identified. See article for details.
  #     $Cecon            - Conservative estimate of which elements of B are identified. See article for details.
    
  # Runs plain LLC algorithm
    n<-ncol(DATA[[1]]$Cx)
    #parameters for pc and ind and dep thresolds
    pc_significance <- 0.05
    pc_maxsetsize <- maxsetsize
    pind <- alpha
    pdep <- alpha

    #intialize the equation systems
    eqs<-list()
    for ( node in 1:n ) {
      eqs[[node]]<-list(n=n,P=NA)          # named list with attributes n = number of nodes and P = NA
      eqs[[node]]$K<-array(0, c(0,n-1) )   # Matrix of n-1 columns with 0 rows. shape is c(0,n-1)
      eqs[[node]]$k<-c()
    }

    #get the weigths
    #Each experiment has a weight
    if ( is.infinite(DATA[[1]]$N) ) { #handling infinite sample limit data
      w<-rep(1,length(DATA))
    } else {
      w<-rep(NA,length(DATA))
      for ( d in 1:length(DATA) ) {
        w[d]<-DATA[[d]]$N
      }
      w<-w/sum(w)*length(DATA) #Otherwise weight the equations according to sample sizes
                  #What would be more optimal?
    }

    #Fulfilling the satisfied pairs when processing the data
    # PC[u,i] is true if u observed and i intervened.
    # This corresponds to the pair condition for the ordered pair (i,u) being satisfied following the notation of the paper. 
    PC<-array(FALSE,c(n,n))
    diag(PC)<-TRUE

    # COV[u,i] is true if u and i are observed.
    COV<-array(FALSE,c(n,n))
    # Those two arrays save constraints on the disturbance terms for rule 1 and 3 respectively
    E1<-array(NA,c(n,n,length(DATA)))
    E3<-array(NA,c(n,n,length(DATA)))

    # For each experimental covariance matrix
    for ( d in 1:length(DATA) ) {
      #first decorrelate the data, if the variables were not randomized independently
      #with finite number of samples this is also useful.
      #Lemma 5 in paper
      DATA[[d]]$Cx<-Cx<-decorrelate(DATA[[d]]$Cx,DATA[[d]]$e)
    
      #reading the intervention vector
      J<-which(DATA[[d]]$e==1)
      U<-which(DATA[[d]]$e==0)
      #this now works for J = []
      COV[U,U]<-TRUE

      #create the equations to the subsystem!
      #This creates a system of equations for each node of the from T_uu * b_u = t_u (notation of paper)
      #Here the notation eqs[[u]]$K * b = eqs[[u]]$k is used.
      for ( i in J ) {
        for ( u in U ) {
          PC[u,i]<-TRUE
          # COV[u,U]<-TRUE
         # kvec<-rep(0,n) #kvec[j] is the coefficient over b_uj thus t(i->j)

         # Cx[,i] is i-th column of Cx and consists of t(i->1||J), t(i->2||J), t(i->3||J) .. t(i->n||J). (see 18 in Paper)
         # w[d] is the weight of experiment d

         #DOES THIS WORK FOR OTHER THAN SINGLE INTERVENTION DATA???
          kvec<-w[d]*Cx[,i] #since decorrelation is used also puts coefficient 1*b_ui
                            #and the coefficient for any other intervened variable is zero
                            #here multiplying by the samplesize weight of the equation
                            #

          # This adds row t(i->1||J), t(i->2||J), t(i->3||J) .. t(i->n||J) to K of K * b = k
          # Since decorr is used, we already have t(i->i||J) = 1 and t(i->j||J) = 0 for j in J
          # Taking out diagonal element: e.g. u = 3, we know that b_3,3=0, so we take out t(i->3||J) from the added row above
          eqs[[u]]$K<-rbind(eqs[[u]]$K,kvec[-u] ) #taking out the diagonal element

          # This adds t(i->u||J) to k of K * b = k
          eqs[[u]]$k<-c(eqs[[u]]$k,kvec[u])
          #browser()
        }
      }

      ## Faithfulness constraints // added from 2010 code
      ## Each rule (enforce.rule1, enforce.rule2 and enforce.rule3) return Z and E, where
      ## Z - A (nxn) matrix recording which elements of B should be zero.
      ##     An entry Z[i,j] > 0 means B[i,j] = 0 by faithfulness.
      ## E - Which elements of Ce should be zero
      ##     An entry E[i,j] > 0 means Ce[i,j] = 0 by faithfulness.
      if ( length(rules) > 0 ) {
        #keeping book how many times each element is put to zero
        Z<-array(0,c(n,n))
        #keep track of E (new compared to code 2010, which does not use E)
        E<-array(0,c(n,n))

        if ( any(rules == 1 ) ) {
          pc_result <- pc_skeleton( DATA[[d]]$Cx, DATA[[d]]$N,
                                    significance=pc_significance,
                                    maxsetsize=pc_maxsetsize)

          # if ( verbose > 0 ) { cat('PC result:\n'); print(pc_result$skel) }
          rule1 <- enforce.rule1( pc_result$skel, which( DATA[[d]]$e == 1 ),
                                  sepset=pc_result$sepset, sign=pc_result$sign,
                                  verbose=0 )
          Z <- Z + rule1$Z
          # Keep in mind E[u,i] > 0 if by faithfulness Ce[u,i] = 0.
          # Save inference of 0 covariance between such pairs u,i
          E <- rule1$E
          E[ E <= 0 ] <- NA
          E[ E > 0 ] <- 0
          E1[,,d] <- E
        }

        if ( any(rules == 2 ) ) {
          # Note rule 2 never adds constraints to disturbance covariance matrix
          rule2 <- enforce.rule2( DATA[[d]]$Cx, DATA[[d]]$N, which( DATA[[d]]$e == 1 ), which( DATA[[d]]$e == 0 ),
                                  pind=pind, pdep=pdep, verbose=0 )
          Z <- Z + rule2$Z
        }

        if ( any(rules == 3 ) ) {
          rule3 <- enforce.rule3( DATA[[d]]$Cx, DATA[[d]]$N, which( DATA[[d]]$e == 1 ), which( DATA[[d]]$e == 0 ),
                                  pind=pind, pdep=pdep, verbose=0 )
          Z <- Z + rule3$Z
          E <- rule3$E
          E[ E <= 0 ] <- NA
          E[ E > 0 ] <- 0
          E3[,,d] <- E
        }

        # Here we add equations
        for ( i in 1:nrow(Z) ) {
          for ( j in (1:ncol(Z))[-i] ) {
            if ( Z[i,j] > 0 ) {
              #This is different than the code from 2010.
              #This code has a separate sytem of equations with n-1 columns for each node in (1..n),
              #whereas the 2010 code had one sytem of equations with n(n-1) columns.

              #Z[i,j] > 0 means B[i,j] = 0 (j->i) by faithfulness.
              #Add the equation bij = 0 with weight w[d] (same as non-faithfulness related equations)
              #To the system of equations belonging to node i
              #System of equations for node i is eqs[[i]]$K * b = eqs[[i]]$k
              #Add ci,1 * bi,1 + ci,2 * bi,2 + .. + ci,i-1 * bi,1-1 + ci,i+1 * bi,1+1 + .. + ci,n * bi,n = 0 
              #Where ci,l = 1 iff l=j, else 0
              #Right hand side of equation
              eqs[[i]]$k<-c(eqs[[i]]$k, 0)
              #Left hand side for equations
              a <- rep(0, eqs[[i]]$n)
              a[j] <- w[d]
              eqs[[i]]$K<-rbind(eqs[[i]]$K, a[-i])
              #Satisfy the ordered pair condition for pair (j, i)
              PC[i, j]<-TRUE
              # eqs<-addZeroEquation( eqs, i, j, w=Z[i,j] ) # original
            }
          }
        }
      }
    }

    B<-array(NA,c(n,n))
    diag(B)<-0

    ###Solving the system using different regularizations.
    for ( node in 1:n ) {
      #cat(node,'\n')
      if ( nrow(eqs[[node]]$K) == 0 ) { #no equations for this subsystem
        B[node,-node]<-0
        next
      }
      if ( penalty == 'none' ) {
        B[node,-node]<-estimate.none(eqs[[node]])
      } else if ( penalty == 'L0' ) {
        B[node,-node]<-estimate.L0(eqs[[node]], lambda=lambda)
      } else if ( penalty == 'L1' ) {
        B[node,-node]<-estimate.L1(eqs[[node]], lambda=lambda)
      } else if ( penalty == 'L2' ) {
        B[node,-node]<-estimate.L2(eqs[[node]], lambda=lambda)
      }
    }

    ############Then should determine Ce
    Ces<-array(NA,c(n,n,length(DATA)))
    for ( d in 1:length(DATA) ) {
      J<-which(DATA[[d]]$e==1)
      U<-diag(1*DATA[[d]]$e==0) #note this is a matrix
      #note that the data has been decorrelated
      A<-diag(n)-U%*%B
      #Estimate for (Ce[U,U])
      Ces[,,d]<-A%*%DATA[[d]]$Cx%*%t(A)
      Ces[J,,d]<-NA
      Ces[,J,d]<-NA
    }
    #browser()
    Ce<-array(NA,c(n,n))
    for ( i in 1:n) {
      for( j in 1:n) {
        #could use some sort of regularization here as well
        # Ce[i,j]<-sum(w*Ces[i,j,],na.rm=TRUE)/sum(w[ !is.na(Ces[i,j,])])
        # This is modifies the nonfaithful approximation of Ce (commented out above).
        # The denominator of the weighted sum is adjusted to account for inferences from faithfulness like Ce[u,i] = 0
        # Note that the numerator remains the same (adding 0).
        Ce[i,j]<-sum(w*Ces[i,j,],na.rm=TRUE)/(sum(w[ !is.na(Ces[i,j,])]) + sum(w[ !is.na(E1[i,j,])]) + sum(w[ !is.na(E3[i,j,])]))
      }
    }


    #determine a quantatity in the 
    Bnull<-array(NA,c(n,n))
    diag(Bnull)<-0 #all diagonal are determined
    Bvar<-Bnull
    Bcon<-array(NA,c(n,n)) #conservative estimate of the determined coefficients
    diag(Bcon)<-TRUE #all diagonal are determined

    for ( i in 1:n) {

      #Conservative criterion:
      #the i:th row is determined if all PCs (*,i) are determined
      Bcon[i,-i]<-all(PC[i,-i])

      #Looks at the norms of vectors projected to the null space
      #first project ei to the space spanned by columns of K
      #then if the ei - proj(ei) is ei projected to the null space
      #if it is is a zero vector, then the parameter is determined
      if ( !null || nrow(eqs[[i]]$K) == 0 ) { #no equations for this subsystem
                                             #or not requesting the identifiability information
        Bnull[i,-i]<-Inf
        Bvar[i,-i]<-Inf
        next
      }

      determ<- Null(eqs[[i]]$K,0.01)         
      Bnull[i,-i]<-sqrt( rowSums(determ^2) ) #putting in the norms

      Kpinv<-mpinv(eqs[[i]]$K)

      #the columns of this matrix are the residuals when ei is project to the space
      #spanned by K.
      #determ_alt<-diag(ncol(eqs[[i]]$K))-Kpinv%*%eqs[[i]]$K 
      #Bnorm[i,-i]<-sqrt( rowSums(determ_alt^2) ) #putting in the norms
      #Bmax[i,-i]<- apply( abs( t(Kpinv) ),2, max)

      #here putting N(0,100) prior for the parameters ( uniformative)
      #error are from N(0,0.1)
      #then calculating the posterior variances appearing at the inverse of the information matrix

      Bvar[i,-i]<-abs( diag( mpinv( 1/(2*1e-1)*t( eqs[[i]]$K ) %*% eqs[[i]]$K + 1/(2*1e2)*diag(ncol(eqs[[i]]$K) ) ))) #sqrt( colSums(eqs[[i]]$K^2) )/max(svd(eqs[[i]]$K)$d)
      #BN<-sqrt( rowSums(determ^2) )/max(svd(determ)$d)

      #browser()

    }
    #Then need to estimate which parameters of the Error Covariance are determined
    #this is directly from the article
    Cecon<-array(NA,c(n,n))
    for ( i in 1:n) {
      for ( j in 1:n) {
         Cecon[i,j]<-all(Bcon[i,]) & all(Bcon[j,]) & COV[i,j]
      }
    }

    #Return all the info.
    list(B=B,Ce=Ce,PC=PC,COV=COV,Bcon=Bcon,Bnull=Bnull,Cecon=Cecon,Bvar=Bvar)
}