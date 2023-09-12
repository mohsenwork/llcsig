

pc_skeleton <- function(C, N, significance=.05, maxsetsize=-1) {
  # Runs the first stage of PC, finds the skeleton.
  #
  # THIS CODE IS INSPIRED BY AND LOOSELY BASED ON PROF. SPIRTES'S PC CODE
  # THAT WAS USED WITH PERMISSION IN THE CODE OF UAI09:"Bayesian Discovery
  # of linear Acyclic Causal Models".
  #
  #INPUT:
  # C            - Covariance matrix of the data.
  # N            - The number of samples used for creating C.
  # significance - Significance level for the tests.
  # maxsetsize   - The maximum size of a set that is used to condition on.
  #                -1 for no restrictions.
  #OUTPUT:
  # An list with
  # $skel        - Skeleton. 1 where edges.
  # $sepset      - List of separating sets.
  # $sign        - List of significances, where variables where independent.

  n <- nrow(C)
  skel <- array(1,c(n,n))
  diag(skel) <- 0

  #keeping track of separating sets and significance levels
  sepset <- vector(length = n, mode = "list")
  sig <- vector( length = n, mode = "list") #added by antti
  for (i in 1:n) {
    sepset[[i]] <- vector(length = n, mode = "list")
    sig[[i]] <- vector( length = n, mode = "list")
  }

  #here the first marginal test
  for (x in 1:(n-1)) {
    for (y in (x+1):n) {
      pval <- isIndependent(C, x, y, given=NULL, N=N)
      if ( pval > significance ) {
          skel[x,y] <- skel[y,x] <- 0
          #separated by the emptyset, also pval saved
          sepset[[x]][[y]] <- sepset[[y]][[x]] <- numeric(0)
          sig[[x]][[y]] <- sig[[y]][[x]] <- pval
      }
    }
  }

  setsize <- 1
  while ( ( is.na(maxsetsize) || maxsetsize == -1 || setsize <= maxsetsize ) && #only go to certain setsizes
          max(rowSums(skel)) >= setsize ) { #as long as some row has at least i possible neighbours

    for (x in 1:n) { #going through all nodes
      if ( sum(skel[x,]) <= setsize)  next #dont consider nodes that have less than i possible neighbours

      for (y in which(skel[x,] == 1) ) { #y is one of the neighbours
        if (sum(skel[x,]) <= setsize) break
        # The combn function is obnoxiously inconsistent.  It behaves
        # differently when the first argument has length one and
        # produces a non-matrix output when there is only one combination
        # in the output.

        tf <- setdiff(which(skel[x,] == 1), y) #taking y out of the neighbours of x

        if (length(tf) == setsize) { 
          sepsets <- matrix(tf) 
        } else {
          sepsets <- combn(tf, setsize)
        }

        for ( i in 1:ncol(sepsets) ) {
          pval <- isIndependent(C, x, y, given=sepsets[,i], N=N )
          if ( pval > significance ) {
              skel[x,y] <- skel[y,x] <- 0

              #saving the separating set and significance level
              sepset[[x]][[y]] <- sepset[[y]][[x]] <- sepsets[,i]
              sig[[x]][[y]] <- sig[[y]][[x]] <- pval

              break
          }#if
        }#for i
      }#for y
    }#for x
    setsize = setsize + 1 #consider next one element bigger sets
  }#while 

  list(skel=skel, sepset=sepset, sign=sig)
}

