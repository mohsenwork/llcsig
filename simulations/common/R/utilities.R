# D-separation for cyclic mixed graphs code from paper:
# "Constraint-Based Causal Discovery: Conflict Resolution with Answer Set Programming"
# By Antti Hyttinen, Frederick Eberhardt, Matti Järvisalo
# The code below has the same logic as the original, with some non-functional changes.
# Keep in mind that indexing in Python starts at zero while indexing in R starts at one.


directed_reachable<-function(x, y, C, J, D, B, U) {
  #Implements a d-separation oracle.
  #x,y    variables considered
  #C      vector of variable indexes in the conditioning set
  #J      vector of variable indexes in the intervention set  
  #D      binary matrix of directed edges, G[i,j] = 1 means j->i is present
  #B      symmetric binary matrix of bidirected edges <-> , for latent confounding
  #U      symmetric binary matrix of undirected edges --- , should be empty at beginning

  #Note that this is not currently the most efficient implementation.
  #-One could use matrix operations
  #-One could calculate more of the relations all at once.
    
  n<-nrow(D)

  D[J,]<-0
  B[J,]<-0
  B[,J]<-0
  
  ###############################################################
  
  HH<-B #symmetric head head paths
  
  if ( is.null(U) ) { #there should generally be no head-head paths in the beginning 
    TT<-array(0,c(n,n)) #this is for tail tail paths
  } else {
    TT<-U
  }
  TH<-D # #notice here TH[x,y] = 1 iff path y->x
  
  #TH and HH self paths do not affect the d-connectedness so they can be ignored
  diag(TH)<-(-1)
  diag(HH)<-(-1)  
  
  for ( node in 1:n ) {    #doing either conditioning or marginalizing to all variables
    
    if (node == x || node == y) next;#skip variables that are in the final graph
    
    #gather up all different kinds of parents, or connected nodes
    thpa<-which( TH[node,] == 1)
    htpa<-which( TH[,node] == 1) #aka children
    hhpa<-which( HH[node,] == 1)
    ttpa<-which( TT[node,] == 1)
    
    if ( !(node %in% C) ) {
      #the marginalization operation is more difficult
      #i-->node-->j
      for ( i in thpa ) {
        for ( j in htpa ) {
          TH[j,i]<-1  
        }
      }
      
      #i--> node --- j
      for ( i in thpa ) {
        for ( j in ttpa ) {
          TT[i,j]<-TT[j,i]<-1  
        }
      }
      #i --> node <-- is not ok
      #i --> node <-> is not ok
      #####################################
      
      #i <-> node --> j
      for ( i in hhpa ) {
        for ( j in htpa ) {
          HH[j,i]<-HH[i,j]<-1  
        }
      }
      
      #i <-> node --- j
      for ( i in hhpa ) {
        for ( j in ttpa ) {
          TH[i,j]<-1  
        }
      }
      #i <-> node <-- is not ok
      #i <-> node <-> is not ok
      ######################################
      
      # i --- node --- j
      #tail tail parents connected
      for ( i in ttpa ) { #connects node to itself as well so tt-self cycle is inherited
        for ( j in ttpa ) {
          TT[i,j]<-TT[j,i]<-1  
        }
      }
      #i --- node --> j
      for ( i in ttpa ) { #connects node to itself as well so tt-self cycle is inherited
        for ( j in htpa ) {
          TH[j,i]<-1  
        }
      }    
      # i --- node <-> j done already
      # i --- node <-- j done already
      ##############################################
      
      # i <-- node --> j
      for ( i in htpa ) { #connects node to itself as well so tt-self cycle is inherited
        for ( j in htpa ) {
          HH[i,j]<-HH[j,i]<-1  
        }
      }    
      #i <-- node <-> j done already
      #for ( i in htpa ) { #connects node to itself as well so tt-self cycle is inherited
      #  for ( j in hhpa ) {
      #    HH[i,j]<-HH[j,i]<-1  
      #  }
      #}    
      # i --- node <-> j done already
      # i --- node <-- j done already
      
      
    } #if node not in C
    if ( node %in% C || TT[node,node] == 1 ) {
        #notice the simplicity here!
        #an unconditioned node with a selfloop actually allows through all traffic!!!
      #only three options
      # i--> node <--j
      for ( i in thpa ) {
        for ( j in thpa ) { #notice that this connects that parents to them selves as well
          TT[i,j]<-TT[j,i]<-1  
        }
      } 
      # i<-> node <->j
      #hh parents need to be connected by head head nodes
      for ( i in hhpa ) {
        for ( j in hhpa ) {
          HH[i,j]<-HH[j,i]<-1  
        }
      } 
      # i<-> node <--j
      #connecting hh parent to th parent
      for ( i in hhpa ) {
        for (j in thpa ) {
          TH[i,j]<-1
        }
      }
    } #if node in C or a tail to tail path
    #only tailtail cycles are relevant to d-connection
    diag(TH)<-(-1)
    diag(HH)<-(-1)
    
    #now take the node away
    TH[node,]<-TH[,node]<-HH[,node]<-HH[node,]<-TT[,node]<-TT[node,]<-(-1)
  }
  #the nodes are connected if any of the paths is present in the end, where
  #all variables 
  HH[x,y] ==1 | TH[x,y] ==1 |TH[y,x] == 1 | HH[y,x] == 1 | TT[y,x] == 1
}