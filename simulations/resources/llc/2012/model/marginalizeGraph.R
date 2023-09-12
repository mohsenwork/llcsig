 marginalizeGraph<-function(G, observed ) {
#Marginalized a graph containing the direct causal relationships, such that causal relationships are preserved.
# x->z->y => x->y, x<-z->y => x y
#INPUT:
# G   - original graph
# observed  - vector of truth value wheter a variable is observed
#OUTPUT:
# G - graph with some variables marginalized
#

   while( any( !observed ) ) {
     m<-which(!observed)[1]
     for ( i in which(G[m,]==1) ) { #edge from i to m
       for ( j in which(G[,m]==1) ) { #edge from m to j
	  if ( j != i ) G[j,i]<-1 #add edge
	}
     }
     G<-G[-m,-m] #now we can take out the variable
     observed<-observed[-m]
   }
   G
 }

