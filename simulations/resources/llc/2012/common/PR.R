PR<-function(Z, golden, file=NA) {
  #Function drawing the Precision-Recall curve and returning the area under it.
  #INPUT:
  # Z  - Matrix of Z scores where Z[i,j] is the score of edge j -> i.
  # golden - Golden standard graph, golden$G[i,j]=1 if j->i exists in the true graph.
  #          Either from dream.golden or from randomModel.
  #file - file name to draw the graph into. NA if no curve plotted.
  #OUTPUT:
  #AUPR - Area under the PR curve.

  #see ROC.R for further comments on the code.
  
  Z<-abs(Z)
  Z<-Z#+1e-10*rnorm(length(Z))#adding some noise to fix the order

  n<-nrow(golden$G)

  IND<-array(TRUE,c(n,n))
  diag(IND)<-FALSE

  z<-as.vector(Z[IND])
  g<-as.vector(golden$G[IND])

  #thresold<-seq(from=min(z)-0.1,to=max(z)+0.1,by=0.1)
  threshold<-c(-Inf,sort(z),Inf) #use these thresolds
  PRECISION<-rep(NA,length(threshold))
  RECALL<-rep(NA,length(threshold))

  
  for ( i in 1:length(threshold) ) {
    TP<-sum( g[ z >= threshold[i] ] == 1 ) #true positive, correct edge
    FP<-sum( g[ z >= threshold[i] ] == 0 ) #false positive, additional edge
    TN<-sum( g[ z < threshold[i] ] == 0 ) #true negative, correct zero
    FN<-sum( g[ z < threshold[i] ] == 1 ) #false negative, missed edge
    PRECISION[i]<- (TP+1e-20)/(TP+FP+1e-20) #the rates
    RECALL[i]<- TP/(TP+FN)
    #cat(TP,FP,TN,FN,PRECISION[i],RECALL[i],'\n')
  }

  if ( !is.na(file) ) {
    pdf(6,6,file=file)
    par(lwd=3,mai=c(0.8,0.8,0.1,0),cex=1.25)
    plot(c(-1),c(-1),xlim=c(0,1),ylim=c(0,1),col='black',type='l',
      ylab='', xlab='',main='',axes=FALSE,lty=2) #this is the line
    box()
    lines(RECALL,PRECISION,col='red')

    axis(1,at=c(0,.25,.5,.75,1),labels=c('0','.25','.5','.75','1.0'),lwd.ticks=3,lwd=3)
    axis(1,at=c(.5),labels=c('RECALL'),lwd.ticks=0,lwd=0,line=1)
    axis(2,at=c(0,.25,.5,.75,1),labels=c('0','.25','.5','.75','1.0'),lwd.ticks=3,lwd=3)
    axis(2,at=c(.5),labels=c('PRECISION'),lwd.ticks=0,lwd=0,line=1)
    dev.off()
  }
  #calculate AUC and return it!
  AUC<-0
  for ( i in 1:(length(RECALL)-1) ) {
     AUC<-AUC - (PRECISION[i+1]+PRECISION[i])*(RECALL[i+1]-RECALL[i])/2
  }

  #old way of calculating it
#  AUC<-sum(abs(diff(RECALL))*PRECISION[2:length(PRECISION)]) #area under the curve
    
  AUC
}