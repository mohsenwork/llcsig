ROC<-function(Z, golden, file=NA) {
  #Function drawing the ROC curve and returning the area under it.
  #INPUT:
  # Z  - Matrix of Z scores where Z[i,j] is the score of edge j -> i.
  # golden - Golden standard graph, golden$G[i,j]=1 if j->i exists in the true graph.
  #          Either from dream.golden or from randomModel.
  #file - file name to draw the graph into. NA if no curve plotted.
  #OUTPUT:
  #AUROC - Area under the ROC curve.

  #Make sure we have the absolute values of the scores.

  Z<-abs(Z)
  #Z<-Z#+1e-10*rnorm(length(Z))#adding some noise to fix the order

  n<-nrow(golden$G)

  #indexing matrix to help to deal with the diagonal
  IND<-array(TRUE,c(n,n))
  diag(IND)<-FALSE

  z<-as.vector(Z[IND])
  g<-as.vector(golden$G[IND])

  #thresold<-seq(from=min(z)-0.1,to=max(z)+0.1,by=0.1)
  threshold<-c(-Inf,sort(z),Inf) #use these thresolds
  TPR<-rep(NA,length(threshold))
  FPR<-rep(NA,length(threshold))
  
  #calculate the rates for different values of the threshold
  for ( i in 1:length(threshold) ) {
    TP<-sum( g[ z >= threshold[i] ] == 1 ) #true positive, correct edge
    FP<-sum( g[ z >= threshold[i] ] == 0 ) #false positive, additional edge
    TN<-sum( g[ z < threshold[i] ] == 0 ) #true negative, correct zero
    FN<-sum( g[ z < threshold[i] ] == 1 ) #false negative, missed edge

    TPR[i]<- TP/(TP+FN) #the rates
    FPR[i]<- FP/(FP+TN)
  }

  if ( !is.na(file) ) { #plot the curver if wanted.
    pdf(6,6,file=file,title='ROC')
    par(lwd=3,mai=c(0.8,0.8,0.1,0),cex=1.25)
    plot(c(0,1),c(0,1),xlim=c(0,1),ylim=c(0,1),col='black',type='l',
      ylab='', xlab='',main='',lwd=3,axes=FALSE,lty=2,cex=5) #this is the line
    box()
    lines(FPR,TPR,col='red')
    axis(1,at=c(0,.25,.5,.75,1),labels=c('0','.25','.5','.75','1.0'),lwd.ticks=3,lwd=3)
    axis(1,at=c(.5),labels=c('FPR'),lwd.ticks=0,lwd=0,line=1)
    axis(2,at=c(0,.25,.5,.75,1),labels=c('0','.25','.5','.75','1.0'),lwd.ticks=3,lwd=3)
    axis(2,at=c(.5),labels=c('TPR'),lwd.ticks=0,lwd=0,line=1)
    dev.off()
  }

  #calculate Area Under Curve and return it!
  AUC<-sum(abs(diff(FPR))*TPR[2:length(TPR)])
    
  AUC
}