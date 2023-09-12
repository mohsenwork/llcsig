llc_bootstrap<-function(DATA, n, maxsetsize, rules, penalty='none', lambda=0, alpha=0.05, bootstraps=30 ) {
  #Runs the bootstrapped LLC algorithm, returning also the z-scores. Basically bootstraps "llc_plain".
  #INPUT:
  # DATA        - data in the from from createDataSet.
  # n           - number of nodes
  # penalty     - 'none', 'L0', 'L1' or 'L2' - type of penalization/regularization used
  #               in the solving of the linear equation systems
  # lambda      - regularization parameter
  # alpha       - specifies the significance level for the conditional independence tests
  # bootstraps  - number times to run llc_plain

  #OUTPUT:
  # list with
  #     $B      - Mean of the estimated direct effects over the bootstraps.
  #     $Bsd    - Standard deviation of the estimated direct effects over the bootstraps.
  #     $Ce     - Mean of estimated elements of the covariance matrix.
  #     $Cesd   - Standard deviations of estimated elements of the covariance matrix.
  #     $Bz     - Z-scores for the elements of B. Absolute Z score is a measure of confidence on the existence of the edge.
  #     $Cez    - Z-scores for the elements of Ce. Absolute Z score is a measure of confidence on the existence of the edge.
  #     $PC     - Which pair conditions were satisfied by the experiments?
  #     $COV    -  Which cov conditions were satisfied by the experiments? (symmetric matrix)
  #     $Bcon   - Conservative estimate of which elements of B are identified. See article for details.
  #     $Bnull  - Norms of the parameters projected to the null-space. (As means of the norms over the bootstraps).
  #               Values close to 0 mean that the element of B is identified. See article for details.
  #     $Cecon  - Conservative estimate of which elements of B are identified. See article for details.


  # n<-nrow(DATA[[1]]$Cx)
  # We dont actuall need Cx, just the number of nodes.

  B<-array(NA, c(n,n,bootstraps) )
  Ce<-array(NA, c(n,n,bootstraps) )

  Bnull<-array(NA, c(n,n,bootstraps) )


  for ( b in 1:bootstraps ) {
    #resampling the data, aka bootsrrapping
    D<-resampleData(DATA)


    R<-llc_plain(D,maxsetsize=maxsetsize,rules=rules,penalty=penalty,lambda=lambda,null=FALSE, alpha=alpha)
                                                    #do not request the identifiability information

    #saving the results for which we want to examine the distribution
    B[,,b]<-R$B
    Ce[,,b]<-R$Ce
    Bnull[,,b]<-R$Bnull
  }

  #calculation of the different means
  Bmean<-apply(B,c(1,2),mean)
  Bsd<-apply(B,c(1,2),sd)

  Cemean<-apply(Ce,c(1,2),mean)
  Cesd<-apply(Ce,c(1,2),sd)

  Bnullmean<-apply(Bnull,c(1,2),mean)

  
  #Some results appear in R, which is just the result of the last bootstrap round.
  
  #Note the calculation of the Z-scores here. 
  
  list(B=Bmean,Bsd=Bsd,Ce=Cemean,Cesd=Cesd,Bz=Bmean/(Bsd+1e-10),Cez=Cemean/(Cesd+1e-10),
       PC=R$PC,COV=R$COV,Bcon=R$Bcon,Bnull=Bnullmean,Cecon=R$Cecon)
}