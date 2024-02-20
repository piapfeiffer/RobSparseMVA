#' Robust Sparse CCA via alternating regressions
#'
#' Function to perform Robust Sparse CCA  using alternating regressions from
#' paper by Wilms, I. and Croux, C. (2016), "Robust sparse canonical correlation analysis"
#'
#' Implementation by Wilms, I.
#'
#'
#' @param X                   An (nxp) data matrix.
#' @param Y                   An (nxq) data matrix
#' @param lambdaAseq          grid of sparsity parameters for lambdaA
#' @param lambdaBseq          grid of sparsity parameters for lambdaB
#' @param rank                number of canonical vector pairs to extract
#' @param max.iter            maximum number of iterations
#' @param conv                tolerance value for convergence
#' @param ncores              number of cores to use
#' @param tol                 tolerance parameter sparseLTS function
#' @param mode                mode argument sparseLTS function, alternating regressions fit
#' @param alpha.sparselts     alpha argument sparse LTS function
#' @param lambda.start        sparsity parameter used in inital sparse LTS fit
#' @param lambdaA.final       sparsity parameters for lambdaA used in final sparse LTS fit (to express canonical vectors in terms of original data sets)
#' @param lambdaB.final       sparsity parameters for lambdaB used in final sparse LTS fit (to express canonical vectors in terms of original data sets)
#' @param lambdamode.start    mode argument sparseLTS function, inital sparse LTS fit
#' @param lambdamode.final    mode argument sparseLTS function, final sparse LTS fit
#'
#' @references Wilms, I., & Croux, C. (2016). Robust sparse canonical correlation analysis. BMC systems biology, 10(1), 1-13.
#'
#' @returns ALPHA:            (pxr) estimated canonical vectors correpsonding to the first data set
#' @returns BETA:             (qxr) estimated canonical vectors correpsonding to the second data set
#' @returns cancors:          r estimated canonical correlations
#' @returns U_ALL:            (nxr) estimated canonical variates corresponding to the first data set
#' @returns V_ALL:            (nxr) estimated canonical variates corresponding to the second data set
#' @returns lamdbaA:          value of the sparsity parameter lambdaA
#' @returns lamdbaB:          value of the sparsity parameter lambdaB
#' @returns it:               number of iterations
#' @export
ccaAR<-function(X,Y,lambdaAseq=seq(from=1,to=0,by=-0.01),
                          lambdaBseq=seq(from=1,to=0,by=-0.01),
                          rank,max.iter=50,conv=5*10^-2,
                          ncores=1, tol =10^-8,mode="fraction",
                          alpha.sparselts=0.75,lambda.start=0.001,
                          lambdaA.final=lambdaAseq,lambdaB.final=lambdaBseq,
                          lambdamode.start="fraction",
                          lambdamode.final="fraction"){

  ### STORE RESULTS
  ALPHA_ALL<-matrix(NA,ncol=rank,nrow=ncol(X))
  BETA_ALL<-matrix(NA,ncol=rank,nrow=ncol(Y))
  U_ALL<-matrix(NA,ncol=rank,nrow=nrow(X))
  V_ALL<-matrix(NA,ncol=rank,nrow=nrow(Y))
  cancors<-matrix(NA,ncol=rank,nrow=1)
  lambdaA_ALL<-matrix(NA,ncol=rank,nrow=max.iter)
  lambdaB_ALL<-matrix(NA,ncol=rank,nrow=max.iter)
  iterations<-matrix(NA,ncol=rank,nrow=1)
  obj.it<-matrix(NA,ncol=rank,nrow=max.iter+1)

  ### START CODE

  # Sequentially extract the r canonical vector pairs
  for (i.r in 1:rank){

    if (i.r==1){ # for r=1: start from original data sets
      X_data<-X
      Y_data<-Y
    }

    # STARTING VALUE: PCA-based
    decomp<-rARPACK::eigs(SpatialNP::SCov(X_data),k=1)
    B.PRINCOMP<-X_data%*%decomp$vectors[,1]


      LASSOFIT_init<-robustHD::sparseLTS(x=Y_data,y=c(B.PRINCOMP),lambda=lambda.start,mode = lambdamode.start,intercept=F,normalize=T,initial="sparse", tol = tol,ncores=ncores,alpha=alpha.sparselts)
      B.STARTING<-matrix(LASSOFIT_init$coefficients,ncol=1)
      B.STARTING<-apply(B.STARTING,2,NORMALIZATION_UNIT)
      A.STARTING<-matrix(decomp$vectors[,1],ncol=1)


    # CONVERGENCE CRITERION
    obj.initial<-mean((X_data%*%matrix(A.STARTING,ncol=1)-Y_data%*%matrix(B.STARTING,ncol=1))^2)
    obj.it[1,i.r]<-obj.initial

    # INITIALIZE CONVERGENCE PARAMETERS
    it<-1
    diff.obj<-conv*10

    # FROM i until convergence: canonical vectors
    while( (it<max.iter) & (diff.obj>conv) ){

      # Estimating A conditional on B
        FIT.A<-RobustSparse.alternating(Xreg=X_data,Yreg=Y_data%*%B.STARTING,lambdaseq=lambdaAseq, tol = tol,ncores=ncores,mode=mode,alpha.sparselts=alpha.sparselts)
        AHAT_FINAL<-FIT.A$COEF_FINAL
        lambdaA_ALL[it,i.r]<-FIT.A$LAMBDA_FINAL


      # Estimating B conditional on A
        FIT.B<-RobustSparse.alternating(Xreg=Y_data,Yreg=X_data%*%AHAT_FINAL,lambdaseq=lambdaBseq, tol = tol,ncores=ncores,mode=mode,alpha.sparselts=alpha.sparselts)
        BHAT_FINAL<-FIT.B$COEF_FINAL
        lambdaB_ALL[it,i.r]<-FIT.B$LAMBDA_FINAL

      # Check convergence
      obj.new<-mean((X_data%*%matrix(AHAT_FINAL,ncol=1)-Y_data%*%matrix(BHAT_FINAL,ncol=1))^2)
      obj.it[it+1,i.r]<-obj.new
      diff.obj<-abs(obj.new-obj.initial)/obj.initial

      # Updated starting values
      B.STARTING<-BHAT_FINAL
      A.STARTING<-AHAT_FINAL
      obj.initial<-obj.new
      it<-it+1

    } # end while-loop

    # Number of ITERATIONS
    iterations[1,i.r]<-it

    # CANONICAL VECTORS VARIATES after convergence
    Uhat<-X_data%*%AHAT_FINAL
    Vhat<-Y_data%*%BHAT_FINAL

    # Express canonical vectors in terms of ORIGINAL DATA MATRICES

    if (i.r==1){ # FIRST DIMENSION

      # Final estimates of canonical vectors, variates and canonical correlation
      ALPHA_ALL[,i.r]<-AHAT_FINAL
      BETA_ALL[,i.r]<-BHAT_FINAL
      U_ALL[,i.r]<-Uhat
      V_ALL[,i.r]<-Vhat
      # cancors[1,i.r]<-abs(robustbase::covMcd(cbind(Uhat,Vhat),alpha=0.75,cor=T)$cor[2,1])
      cancors[1,i.r]<-abs(cor(Uhat,Vhat))

      # Deflated data matrices
      for (i in 1:ncol(X_data)){
        if (all(X_data[,i]==0)){
          X_data[,i]=0
        } else {
          if (sum(abs(Uhat)-abs(X_data[,i]))==0){
            X_data[,i]<-0
          } else {
            X_data[,i]<-robustbase::ltsReg(x=Uhat,y=c(X_data[,i]),intercept=F)$resid
          }
        }
      }

      for (i in 1:ncol(Y_data)){
        if (all(Y_data[,i]==0)){
          Y_data[,i]=0
        } else {
          if (sum(abs(Vhat)-abs(Y_data[,i]))==0){
            Y_data[,i]<-0
          } else
          {Y_data[,i]<-robustbase::ltsReg(x=Vhat,y=c(Y_data[,i]),intercept=F)$resid
          }
        }
      }

      lambdaA_FINAL<-FIT.A$LAMBDA_FINAL
      lambdaB_FINAL<-FIT.B$LAMBDA_FINAL

    } else { # HIGHER ORDER DIMENSIONS

      # A expressed in terms of original data set X

        FIT.Aorig<-RobustSparse.alternating(Yreg=Uhat,Xreg=X,lambdaseq=lambdaA.final, tol = tol,ncores=ncores,mode=lambdamode.final,alpha.sparselts=alpha.sparselts)
        ALPHAhat<-FIT.Aorig$COEF_FINAL
        lambdaA_FINAL<-FIT.Aorig$LAMBDA_FINAL

      # B expressed in terms of original data set Y

        FIT.Borig<-RobustSparse.alternating(Yreg=Vhat,Xreg=Y,lambdaseq=lambdaB.final, tol = tol,ncores=ncores,mode=lambdamode.final,alpha.sparselts=alpha.sparselts)
        BETAhat<-FIT.Borig$COEF_FINAL
        lambdaB_FINAL<-FIT.Borig$LAMBDA_FINAL

      # Final estimates of canonical vectors, variates and canonical correlation
      ALPHA_ALL[,i.r]<-ALPHAhat
      BETA_ALL[,i.r]<-BETAhat
      Uhat<-X%*%ALPHAhat
      Vhat<-Y%*%BETAhat
      U_ALL[,i.r]<-Uhat
      V_ALL[,i.r]<-Vhat
      cancors[1,i.r]<-abs(cor(Uhat,Vhat))
      #cancors[1,i.r]<-abs(robustbase::covMcd(cbind(Uhat,Vhat),alpha=0.75,cor=T)$cor[2,1])

      # Deflated data matrices: regress original data sets on all previously found canonical variates
      for (i in 1:ncol(X_data)){
        if (all(X_data[,i]==0)){
          X_data[,i]=0
        } else {
          X_data[,i]<-robustbase::ltsReg(x=U_ALL[,1:i.r],y=c(X_data[,i]),intercept=F)$resid
        }

      }

      for (i in 1:ncol(Y_data)){
        if (all(Y_data[,i]==0)){
          Y_data[,i]=0
        } else {
          Y_data[,i]<-robustbase::ltsReg(x=V_ALL[,1:i.r],y=c(Y_data[,i]),
                                         intercept=F)$resid
        }

      }
    }

  } # END FOR-LOOP

  ## OUTPUT
  out<-list(ALPHA=ALPHA_ALL,BETA=BETA_ALL,U_ALL=U_ALL,V_ALL=V_ALL,cancors=cancors,lambdaA=lambdaA_FINAL,lambdaB=lambdaB_FINAL,it=iterations)

}

## Robust Sparse CCA AUXILIARY FUNCTIONS
RobustSparse.alternating<-function(Xreg,Yreg,lambdaseq, tol,ncores,mode,alpha.sparselts=0.75){
  # AUX FUNCTION

  ### Function to perform robust sparse alternating regression

  ### INPUT
  #Xreg               : design matrix
  #Yreg               : response
  #lambdaseq          : sequence of sparsity parameters
  # tol               : tolerance parameter sparseLTS function
  # ncores            : number of cores to use
  # mode              : mode argument sparseLTS function
  # alpha.sparselts   : alpha argument sparse LTS function

  ### OUTPUT
  #COEF_FINAL         : estimated coefficients
  #LAMBDA_FINAL       : optimal sparsity parameter

  ## Standardize
  Xreg_st<-(Xreg-matrix(apply(Xreg,2,median),ncol=ncol(Xreg),nrow=nrow(Xreg),byrow=T))/matrix(apply(Xreg,2,mad),ncol=ncol(Xreg),nrow=nrow(Xreg),byrow=T)
  for (i.variable in 1:ncol(Xreg)){
    if (is.na(apply(Xreg_st,2,sum)[i.variable])==T) {
      Xreg_st[,i.variable]<-0}
  }

  ## Sparse LTS Fit
  SPARSELTS<-robustHD::sparseLTS(x=Xreg_st,y=c(Yreg),lambda=lambdaseq,mode = mode,intercept=F,normalize=T,initial="sparse", tol =tol,ncores=ncores,alpha=alpha.sparselts)
  SPARSELTS$coefficients<-as.matrix(SPARSELTS$coefficients)
  NONZERO_SOLUTION<-apply(SPARSELTS$coefficients,2,NONZERO)
  COEF_REDUCED<-matrix(SPARSELTS$coefficients[,which(NONZERO_SOLUTION!=0)],nrow=nrow(SPARSELTS$coefficients))
  LAMBDA_REDUCED<-lambdaseq[which(NONZERO_SOLUTION!=0)]

  if (ncol(COEF_REDUCED)>1){
    COEF_FINAL<-matrix(COEF_REDUCED[,which.min(SPARSELTS$crit$values[which(NONZERO_SOLUTION!=0),1])],ncol=1)
    COEF_FINAL[which(apply(Xreg,2,mad)!=0),]<-COEF_FINAL[which(apply(Xreg,2,mad)!=0),]/apply(Xreg,2,mad)[which(apply(Xreg,2,mad)!=0)]
    COEF_FINAL<-apply(COEF_FINAL,2,NORMALIZATION_UNIT)
    LAMBDA_FINAL<-LAMBDA_REDUCED[which.min(SPARSELTS$crit$values[which(NONZERO_SOLUTION!=0),1])]
  } else {
    COEF_FINAL<-COEF_REDUCED
    COEF_FINAL[which(apply(Xreg,2,mad)!=0),]<-COEF_FINAL[which(apply(Xreg,2,mad)!=0),]/apply(Xreg,2,mad)[which(apply(Xreg,2,mad)!=0)]
    COEF_FINAL<-apply(COEF_FINAL,2,NORMALIZATION_UNIT)
    LAMBDA_FINAL<-LAMBDA_REDUCED
  }


  ## OUTPUT
  out<-list(COEF_FINAL=COEF_FINAL,LAMBDA_FINAL=LAMBDA_FINAL)
}

NORMALIZATION_UNIT<-function(U){
  # AUX FUNCTION
  # output: normalized vector U
  length.U<-as.numeric(sqrt(t(U)%*%U))
  if(length.U==0){length.U<-1}
  Uunit<-U/length.U
}

BIC<-function(U,Y.data,X.data){
  # AUX FUNCTION
  ### Calculate value of Bayesian information Criterion
  sigmahat<-sum(diag((1/nrow(Y.data))*(t(Y.data-X.data*U[1,])%*%(Y.data-X.data*U[1,]))))
  BIC<-nrow(Y.data)*log(sigmahat) + length(which(U!=0))*log(nrow(Y.data))
  return(BIC)
}

NONZERO<-function(U){
  # AUX FUNCTION
  # output: number of nonzero components
  return(length(which(U!=0)))
}
