# linear basis
linearBasis <- function(X0,p){
  N = size(X0)[1]
  n1 = size(p)[1]
  n2 = size(p)[2]
  X = matrix(1,N,n1)
  for(j in 1:n1){
    for(k in 1:n2){
      X[,j] = X[,j,drop=FALSE]*X0[,k,drop=FALSE]^p[j,k]
    }
  }
  return(X)
}

# linear regression
linearRegression <- function(X0,p,y,X0i){
  X = linearBasis(X0,p)
  pLin = qr.solve(X,y)
  Xi = linearBasis(X0i,p)
  yi = Xi%*%pLin
}

# regression RMSE
lrRRMSE <- function(X0,p,y,nfold,seed,runXval){
  # training error
  yi = linearRegression(X0,p,y,X0)
  trainingRRMSE = sqrt(mean((yi-y)^2))/std(y)
  # cross-validation error
  if(runXval){
    N = size(y)[1]
    #s = get.seed()
    set.seed(seed)
    id0 = matrix(sample(1:N,N),N,1)
    #set.seed(s)
    Nout = N/nfold
    xvalMSE = 0
    for(k in 1:nfold){
      idout = id0[round((k-1)*Nout+1):round(k*Nout)]
      yi = linearRegression(X0[-idout,,drop=FALSE],p,y[-idout,,drop=FALSE],X0[idout,,drop=FALSE])
      xvalMSE = xvalMSE + mean((yi-y[idout,,drop=FALSE])^2)
    }
    xvalRRMSE = sqrt(xvalMSE/nfold)/std(y)
  } else {
    xvalRRMSE = NA
  }
  rrmse = list(training=trainingRRMSE,xval=xvalRRMSE)
  return(rrmse)
}
