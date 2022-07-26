verifRankHistNoPlot <- function(forecasts, observations) 
{
#
# copyright 2006-present, University of Washington. All rights reserved.
# for terms of use, see the LICENSE file
#
#rank <- apply((forecasts > observations),1,sum)+1
 rank <- apply(cbind(observations,forecasts),1,function(x)
               rank(x,ties="random")[1] )
 k <- ncol(forecasts)
 hst=hist(rank, breaks = 0:(k+1), prob = TRUE, xaxt = "n",
     xlab = "", ylab = "", main = "Verification Rank Histogram",plot=FALSE)
 #axis(1, at = seq(.5, to = k+.5, by = 1), labels = 1:(k+1))
 #abline(h=1/(k+1), lty=2)
 #invisible(rank)
 return(hst$counts)
}

ensembleRankHist <- function(ensemble,reference,printing)
{
  N = size(ensemble)[1]
  nMember = size(ensemble)[2]
  #count = matrix(0,1,nMember)
  #for(n in 1:N){
  #  ens = sort(ensemble[n,])
  #  ref = reference[n]
  #  if(ref<min(ens)){ref=min(ens)}
  #  if(ref>max(ens)){ref=max(ens)}
  #  id = interp1(ens,1:nMember,xi=ref,method="nearest")
  #  count[id] = count[id]+1
  #}
  #if(printing){print(count)}
  #count = Rankhist(ensemble,reference,reduce.bins=1,handle.na = "na.fail")
  count = verifRankHistNoPlot(ensemble,reference)
  #print(count)
  freq = count/N
  return(freq)
}

ensembleRankHistStd <- function(ensemble,reference,printing)
{
  freq = ensembleRankHist(ensemble,reference,printing)
  return(std(freq)/mean(freq))
}

curveRMSEbyProxy <- function(proxy,Exval,ensemble){
  nbin = 4
  proxy1 = linspace(min(proxy[,1]),max(proxy[,1]),nbin)
  proxy2 = linspace(min(proxy[,2]),max(proxy[,2]),nbin)
  binnedRMSExval = matrix(NA,nbin*nbin,1)
  binnedSTDensemble = matrix(NA,nbin*nbin,1)
  binNR = matrix(1:(nbin*nbin),nbin*nbin,1)
  binNRi = interp2(proxy1,proxy2,matrix(binNR,nbin,nbin),proxy[,1,drop=FALSE],proxy[,2,drop=FALSE],method='nearest')
  for(k in 1:(nbin*nbin)){
    idin = (binNRi == k)
    if(sum(idin)>20){
      binnedRMSExval[k] = sqrt(mean(Exval[idin]^2))
      binnedSTDensemble[k] = sqrt(mean(rowVars(ensemble[idin,,drop=FALSE])))
    }
  }
  id = is.finite(binnedRMSExval)
  curveRMSE = sqrt(mean(((binnedSTDensemble[id]-binnedRMSExval[id])/binnedRMSExval[id])^2))
  return(curveRMSE)
}

EDITfactorByProxy <- function(EDITproxyStation,EDITproxy,ensemble,coefficient,negativeAllowed){
  # normalize
  nProxy = size(EDITproxyStation)[2]
  mu = colMeans(EDITproxyStation)
  st = sqrt(colVars(EDITproxyStation))
  for(dim in 1:nProxy){
    EDITproxy[,dim] = (EDITproxy[,dim]-mu[dim])/st[dim]
  }
  
  # factor
  X = cbind(1,EDITproxy)
  #print(size(X))
  #print(coefficient)
  factorBase = exp(X%*%t(coefficient))
  factorMax = rowMeans(ensemble) / (rowMeans(ensemble)-rowMins(ensemble)+1e-6)
  if(negativeAllowed){
    factor = factorBase
  } else {
    factor = pmin(factorBase,factorMax)
  }

  # apply factoring
  nMember = size(ensemble)[2]
  newEnsemble = NA*ensemble
  ensembleMean = rowMeans(ensemble)
  for(k in 1:nMember){
    newEnsemble[,k] = ensembleMean + factor*(ensemble[,k]-ensembleMean)
  }
  return(newEnsemble)
}

EDITcostfunction <- function(pars,EDITproxyStation,ensemble,reference,Exval,minvals,type,negativeAllowed)
{
  coefficients = t(pars)
  #print(coefficients)

  # impose bounds
  #if(max(abs(coefficients))>1.0){return(1e6)}

  # factor ensemble
  ensembleTest = EDITfactorByProxy(EDITproxyStation,EDITproxyStation,ensemble,coefficients,negativeAllowed)

  # compute cost
  if(type==1)
  {
    stdTest = ensembleRankHistStd(ensembleTest,reference,FALSE)
    return(stdTest)
  }

  if(type==2)
  {
    curveRMSEtest = curveRMSEbyProxy(EDITproxyStation,Exval,ensembleTest)
    return(curveRMSEtest)
  }

  if(type==3)
  {
    stdTest = ensembleRankHistStd(ensembleTest,reference,FALSE)
    curveRMSEtest = curveRMSEbyProxy(EDITproxyStation,Exval,ensembleTest)
    return(sqrt( (stdTest/minvals[1])^2 + (curveRMSEtest/minvals[2])^2 ))
   }

}
