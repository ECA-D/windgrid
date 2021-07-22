library(geosphere)

# kriging lag
krigingLag <- function(x1,x2)
{
  N1 = dim(x1)
  N2 = dim(x2)
  X1 = meshgrid(x2[,1,drop=FALSE],x1[,1,drop=FALSE])
  X2 = meshgrid(x2[,2,drop=FALSE],x1[,2,drop=FALSE])
  n1 = size(X1$Y)[1]
  n2 = size(X1$Y)[2]
  lon1 = matrix(X1$Y,n1*n2,1)
  lon2 = matrix(X1$X,n1*n2,1)
  lat1 = matrix(X2$Y,n1*n2,1)
  lat2 = matrix(X2$X,n1*n2,1)
  H = abs(distCosine(cbind(lon1,lat1),cbind(lon2,lat2)))
  #H = distVincentyEllipsoid(cbind(lon1,lat1),cbind(lon2,lat2))
  H = matrix(H,n1,n2)
  return(H)
}

ensembleRankHistStd <- function(ensemble,reference,printing)
{
  N = size(ensemble)[1]
  nMember = size(ensemble)[2]
  #nReduced = max(8,ceil(N/50))
  nReduced = 8
  count = matrix(0,1,nReduced)
  for(n in 1:N){
    ensfull = ensemble[n,]
    ens = interp1((0.5:(nMember-0.5))/nMember,ensfull,(0.5:(nReduced-0.5))/nReduced)
    ref = reference[n]
    if(ref<min(ens)){ref=min(ens)}
    if(ref>max(ens)){ref=max(ens)}
    id = interp1(ens,1:nReduced,xi=ref,method="nearest")
    count[id] = count[id]+1
  }
  if(printing){print(count)}
  freq = count/N
  return(std(freq))
}

ensembleDiversityFactoring <- function(ensemble,factor)
{
  nMember = size(ensemble)[2]
  ensembleMean = NA*ensemble
  mu = rowMeans(ensemble)
  for(k in 1:nMember){
    ensembleMean[,k] = mu
  }
  newEnsemble = ensembleMean + factor*(ensemble-ensembleMean)
  return(newEnsemble)
}

ensembleDiversityFactoringNdim <- function(ensemble,dist2nearestAnalyse,dist2nearestPredict,factor)
{
  # determine factor based on distance to nearest neighbour
  factorInterp = factorInt(dist2nearestAnalyse,dist2nearestPredict,factor)
  
  # apply factoring
  nMember = size(ensemble)[2]
  newEnsemble = NA*ensemble
  ensembleMean = rowMeans(ensemble)
  for(k in 1:nMember){
    newEnsemble[,k] = ensembleMean + factorInterp*(ensemble[,k]-ensembleMean)
  }
  return(newEnsemble)
}

ensembleDiversityFactoringNdim2 <- function(ensemble,dist2nearest,binDist,factor)
{
  # determine factor based on distance to nearest neighbour
  ndim = length(factor)
  distMin = min(min(binDist),min(dist2nearest))
  distMax = max(max(binDist),max(dist2nearest))
  binDist = c(distMin,binDist,distMax)
  factor = c(factor[1],factor,factor[ndim])
  factorInterp = interp1(log(binDist),log(factor),log(dist2nearest),method='linear')
  
  # apply factoring
  nMember = size(ensemble)[2]
  newEnsemble = NA*ensemble
  ensembleMean = rowMeans(ensemble)
  for(k in 1:nMember){
    newEnsemble[,k] = ensembleMean + factorInterp*(ensemble[,k]-ensembleMean)
  }
  return(newEnsemble)
}

factorInt <- function(dist2nearestAnalyse,dist2nearestPredict,factor){
  ndim = size(factor)[2]
  if(ndim==1){
    factor = cbind(factor,factor)
    ndim = 2
  }
  distP = linspace(0,1,ndim)
  distQ = unname(quantile(dist2nearestAnalyse,probs=distP))
  factor = c(factor[1],factor[1:ndim],factor[ndim])
  minDist = min(distQ[1],min(dist2nearestPredict))
  maxDist = max(distQ[ndim],max(dist2nearestPredict))
  distQ = c(minDist,distQ[1:ndim],maxDist)
  factorInterp = exp(interp1(log(distQ),log(factor),log(dist2nearestPredict),method='linear'))
  return(factorInterp)
}
