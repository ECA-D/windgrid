#dev.off()

library(pracma)
require(Hmisc)
library(readr)
library(matrixStats)
library(tictoc)
library(ensembleBMA)
library(stats)
library(utils)
library(Matrix)
library(spam)
library(qlcMatrix)
library(sp)
library(yaml)
suppressWarnings(library(randomFunctions)) # Prevent trivial warnings

#setwd('~/Jouke/WindGrid/scripts/FullUpdates')

##########################################
#                                        #
#        Get gridding paramters          #
#                                        #
##########################################

args <- getArgs()
iyaml <- yaml.load_file(args$yaml)

########################### Set Arguments ############################
gen.args <- c("stn_db","lonmin","lonmax","latmin","latmax","resolution",
	      "nEnsemble","verbose","txtload_source","NstationMax",
	      "degTreshold","ivar")

gridding.args <- c("gpr_source","edit_source","relnoisey","lagTreshold","thetaMin",
                   "thetaMax","nThetaBrute","plots","useBackground0",
                   "corrPower","randomSeed","nfold","allvars","useSparseMatrices",
                   "xvalSeed","EDIToptimMethod","EDIToptimMaxIter","writeRawEnsemble",
                   "writeTunedEnsemble")

for(x in gen.args){ assign(x,iyaml[[x]]) }
for(x in gridding.args){ assign(x,iyaml$gridding[[x]]) }

source(txtload_source)
source(gpr_source)
source(edit_source)



## Define outputDir directory
if(is.null(args$scratch)){
    outputDir <- iyaml$scratch
} else {
    outputDir <- args$scratch
}

## Define inputDir directory
## This script uses temporary output from preprocess and produces
## other temporary output, so inputDir is the same as outputDir 
if(is.null(args$scratch)){
    inputDir <- iyaml$scratch
} else {
    inputDir <- args$scratch
}

if(verbose>2){tic()}
if(verbose>2){
  print("useSparseMatrices:")
  print(useSparseMatrices)
  print("lagTreshold:")
  print(lagTreshold)
  print("EDIToptimMaxIter:")
  print(EDIToptimMaxIter)
  print("inputDir:")
  print(inputDir)
  print("outputDir:")
  print(outputDir)
}

# Define years
if(is.null(args$year)){
   YEAR <- 1950
} else {
    yrsplit <- strsplit(args$year,";")[[1]]
    YEAR <- as.numeric(yrsplit)
    YEAR <- if(length(YEAR)==1) YEAR else YEAR[1]:YEAR[2]
}

# Define months
if(is.null(args$month)){
   MONTH <- 1:12
} else {
    mntsplit <- strsplit(args$month,";")[[1]]
    MONTH <- as.numeric(mntsplit)
    MONTH <- if(length(MONTH)==1) MONTH else MONTH[1]:MONTH[2]
}

nlineblockwise = 100

######################
### LOOP OVER YEAR ###
######################
for(year in YEAR){
  
  ystr = toString(year)
  
  backgroundgrid = readRDS(file = paste(inputDir,"backgroundgrid_",ystr,".rds",sep=""))

  ### load DEM ###
  demGrid = readRDS(file = paste(outputDir,"grid_",toString(year),".rds",sep=""))
  TPI = matrix(0,length(demGrid$TPI),1)
  TPI[is.finite(demGrid$TPI)] = demGrid$TPI[is.finite(demGrid$TPI)]
  logTPI = log(abs(TPI)+0.01)
  if(plots){
    Nlon = length(backgroundgrid$lon)
    Nlat = length(backgroundgrid$lat)
    png(filename=paste('figures/grid_logtpi_',toString(year),'.png',sep=''))
        filled.contour(backgroundgrid$lon,backgroundgrid$lat,t(matrix(logTPI,Nlat,Nlon)),color.palette=colorRampPalette(c("blue","yellow","red")))
        dev.off()
  }
  
  #######################
  ### LOOP OVER MONTH ###
  #######################
  for(month in MONTH){
    
    if(month<10){
      mstr = paste("0",toString(month),sep="")
    } else {
      mstr = toString(month)
    }
    
    useBackground = useBackground0
    test = max(is.finite(backgroundgrid$LOGBACKGROUND_MEMBER[,month,1]))
    if(test<1){
      useBackground = FALSE
    }
    
    if(verbose>0){print(paste('== year',ystr,'-> month',mstr,'=='))}
    nDays = monthDays(as.Date(paste(ystr,'-',mstr,'-01',sep='')))
    
    ### load station data ###
    constantStations = TRUE
    if(verbose>1){print('Reading station data for current month ...')}
    stationData = readStationData(year,month,stn_db,allvars,lonmin,lonmax,latmin,latmax,ivar)
    
    NNstation = size(stationData$speed)[1]
    
    ### IF constant stations THEN compute constant monthly analysis and prediction lags ###
    if(constantStations){
      xStation = cbind(stationData$lon[1:NNstation,1,drop=FALSE],stationData$lat[1:NNstation,1,drop=FALSE])
      xGrid = cbind(backgroundgrid$LON,backgroundgrid$LAT)
      
      if(verbose>1){print('Computing constant monthly analysis lags')}
      if(useSparseMatrices){
        monthLagAnalyseSparse = krigingLagSparseFaster(xStation,xStation,lagTreshold,nlineblockwise)
        if(verbose>2){print(object.size(monthLagAnalyseSparse))}
      } else {
      	monthLagAnalyse = krigingLag(xStation,xStation)
        if(verbose>2){print(object.size(monthLagAnalyse))}
      }

      Ngrid = size(xGrid)[1]
      nblock = ceiling(Ngrid/nlineblockwise)
      idblock = seq(0,Ngrid,by=nlineblockwise)
      if(length(idblock)<=nblock){
        idblock = cbind(idblock,Ngrid)
      }
      if(verbose>1){print('Computing constant monthly prediction mask')}
      MaskInGrid = matrix(NA,Ngrid,1)
      for(k in 1:nblock){
        idnow = (idblock[k]+1):idblock[k+1]
        lag = krigingLag(xGrid[idnow,,drop=FALSE],xStation)
        MaskInGrid[idnow] = (rowMins(lag)<=lagTreshold)
      }
      MaskInGrid4plot = round(MaskInGrid)
      nInMask = sum(MaskInGrid4plot)
      
      if(verbose>1){print('Computing constant monthly prediction lags')}
      if(useSparseMatrices){
	monthLagPredictSparse = krigingLagSparseFaster(xGrid[MaskInGrid,,drop=FALSE],xStation,lagTreshold,nlineblockwise)
        if(verbose>2){print(object.size(monthLagPredictSparse))}
      } else {
      	monthLagPredict = krigingLag(xGrid[MaskInGrid,,drop=FALSE],xStation) # can run  out of memory on local machine
        if(verbose>2){print(object.size(monthLagPredict))}
      }
      
    } else {
      
      stop('this part of the code should have become redundant')
      
    }
    
    ### get grid storage sizes ###
    nMember = size(backgroundgrid$LOGBACKGROUND_MEMBER)[3]
    nDays = monthDays(as.Date(paste(ystr,'-',mstr,'-01',sep='')))
    
    ### pre-allocate daily grids and performance summaries ###
    grid = list(lon=backgroundgrid$lon,lat=backgroundgrid$lat)
    grid$LON = backgroundgrid$LON
    grid$LAT = backgroundgrid$LAT
    grid$WINDSPEED_ENSEMBLE_RAW = array(NA,c(Ngrid,nMember,nDays))
    grid$WINDSPEED_ENSEMBLE_TUNED = array(NA,c(Ngrid,nMember,nDays))
    grid$WINDSPEED_ENSEMBLE_MEAN = matrix(NA,Ngrid,nDays)
    grid$WINDSPEED_ENSEMBLE_RAW_SPREAD_QUANT = matrix(NA,Ngrid,nDays)
    grid$WINDSPEED_ENSEMBLE_TUNED_SPREAD_QUANT = matrix(NA,Ngrid,nDays)
    grid$WINDSPEED_ENSEMBLE_RAW_SPREAD_STD = matrix(NA,Ngrid,nDays)
    grid$WINDSPEED_ENSEMBLE_TUNED_SPREAD_STD = matrix(NA,Ngrid,nDays)
    grid$rmse_residual = matrix(NA,1,nDays)
    grid$rmse_xval = matrix(NA,1,nDays)
    grid$rmse_spread = matrix(NA,1,nDays)
    grid$tuning_std_raw = matrix(NA,1,nDays)
    grid$tuning_std_tuned = matrix(NA,1,nDays)
    grid$tuning_curveRMSE_raw = matrix(NA,1,nDays)
    grid$tuning_curveRMSE_tuned = matrix(NA,1,nDays)
    grid$tuning_coeff = matrix(NA,3,nDays)

    #######################################
    ### loop over days to compute grids ###
    #######################################
    for(day in 1:nDays){
      
      if(day<10){
        dstr = paste("0",toString(day),sep="")
      } else {
        dstr = toString(day)
      }
      if(verbose>2){print(paste(ystr,mstr,dstr))}
      
      ### Compute anomalies ###
      stationBackgroLog = matrix(NA,NNstation,nMember)
      Nlon = length(backgroundgrid$lon)
      Nlat = length(backgroundgrid$lat)
      for(member in 1:nMember){
        if(useBackground){
          stationBackgroLog[,member] = interp2(backgroundgrid$lat,backgroundgrid$lon,t(matrix(backgroundgrid$LOGBACKGROUND_MEMBER[,month,member,drop=FALSE],Nlat,Nlon)),stationData$lat[,day,drop=FALSE],stationData$lon[,day,drop=FALSE])
        } else {
          if(verbose>2){print('>  Computing todays anomalies (not using background)')}
          stop('This part of the code should have become redundant')
        }
      }
      stationAnomalyLog = log(stationData$speed[,day]) - stationBackgroLog
      
      ### select today's finite data ###
      idFinSpeed = is.finite(rowSums(stationAnomalyLog))
      xStation = cbind(stationData$lon[idFinSpeed,day,drop=FALSE],stationData$lat[idFinSpeed,day,drop=FALSE])
      yStation = stationAnomalyLog[idFinSpeed,,drop=FALSE]
      n = size(yStation)[1]
      
      ### kriging: normalise ###
      mu = colMeans(yStation)
      st = sqrt(colVars(yStation))
      for(member in 1:nMember){
        yStation[,member] = (yStation[,member] - mu[member]) / st[member]
      }
      
      if(constantStations){
        
        ### maximum likelihood estimate ###
        if(verbose>2){print('>  Todays maximum likelihood estimate')}
        thetaBrute = exp(linspace(log(thetaMin),log(thetaMax),n=nThetaBrute))
        negLogLike = NA*thetaBrute
        yStationMean = rowMeans(yStation)
        for(iterBrute in 1:nThetaBrute){
          if(useSparseMatrices){
            cutoff = exp(-.5*(lagTreshold^corrPower)/(thetaBrute[iterBrute]^corrPower))
            P = monthLagAnalyseSparse
            P@x = exp(-.5*((lagTreshold-monthLagAnalyseSparse@x)^corrPower)/(thetaBrute[iterBrute]^corrPower))
            P@x = (P@x-cutoff) / (1.-cutoff)
            P = P[idFinSpeed,idFinSpeed,drop=FALSE]
            R = sparseMatrix(i=1:sum(idFinSpeed),j=1:sum(idFinSpeed),x=relnoisey^2)
            A = R + P
            eigs = eigen(A,only.values=TRUE)
            negLogLike[iterBrute] = sum(log(eigs$values)) + t(yStationMean)%*%qr.solve(A,yStationMean)
          } else {
            P = exp(-.5*((monthLagAnalyse[idFinSpeed,idFinSpeed,drop=FALSE])^corrPower)/(thetaBrute[iterBrute]^corrPower))
            R = relnoisey^2 * eye(size(P)[1])
            A = R + P
            negLogLike[iterBrute] = sum(log(eig(A))) + t(yStationMean)%*%qr.solve(A,yStationMean)
          }
        }
        if(plots){
          png(filename=paste("figures/MLE_",ystr,"_",mstr,"_",dstr,".png",sep=""))
          semilogx(thetaBrute,negLogLike,type='b',xlab='Theta',ylab='neg log like')
          dev.off()
        }
        idfin = is.finite(negLogLike)
        thetaBrute = thetaBrute[idfin]
        negLogLike = negLogLike[idfin]
        idmin = which.min(negLogLike)
        theta = thetaBrute[idmin]
        if(verbose>2){print(paste('   > thetaMLE =',toString(theta)))}       

        # kriging: analysis & prediction matrices
        if(verbose>1){print('Todays kriging analysis and prediction')}
        if(useSparseMatrices){
          cutoff = exp(-.5*(lagTreshold^corrPower)/(theta^corrPower))
          P = monthLagAnalyseSparse
          P@x = exp(-.5*((lagTreshold-monthLagAnalyseSparse@x)^corrPower)/(theta^corrPower))
          P@x = (P@x-cutoff) / (1.-cutoff)
          P = P[idFinSpeed,idFinSpeed,drop=FALSE]
          R = sparseMatrix(i=1:sum(idFinSpeed),j=1:sum(idFinSpeed),x=relnoisey^2)
          A = R + P
          b = monthLagPredictSparse
          b@x = exp(-.5*((lagTreshold-monthLagPredictSparse@x)^corrPower)/(theta^corrPower))
          b@x = (b@x-cutoff) / (1.-cutoff)
          b = b[,idFinSpeed,drop=FALSE]
          Psparse = P
          Rsparse = R
          Asparse = A
          bSparse = b
        } else {
          P = exp(-.5*((monthLagAnalyse[idFinSpeed,idFinSpeed,drop=FALSE])^corrPower)/(theta^corrPower))
          R = relnoisey^2 * eye(size(P)[1])
          A = R + P
          b = exp(-.5*((monthLagPredict[,idFinSpeed,drop=FALSE])^corrPower)/(theta^corrPower))
          Pfull = P
          Rfull = R
          Afull = A
          bFull = b
        }
        
        # kriging prediction
        yin = b %*% qr.solve(A,yStation)

        # kriging x-validation prediction
        N = size(yStation)[1]
        yinXval = matrix(NA,N,nMember)
        set.seed(xvalSeed)
        id0 = matrix(sample(1:N,N),N,1)
        Nout = N/nfold
        if(useSparseMatrices){
          dist2allSparse = monthLagAnalyseSparse[idFinSpeed,idFinSpeed,drop=FALSE]
          dist2allSparse[dist2allSparse==lagTreshold] = 0.
        } else {
          dist2all = monthLagAnalyse[idFinSpeed,idFinSpeed,drop=FALSE]
          dist2all[dist2all==0] = Inf
        }
        dist2nearestAnalyse = matrix(NA,N,1)

        for(k in 1:nfold){
          idout = id0[ceil((k-1)*Nout+1):ceil(k*Nout)]
          yinXval[idout,] = as.matrix(P[idout,-idout] %*% qr.solve(A[-idout,-idout],yStation[-idout,,drop=FALSE]))
          if(useSparseMatrices){
            dist2nearestAnalyse[idout,] = lagTreshold - rowMax(dist2allSparse[idout,-idout,drop=FALSE])
	  } else {
            dist2nearestAnalyse[idout,] = rowMins(dist2all[idout,-idout,drop=FALSE])        
          }
        end
        dist2nearestAnalyse[dist2nearestAnalyse>lagTreshold] = lagTreshold
    }
        
      } else {
        
        stop('this part of the code should have become redundant')
        
      }
      
      ### kriging: de-normalise ###
      gridAnomalyLog = NA*yin
      xvalAnomalyLog = NA*yinXval
      for(member in 1:nMember){
        gridAnomalyLog[,member] = mu[member] + st[member]*yin[,member,drop=FALSE]
        xvalAnomalyLog[,member] = mu[member] + st[member]*yinXval[,member,drop=FALSE]
      }
      
      if(verbose>2){print('Building and analysing ensembles')}
      
      ### build and analyse ensemble ###
      if(useBackground){
        gridToday = exp( gridAnomalyLog + backgroundgrid$LOGBACKGROUND_MEMBER[MaskInGrid,month,] )
        xvalToday = exp( xvalAnomalyLog + stationBackgroLog[idFinSpeed,,drop=FALSE] )
      } else {
        stop('This part of the code should have become redundant')
        gridToday = exp( gridAnomalyLog )
      }
      if(useSparseMatrices){
        gridToday = as.array(gridToday)
        xvalToday = as.array(xvalToday)
      }

      ### statistics: raw ensemble ###
      gridTodayEnsembleRaw = matrix(NA,Ngrid,nMember)
      gridTodayEnsembleRaw[MaskInGrid,] = gridToday
      gridTodayMean = matrix(NA,Ngrid,1)
      gridTodayMean[MaskInGrid] = rowMeans(gridToday)
      gridTodaySpreadRawQuant = matrix(NA,Ngrid,1)
      gridTodaySpreadRawStd = matrix(NA,Ngrid,1)
      gridTodaySpreadRawQuant[MaskInGrid] = rowQuantiles(gridToday,probs=.95) - rowQuantiles(gridToday,probs=.05)
      gridTodaySpreadRawStd[MaskInGrid] = sqrt(rowVars(gridToday))

      ### prediction errors ###
      stationPred = interp2(backgroundgrid$lat,backgroundgrid$lon,t(matrix(gridTodayMean,Nlat,Nlon)),stationData$lat[idFinSpeed,day,drop=FALSE],stationData$lon[idFinSpeed,day,drop=FALSE])
      Eresidual = stationPred - stationData$speed[idFinSpeed,day]
      RMSEresidual = sqrt(mean(Eresidual^2))
      if(verbose>2){print(paste('>  RMSE residual =',toString(RMSEresidual)))}
      
      Exval = rowMeans(xvalToday) - stationData$speed[idFinSpeed,day]
      RMSExval = sqrt(mean(Exval^2))
      if(verbose>2){print(paste('>  RMSE xval =',toString(RMSExval)))}

      statPred = interp2(backgroundgrid$lat,backgroundgrid$lon,t(matrix(gridTodaySpreadRawStd,Nlat,Nlon)),stationData$lat[idFinSpeed,day,drop=FALSE],stationData$lon[idFinSpeed,day,drop=FALSE])
      RMSEspread = sqrt(mean(statPred^2))
      if(verbose>2){print(paste('>  RMSE spread =',toString(RMSEspread)))}

      grid$rmse_residual[day] = RMSEresidual
      grid$rmse_xval[day] = RMSExval
      grid$rmse_spread[day] = RMSEspread

      ############
      ### EDIT ###
      ############

      # EDIT: raw ensemble
      if(verbose>1){print('Running EDIT ...')}
      if(useSparseMatrices){
        dist2allSparse = monthLagPredictSparse[,idFinSpeed,drop=FALSE]
      	dist2nearestPredict = lagTreshold - rowMax(dist2allSparse)
        dist2nearestPredict[dist2nearestPredict>lagTreshold] = lagTreshold
        dist2nearestPredictSparse = dist2nearestPredict
      } else {
      	dist2all = monthLagPredict[,idFinSpeed,drop=FALSE]
      	dist2nearestPredict = rowMins(dist2all)
        dist2nearestPredict[dist2nearestPredict>lagTreshold] = lagTreshold
        dist2nearestPredictFull = dist2nearestPredict
      }
      stationLogTPI = interp2(backgroundgrid$lat,backgroundgrid$lon,t(matrix(logTPI,Nlat,Nlon)),stationData$lat[idFinSpeed,day,drop=FALSE],stationData$lon[idFinSpeed,day,drop=FALSE])
      gridLogTPI = logTPI[MaskInGrid,,drop=FALSE]
      EDITproxyStation = cbind(log(dist2nearestAnalyse),stationLogTPI)
      EDITproxyGrid = cbind(log(dist2nearestPredict),gridLogTPI)
      quantileStd = ensembleRankHistStd(xvalToday,stationData$speed[idFinSpeed,day,drop=FALSE],FALSE)
      preEDITrankHist = ensembleRankHist(xvalToday,stationData$speed[idFinSpeed,day,drop=FALSE],FALSE)
      if(verbose>2){print(paste('>  EDIT: Raw Xval rank hist RStd =',toString(quantileStd)))}
      curveRMSEraw = curveRMSEbyProxy(EDITproxyStation,Exval,xvalToday)
      if(verbose>2){print(paste('>  EDIT: Raw Xval proxy binning RMSRE =',toString(curveRMSEraw)))}

      # EDIT: optimizer
      optimMaxIter = EDIToptimMaxIter
      reference = stationData$speed[idFinSpeed,day,drop=FALSE]
      par0 = c(0,0*EDITproxyStation[1,,drop=FALSE])
      opt1 = optim(par=par0,fn=EDITcostfunction,EDITproxyStation=EDITproxyStation,ensemble=xvalToday,reference=reference,Exval=Exval,minvals=c(NA,NA),type=1,negativeAllowed=TRUE,control = list(maxit=optimMaxIter),method=EDIToptimMethod)
      opt2 = optim(par=par0,fn=EDITcostfunction,EDITproxyStation=EDITproxyStation,ensemble=xvalToday,reference=reference,Exval=Exval,minvals=c(NA,NA),type=2,negativeAllowed=TRUE,control = list(maxit=optimMaxIter),method=EDIToptimMethod)
      opt3 = optim(par=par0,fn=EDITcostfunction,EDITproxyStation=EDITproxyStation,ensemble=xvalToday,reference=reference,Exval=Exval,minvals=c(opt1$value,opt2$value),type=3,negativeAllowed=TRUE,control = list(maxit=optimMaxIter),method=EDIToptimMethod)
      coeffMinOptimizer = t(opt3$par)
      ensembleTestNegAllowed = EDITfactorByProxy(EDITproxyStation,EDITproxyStation,xvalToday,coeffMinOptimizer,TRUE)
      stdMinOptimizerNegAllowed = ensembleRankHistStd(ensembleTestNegAllowed,reference,FALSE)
      curveRMSEminOptimizerNegAllowed = curveRMSEbyProxy(EDITproxyStation,Exval,ensembleTestNegAllowed)
      ensembleTest = EDITfactorByProxy(EDITproxyStation,EDITproxyStation,xvalToday,coeffMinOptimizer,FALSE)
      stdMinOptimizer = ensembleRankHistStd(ensembleTest,reference,FALSE)
      curveRMSEminOptimizer = curveRMSEbyProxy(EDITproxyStation,Exval,ensembleTest)
      if(verbose>2){print(paste('>  EDIT: Optimizer tuned coefficients =',toString(round(100*coeffMinOptimizer)/100)))}
      if(verbose>2){print(paste('>  EDIT: Optimizer tuned rank hist RStd neg all =',toString(stdMinOptimizerNegAllowed)))}
      if(verbose>2){print(paste('>  EDIT: Optimizer tuned proxy binning RMSRE neg all =',toString(curveRMSEminOptimizerNegAllowed)))}
      if(verbose>2){print(paste('>  EDIT: Optimizer tuned rank hist RStd =',toString(stdMinOptimizer)))}
      if(verbose>2){print(paste('>  EDIT: Optimizer tuned proxy binning RMSRE =',toString(curveRMSEminOptimizer)))}

      grid$tuning_coeff[,day] = t(coeffMinOptimizer)
      grid$tuning_std_raw[day] = quantileStd
      grid$tuning_std_tuned[day] = stdMinOptimizer
      grid$tuning_curveRMSE_raw[day] = curveRMSEraw
      grid$tuning_curveRMSE_tuned[day] = curveRMSEminOptimizer

      # EDIT: tune ensemble based on Pareto optimum
      xvalTodayTuned = EDITfactorByProxy(EDITproxyStation,EDITproxyStation,xvalToday,coeffMinOptimizer,FALSE)
      gridTodayTuned = EDITfactorByProxy(EDITproxyStation,EDITproxyGrid,gridToday,coeffMinOptimizer,FALSE)
      ensMin = min(gridToday[is.finite(gridToday)])
      if(verbose>2){print(ensMin)}
      if(ensMin<0.0){
        if(verbose>2){print('>  WARNING: negative raw ensemble member encountered')}
      }
      ensMinTuned = min(gridTodayTuned[is.finite(gridTodayTuned)])
      if(verbose>2){print(ensMinTuned)}
      if(ensMinTuned<0.0){
        if(verbose>2){print('>  WARNING: negative tuned ensemble member encountered')}
      }

      postEDITrankHist = ensembleRankHist(xvalTodayTuned,stationData$speed[idFinSpeed,day,drop=FALSE],FALSE)

      # EDIT: populate grid
      gridTodayEnsembleTuned = matrix(NA,Ngrid,nMember)
      gridTodayEnsembleTuned[MaskInGrid,] = gridTodayTuned
      gridTodaySpreadTunedQuant = matrix(NA,Ngrid,1)
      gridTodaySpreadTunedStd = matrix(NA,Ngrid,1)
      gridTodaySpreadTunedQuant[MaskInGrid] = rowQuantiles(gridTodayTuned,probs=.95) - rowQuantiles(gridTodayTuned,probs=.05)
      gridTodaySpreadTunedStd[MaskInGrid] = sqrt(rowVars(gridTodayTuned))

      ################
      ### Plotting ###
      ################
      if(plots){

        Nlon = length(backgroundgrid$lon)
        Nlat = length(backgroundgrid$lat)
        png(filename=paste('figures/final_grid_mean_',toString(year),'_',mstr,'_',dstr,'.png',sep=''))
        filled.contour(backgroundgrid$lon,backgroundgrid$lat,t(matrix(gridTodayMean,Nlat,Nlon)),color.palette=colorRampPalette(c("blue","yellow","red")))
        dev.off()
        png(filename=paste('figures/final_grid_log10mean_',toString(year),'_',mstr,'_',dstr,'.png',sep=''))
        filled.contour(backgroundgrid$lon,backgroundgrid$lat,t(matrix(log10(gridTodayMean),Nlat,Nlon)),color.palette=colorRampPalette(c("blue","yellow","red")))
        dev.off()
        png(filename=paste('figures/final_grid_log10spreadRawStd_',toString(year),'_',mstr,'_',dstr,'.png',sep=''))
        filled.contour(backgroundgrid$lon,backgroundgrid$lat,t(matrix(log10(gridTodaySpreadRawStd),Nlat,Nlon)),color.palette=colorRampPalette(c("blue","yellow","red")))
        dev.off()
        png(filename=paste('figures/final_grid_log10spreadEDITedStd_',toString(year),'_',mstr,'_',dstr,'.png',sep=''))
        filled.contour(backgroundgrid$lon,backgroundgrid$lat,t(matrix(log10(gridTodaySpreadTunedStd),Nlat,Nlon)),color.palette=colorRampPalette(c("blue","yellow","red")))
        dev.off()
        idFinSpeedRaw = is.finite(stationData$speedRaw[,day])
        png(filename=paste('figures/final_station2station_',toString(year),'_',mstr,'_',dstr,'.png',sep=''))
        #statPred = interp2(backgroundgrid$lat,backgroundgrid$lon,t(matrix(gridTodayMean,Nlat,Nlon)),stationData$latRaw[idFinSpeedRaw,day,drop=FALSE],stationData$lonRaw[idFinSpeedRaw,day,drop=FALSE])
        #loglog(stationData$speedRaw[idFinSpeedRaw,day],statPred,xlab='raw station observation',ylab='grid @ station')
        #lines(stationData$speedRaw[idFinSpeedRaw,day],stationData$speedRaw[idFinSpeedRaw,day],col='red')
        #dev.off()
        #png(filename=paste('figures/final_station2stationSpreadStd_',toString(year),'_',mstr,'_',dstr,'.png',sep=''))
        #statPred = interp2(backgroundgrid$lat,backgroundgrid$lon,t(matrix(gridTodaySpreadRawStd,Nlat,Nlon)),stationData$latRaw[idFinSpeedRaw,day,drop=FALSE],stationData$lonRaw[idFinSpeedRaw,day,drop=FALSE])
        #loglog(stationData$speedRaw[idFinSpeedRaw,day],statPred,xlab='raw station observation',ylab='grid spread @ station')
        #lines(stationData$speedRaw[idFinSpeedRaw,day],stationData$speedRaw[idFinSpeedRaw,day],col='red')
        #dev.off()
       
        # EDIT plots
        png(filename=paste('figures/ensemble_tuning_rank_histogram_',toString(year),'_',mstr,'_',dstr,'.png',sep=''))
        plot(1:21,preEDITrankHist,col='blue',xlab='Rank',ylab='Frequency',type='b')
        points(1:21,postEDITrankHist,col='red',type='b')
        lines(1:21,0*(1:21)+(1/21),col='black')
        dev.off()

        if(day==nDays){
          png(filename=paste('figures/ensemble_tuning_hist_relStd_',toString(year),'_',mstr,'.png',sep=''))
          hist(grid$tuning_std_tuned/grid$tuning_std_raw,xlab='tuned/raw quantile binning rstd',ylab='count')
          #title('Reduction factor of quantile binning std')
          dev.off()
          png(filename=paste('figures/ensemble_tuning_hist_relRMSE_',toString(year),'_',mstr,'.png',sep=''))
          hist(grid$tuning_curveRMSE_tuned/grid$tuning_curveRMSE_raw,xlab='tuned/raw proxy binning RMSRE',ylab='count')
          #title('Reduction factor of proxy binning RMSRE')
          dev.off()
          png(filename=paste('figures/ensemble_tuning_hist_relVSrel_',toString(year),'_',mstr,'.png',sep=''))
          plot(grid$tuning_curveRMSE_tuned/grid$tuning_curveRMSE_raw,grid$tuning_std_tuned/grid$tuning_std_raw,ylab='tuned/raw quantile binning rstd',xlab='tuned/raw proxy binning RMSRE')
          #title('Reduction factor of proxy binning RMSRE')
          dev.off()
        }
      }

      ### this is where today's values are stored: ###
      grid$WINDSPEED_ENSEMBLE_RAW[,,day] = gridTodayEnsembleRaw
      grid$WINDSPEED_ENSEMBLE_TUNED[,,day] = gridTodayEnsembleTuned
      grid$WINDSPEED_ENSEMBLE_MEAN[,day] = gridTodayMean
      grid$WINDSPEED_ENSEMBLE_RAW_SPREAD_QUANT[,day] = gridTodaySpreadRawQuant
      grid$WINDSPEED_ENSEMBLE_TUNED_SPREAD_QUANT[,day] = gridTodaySpreadTunedQuant
      grid$WINDSPEED_ENSEMBLE_RAW_SPREAD_STD[,day] = gridTodaySpreadRawStd
      grid$WINDSPEED_ENSEMBLE_TUNED_SPREAD_STD[,day] = gridTodaySpreadTunedStd  
    }
    
    if(month<10){
      mstr = paste("0",toString(month),sep="")
    } else {
      mstr = toString(month)
    }

    #####################
    ### write to disk ###
    #####################
    if(!writeRawEnsemble){
      grid$WINDSPEED_ENSEMBLE_RAW = NULL
    }
    if(!writeTunedEnsemble){
      grid$WINDSPEED_ENSEMBLE_TUNED = NULL
    }

    if(verbose>1){print(paste('Writing grid file for year',toString(year),'month',mstr))}
    saveRDS(grid, file = paste(outputDir,"windgrid_ensemble_statistics_",toString(year),"_",mstr,".rds",sep=""))
  }
}

if(verbose>2){toc()}
