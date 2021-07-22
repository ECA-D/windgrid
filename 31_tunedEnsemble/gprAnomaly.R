#dev.off()

library(pracma)
require(Hmisc)
library(readr)
library(matrixStats)
#library(lattice)
#library(ggplot2)
#library(sp)
#library(gstat)


setwd('~/Jouke/WindGrid/scripts/31_tunedEnsemble')

source('../91_utilities/gprFunctions.R')
source('../91_utilities/textload.R')

##############
## settings ##
##############
YEAR = 1950
MONTH = 1
relnoisey = 0.01
lagTreshold = 50e3 # maximum station-to-grid lag for prediction
thetaMin = 20e3
thetaMax = lagTreshold
nThetaBrute = 16
NstationMax = 5000 # 5000
verbose = 1
degTreshold = 0.001
plots = FALSE
useBackground0 = TRUE
corrPower = 1.0
randomSeed = 0
nfold = 20
inputDir = '/data2/Else/Jouke/WindOutput_v2/Forwardsel/'
outputDir = '/data2/Else/Jouke/WindOutput_v2/GPR_gridded/'

for(year in YEAR){
  
  ystr = toString(year)
  
  backgroundgrid = readRDS(file = paste(inputDir,"backgroundgrid_",ystr,".rds",sep=""))
  
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
    
    # load station data
    constantStations = TRUE
    if(verbose>1){print('> Reading station data for current month ...')}
    stationData = readStationData(year,month)
    
    NNstation = size(stationData$speed)[1]
    
    # IF constant stations THEN compute constant monthly analysis and prediction lags
    if(constantStations){
      xStation = cbind(stationData$lon[1:NNstation,1,drop=FALSE],stationData$lat[1:NNstation,1,drop=FALSE])
      xGrid = cbind(backgroundgrid$LON,backgroundgrid$LAT)
      
      if(verbose>1){print('> Computing constant monthly analysis lags')}
      monthLagAnalyse = krigingLag(xStation,xStation)
      
      if(verbose>1){print('> Computing constant monthly prediction mask')}
      Ngrid = size(xGrid)[1]
      MaskInGrid = matrix(NA,Ngrid,1)
      MaskInGrid4plot = matrix(NA,Ngrid,1)
      #tic()
      #monthLagPredict = matrix(NA,0,size(monthLagAnalyse)[2])
      for(k in 1:Ngrid){
        lag = krigingLag(xGrid[k,,drop=FALSE],xStation)
        if(min(lag)<lagTreshold){
          MaskInGrid[k] = TRUE
          MaskInGrid4plot[k] = 1
          #monthLagPredict = rbind(monthLagPredict,lag)
        } else {
          MaskInGrid[k] = FALSE
          MaskInGrid4plot[k] = 0
        }
      }
      nInMask = sum(MaskInGrid4plot)
      
      if(verbose>1){print('> Computing constant monthly prediction lags')}
      monthLagPredict = krigingLag(xGrid[MaskInGrid,,drop=FALSE],xStation) # can run  out of memory on local machine
      #toc()
      
    } else {
      
      stop('this part of the code should have become redundant')
      
    }
    
    # get grid storage sizes
    nMember = size(backgroundgrid$LOGBACKGROUND_MEMBER)[3]
    nDays = monthDays(as.Date(paste(ystr,'-',mstr,'-01',sep='')))
    
    # pre-allocate daily grids
    grid = list(lon=backgroundgrid$lon,lat=backgroundgrid$lat)
    grid$LON = backgroundgrid$LON
    grid$LAT = backgroundgrid$LAT
    grid$WINDSPEED_ENSEMBLE_RAW = array(NA,c(Ngrid,nMember,nDays))
    grid$WINDSPEED_ENSEMBLE_TUNED = array(NA,c(Ngrid,nMember,nDays))
    grid$WINDSPEED_ENSEMBLE_MEAN = matrix(NA,Ngrid,nDays)
    grid$WINDSPEED_ENSEMBLE_RAW_SPREAD_05_95 = matrix(NA,Ngrid,nDays)
    grid$WINDSPEED_ENSEMBLE_TUNED_SPREAD_05_95 = matrix(NA,Ngrid,nDays)
    grid$rmse_residual = matrix(NA,1,nDays)
    grid$rmse_residual_raw = matrix(NA,1,nDays)
    grid$rmse_xval = matrix(NA,1,nDays)
    grid$rmse_spread = matrix(NA,1,nDays)

    # loop over days to compute grids
    for(day in 1:nDays){
      
      if(day<10){
        dstr = paste("0",toString(day),sep="")
      } else {
        dstr = toString(day)
      }
      
      if(verbose>2){print(paste(ystr,mstr,dstr))}
      
 
      stationBackgroLog = matrix(NA,NNstation,nMember)
      Nlon = length(backgroundgrid$lon)
      Nlat = length(backgroundgrid$lat)
      for(member in 1:nMember){
        if(useBackground){
          #if(verbose>2){print('>  Computing todays anomalies (using background)')}
          #testlon = interp2(backgroundgrid$lat,backgroundgrid$lon,t(matrix(backgroundgrid$LON,Nlat,Nlon)),50,10)
          #testlat = interp2(backgroundgrid$lat,backgroundgrid$lon,t(matrix(backgroundgrid$LAT,Nlat,Nlon)),50,10)
          stationBackgroLog[,member] = interp2(backgroundgrid$lat,backgroundgrid$lon,t(matrix(backgroundgrid$LOGBACKGROUND_MEMBER[,month,member,drop=FALSE],Nlat,Nlon)),stationData$lat[,day,drop=FALSE],stationData$lon[,day,drop=FALSE])
        } else {
          if(verbose>2){print('>  Computing todays anomalies (not using background)')}
          stop('This part of the code should have become redundant')
        }
      }
      stationAnomalyLog = log(stationData$speed[,day]) - stationBackgroLog
      
      # select today's finite data
      idFinSpeed = is.finite(rowSums(stationAnomalyLog))
      xStation = cbind(stationData$lon[idFinSpeed,day,drop=FALSE],stationData$lat[idFinSpeed,day,drop=FALSE])
      yStation = stationAnomalyLog[idFinSpeed,,drop=FALSE]
      n = size(yStation)[1]
      
      # kriging: normalise
      mu = colMeans(yStation)
      st = sqrt(colVars(yStation))
      for(member in 1:nMember){
        yStation[,member] = (yStation[,member] - mu[member]) / st[member]
        
      }
      
      if(constantStations){
        
        # maximum likelihood estimate
        if(verbose>2){print('>  Todays maximum likelihood estimate')}
        thetaBrute = exp(linspace(log(thetaMin),log(thetaMax),n=nThetaBrute))
        negLogLike = NA*thetaBrute
        yStationMean = rowMeans(yStation)
        for(iterBrute in 1:nThetaBrute){
          P = exp(-.5*((monthLagAnalyse[idFinSpeed,idFinSpeed,drop=FALSE])^corrPower)/(thetaBrute[iterBrute]^corrPower))
          R = relnoisey^2 * eye(size(P)[1])
          A = R + P
          #invA = inv(A)
          #negLogLike[iterBrute] = log(det(A)) + t(yStationMean)%*%invA%*%yStationMean
          #negLogLike[iterBrute] = sum(log(eig(A))) + t(yStationMean)%*%invA%*%yStationMean
          negLogLike[iterBrute] = sum(log(eig(A))) + t(yStationMean)%*%qr.solve(A,yStationMean)
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
        P = exp(-.5*((monthLagAnalyse[idFinSpeed,idFinSpeed,drop=FALSE])^corrPower)/(theta^corrPower))
        R = relnoisey^2 * eye(size(P)[1])
        A = R + P
        b = exp(-.5*((monthLagPredict[,idFinSpeed,drop=FALSE])^corrPower)/(theta^corrPower))
        
        # kriging prediction
        yin = b %*% qr.solve(A,yStation)

        # kriging x-validation prediction
        N = size(yStation)[1]
        yinXval = matrix(NA,N,nMember)
        id0 = matrix(sample(1:N,N),N,1)
        Nout = N/nfold
        for(k in 1:nfold){
          idout = id0[ceil((k-1)*Nout+1):ceil(k*Nout)]
          yinXval[idout,] = P[idout,-idout] %*% qr.solve(A[-idout,-idout],yStation[-idout,,drop=FALSE])
        end
    }
        
      } else {
        
        stop('this part of the code should have become redundant')
        
      }
      
      # kriging: de-normalise
      gridAnomalyLog = NA*yin
      xvalAnomalyLog = NA*yinXval
      for(member in 1:nMember){
        gridAnomalyLog[,member] = mu[member] + st[member]*yin[,member,drop=FALSE]
        xvalAnomalyLog[,member] = mu[member] + st[member]*yinXval[,member,drop=FALSE]
      }
      
      if(verbose>2){print('>  Building and analysing ensembles')}
      
      # build and analyse ensemble
      if(useBackground){
        gridToday = exp( gridAnomalyLog + backgroundgrid$LOGBACKGROUND_MEMBER[MaskInGrid,month,] )
        xvalToday = exp( xvalAnomalyLog + stationBackgroLog[idFinSpeed,,drop=FALSE] )
      } else {
        stop('This part of the code should have become redundant')
        gridToday = exp( gridAnomalyLog )
      }

      # statistics: raw ensemble
      gridTodayEnsembleRaw = matrix(NA,Ngrid,nMember)
      gridTodayEnsembleRaw[MaskInGrid,] = gridToday
      gridTodayMean = matrix(NA,Ngrid,1)
      gridTodayMean[MaskInGrid] = rowMeans(gridToday)
      gridTodaySpreadRaw = matrix(NA,Ngrid,1)
      gridTodaySpreadRaw[MaskInGrid] = rowQuantiles(gridToday,probs=.95) - rowQuantiles(gridToday,probs=.05)

      # prediction errors

      stationPred = interp2(backgroundgrid$lat,backgroundgrid$lon,t(matrix(gridTodayMean,Nlat,Nlon)),stationData$lat[idFinSpeed,day,drop=FALSE],stationData$lon[idFinSpeed,day,drop=FALSE])
      Eresidual = stationPred - stationData$speed[idFinSpeed,day]
      RMSEresidual = sqrt(mean(Eresidual^2))
      if(verbose>2){print(paste('RMSE residual =',toString(RMSEresidual)))}

      idFinSpeedRaw = is.finite(stationData$speedRaw[,day]) & (stationData$speedRaw[,day]>0)
      stationPred = interp2(backgroundgrid$lat,backgroundgrid$lon,t(matrix(gridTodayMean,Nlat,Nlon)),stationData$latRaw[idFinSpeedRaw,day,drop=FALSE],stationData$lonRaw[idFinSpeedRaw,day,drop=FALSE])
      EresidualRaw = stationPred - stationData$speedRaw[idFinSpeedRaw,day]
      RMSEresidualRaw = sqrt(mean(EresidualRaw[is.finite(EresidualRaw)]^2))
      if(verbose>2){print(paste('RMSE residual raw =',toString(RMSEresidualRaw)))}
      
      Exval = rowMeans(xvalToday) - stationData$speed[idFinSpeed,day]
      RMSExval = sqrt(mean(Exval^2))
      if(verbose>2){print(paste('RMSE xval =',toString(RMSExval)))}

      statPred = interp2(backgroundgrid$lat,backgroundgrid$lon,t(matrix(gridTodaySpreadRaw,Nlat,Nlon)),stationData$lat[idFinSpeed,day,drop=FALSE],stationData$lon[idFinSpeed,day,drop=FALSE])
      RMSEspread = sqrt(mean(statPred^2))
      if(verbose>2){print(paste('RMSE spread =',toString(RMSEspread)))}

      grid$rmse_residual[day] = RMSEresidual
      grid$rmse_residual_raw[day] = RMSEresidualRaw
      grid$rmse_xval[day] = RMSExval
      grid$rmse_spread[day] = RMSEspread

      #png(filename='temp_xval.png')
      #loglog(stationSpeed[idFinSpeed,day],rowMeans(xvalToday))
      #for(k in 1:nMember){
      #  points(stationSpeed[idFinSpeed,day],xvalToday[,k],pch=20,col='green')
      #}
      #lines(stationSpeed[idFinSpeed,day],stationSpeed[idFinSpeed,day],col='red')
      #points(stationSpeed[idFinSpeed,day],rowMeans(xvalToday))
      #dev.off()

      # analyse spread: raw ensemble
      #quantileStd = ensembleRankHistStd(xvalToday,stationData$speed[idFinSpeed,day,drop=FALSE],FALSE)
      #if(verbose>2){print(paste('Raw Xval Quantile Std =',toString(quantileStd)))}

      # tune ensemble diversity
      #nTest = 31
      #factorTest = exp(linspace(log(1/30),log(30),nTest))
      #stdTest = matrix(NA,1,nTest)
      #for(k in 1:nTest){
      #  ensembleTest = ensembleDiversityFactoring(xvalToday,factorTest[k])
      #  stdTest[k] = ensembleRankHistStd(ensembleTest,stationData$speed[idFinSpeed,day,drop=FALSE],FALSE)
      #}
      #idmin = which.min(stdTest)
      #factorMin = factorTest[idmin]
      #quantileStdMin = stdTest[idmin]
      #if(verbose>2){print(paste('Tuned Xval Diversity Factor =',toString(factorMin)))}
      #if(verbose>2){print(paste('Tuned Xval Quantile Std =',toString(quantileStdMin)))}
      
      #xvalTodayTuned = ensembleDiversityFactoring(xvalToday,factorMin)
      #gridTodayTuned = ensembleDiversityFactoring(gridToday,factorMin)

      #tempTest = ensembleRankHistStd(xvalTodayTuned,stationData$speed[idFinSpeed,day,drop=FALSE],TRUE)

      #gridTodayEnsembleTuned = matrix(NA,Ngrid,nMember)
      #gridTodayEnsembleTuned[MaskInGrid,] = gridTodayTuned
      #gridTodaySpreadTuned = matrix(NA,Ngrid,1)
      #gridTodaySpreadTuned[MaskInGrid] = rowQuantiles(gridTodayTuned,probs=.95) - rowQuantiles(gridTodayTuned,probs=.05)

      if(plots){
        Nlon = length(backgroundgrid$lon)
        Nlat = length(backgroundgrid$lat)
        png(filename=paste('figures/final_grid_mean_',toString(year),'_',mstr,'_',dstr,'.png',sep=''))
        filled.contour(backgroundgrid$lon,backgroundgrid$lat,t(matrix(gridTodayMean,Nlat,Nlon)))
        dev.off()
        png(filename=paste('figures/final_grid_log10spread_',toString(year),'_',mstr,'_',dstr,'.png',sep=''))
        filled.contour(backgroundgrid$lon,backgroundgrid$lat,t(matrix(log10(gridTodaySpreadRaw),Nlat,Nlon)))
        dev.off()
        idFinSpeedRaw = is.finite(stationData$speedRaw[,day])
        png(filename=paste('figures/final_station2station_',toString(year),'_',mstr,'_',dstr,'.png',sep=''))
        statPred = interp2(backgroundgrid$lat,backgroundgrid$lon,t(matrix(gridTodayMean,Nlat,Nlon)),stationData$latRaw[idFinSpeedRaw,day,drop=FALSE],stationData$lonRaw[idFinSpeedRaw,day,drop=FALSE])
        loglog(stationData$speedRaw[idFinSpeedRaw,day],statPred,xlab='raw station observation',ylab='grid @ station')
        lines(stationData$speedRaw[idFinSpeedRaw,day],stationData$speedRaw[idFinSpeedRaw,day],col='red')
        dev.off()
        png(filename=paste('figures/final_station2stationSpread_',toString(year),'_',mstr,'_',dstr,'.png',sep=''))
        statPred = interp2(backgroundgrid$lat,backgroundgrid$lon,t(matrix(gridTodaySpreadRaw,Nlat,Nlon)),stationData$latRaw[idFinSpeedRaw,day,drop=FALSE],stationData$lonRaw[idFinSpeedRaw,day,drop=FALSE])
        loglog(stationData$speedRaw[idFinSpeedRaw,day],statPred,xlab='raw station observation',ylab='grid spread @ station')
        lines(stationData$speedRaw[idFinSpeedRaw,day],stationData$speedRaw[idFinSpeedRaw,day],col='red')
        dev.off()
        #png(filename=paste('figures/ensemble_tuning_',toString(year),'_',mstr,'_',dstr,'.png',sep=''))
        #semilogx(factorTest,stdTest,xlab='diversity factor',ylab='xval quantile freq std')
        #points(1,quantileStd,col='blue',pch=2)
        #points(factorMin,quantileStdMin,col='red',pch=2)
        #dev.off()
      }
      
      #these are the pre-allocated matrices:

      #grid$WINDSPEED_ENSEMBLE_RAW = matrix(NA,Ngrid,nMember,nDays)
      #grid$WINDSPEED_ENSEMBLE_TUNED = matrix(NA,Ngrid,nMember,nDays)
      #grid$WINDSPEED_ENSEMBLE_MEAN = matrix(NA,Ngrid,nDays)
      #rid$WINDSPEED_ENSEMBLE_RAW_SPREAD_05_95 = matrix(NA,Ngrid,nDays)
      #grid$WINDSPEED_ENSEMBLE_TUNED_SPREAD_05_95 = matrix(NA,Ngrid,nDays)

      # this is where today's values are stored:

      grid$WINDSPEED_ENSEMBLE_RAW[,,day] = gridTodayEnsembleRaw
      #grid$WINDSPEED_ENSEMBLE_TUNED[,,day] = gridTodayEnsembleTuned
      grid$WINDSPEED_ENSEMBLE_MEAN[,day] = gridTodayMean
      grid$WINDSPEED_ENSEMBLE_RAW_SPREAD_05_95[,day] = gridTodaySpreadRaw
      #grid$WINDSPEED_ENSEMBLE_TUNED_SPREAD_05_95[,day] = gridTodaySpreadTuned

      #stop('debug')
      
    }
    
    if(month<10){
      mstr = paste("0",toString(month),sep="")
    } else {
      mstr = toString(month)
    }
    if(verbose>1){print(paste('Writing grid file for year',toString(year),'month',mstr))}
    saveRDS(grid, file = paste(outputDir,"windgrid_ensemble_statistics_",toString(year),"_",mstr,".rds",sep=""))
    
  }
}
