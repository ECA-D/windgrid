#rm()
#dev.off()

library(pracma)
library(lattice)
#library(ggplot2)
library(matrixStats)

library(yaml)
suppressWarnings(library(randomFunctions)) # Prevent trivial warnings

##########################################
#                                        #
#        Get gridding paramters          #
#                                        #
##########################################

args <- getArgs()
iyaml <- yaml.load_file(args$yaml)

########################### Set Arguments ############################
gen.args <- c("lonmin","lonmax","latmin","latmax","resolution",
	      "nEnsemble","verbose","txtload_source","NstationMax",
	      "degTreshold")

forward.args <- c("lr_source","nIterRefineMin","nIterRefineMax","tol_RRMSE",
	     	  "tol_relax","power_limiter","nfold","seed","plot",
		  "xvalselect","nStationMin","dmax")

for(x in gen.args){ assign(x,iyaml[[x]]) }
for(x in forward.args){ assign(x,iyaml$forward_sel[[x]]) }

source(lr_source)

if(verbose>2){tic()}

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

for(year in YEAR){
  
  ## load data
  grid = readRDS(file = paste(inputDir,"grid_",toString(year),".rds",sep=""))
  station = readRDS(file = paste(inputDir,"station_",toString(year),".rds",sep=""))
  #Nstation = length(station$lon)
  
  ## log, sqrt transforms of variables
#  dmax = 150e3
  station$d2c[(station$d2c>dmax)] = dmax
  grid$D2C[(grid$D2C>dmax)] = dmax
  station$d2c = sqrt(station$d2c)
  grid$D2C = sqrt(grid$D2C)

  station$slope = log(station$slope+0.001)
  grid$SLOPE = log(grid$SLOPE+0.001)

  station$avg_roughness = log(10/(station$avg_roughness+0.001))
  grid$ROUGH = log(10/(grid$ROUGH+0.001))

  station$era5_avg_windspeed = log(station$era5_avg_windspeed)
  grid$ERA5SPEED = log(grid$ERA5SPEED)

  station$avg_windspeed = log(station$avg_windspeed)
  
  # initialize variables for loop
  RRMSE_month = matrix(NA,12,2)
  grid$LOGBACKGROUND_MEMBER = array(NA,dim=c(length(grid$LAT),12,nEnsemble))
  for(iter_month in 1:length(MONTH)){
    
    for(iter_ensemble in 1:nEnsemble){

      if(verbose>0){
        print(paste('== YEAR',toString(year),'MONTH',toString(iter_month),'MEMBER',toString(iter_ensemble),'=='))
      }
      
      # set local settings
      month = MONTH[iter_month]
      nfold = nfold
      seed = seed
      nIterRefineMin = nIterRefineMin
      nIterRefineMax = nIterRefineMax
      tol_RRMSE = tol_RRMSE
      tol_relax = tol_relax
      power_limiter = power_limiter
      
      ###### KEY OPERATION: define I/O ######
      X0month = cbind(station$lat[,month,drop=FALSE],
                      station$lon[,month,drop=FALSE],
                      station$era5_avg_windspeed[,month,drop=FALSE],
                      station$alt[,month,drop=FALSE],
                      station$slope[,month,drop=FALSE],
                      station$tpi[,month,drop=FALSE],
                      station$d2c[,month,drop=FALSE],
                      station$avg_roughness[,month,drop=FALSE])
      X0grid = cbind(grid$LAT,
                     grid$LON,
                     grid$ERA5SPEED[,month,drop=FALSE],
                     grid$ALT,
                     grid$SLOPE,
                     grid$TPI,
                     grid$D2C,
                     grid$ROUGH[,month,drop=FALSE])
      speed = station$avg_windspeed[,month,drop=FALSE]
      #######################################
      
      # remove NANs
      idfin = is.finite(speed)
      X0month = X0month[idfin,,drop=FALSE]
      speed = speed[idfin,,drop=FALSE]
      idfin = is.finite(rowMeans(X0month))
      X0month = X0month[idfin,,drop=FALSE]
      speed = speed[idfin,,drop=FALSE]

      if(length(speed)<nStationMin){
        if(verbose>0){print('Number of stations to small, skipped this month')}
        break
      }
      
      # bootstrap
      set.seed(100*iter_ensemble)
      N = length(speed)
      idBoot = matrix(sample(1:N,N,replace=TRUE),N,1)
      X0month = X0month[idBoot,,drop=FALSE]
      speed = speed[idBoot,,drop=FALSE]
 
      ###### KEY OPERATION: initialize LR ######
      p = rbind(c(0,0,0,0,0,0,0,0),c(0,0,0,0,0,0,1,0),c(0,0,0,0,0,0,0,1))
      ##########################################
      
      ## forward selection
      RRMSE = matrix(NA,nIterRefineMax,2)
      np0 = size(p)[1]
      for(iter_refine in 1:nIterRefineMax){
        if(verbose>1){
          print(paste("== year",toString(year),"month",toString(month),' >> member',toString(iter_ensemble),">> forward selection >> iteration",toString(iter_refine),"of",toString(nIterRefineMax),"=="))
        }
        np = size(p)[1]
        ncovar = size(p)[2]
        
        ###### KEY OPERATION: errors for current design ######
        error = lrRRMSE(X0month,p,speed,nfold,seed,TRUE)
        RRMSE[iter_refine,1] = error$training
        RRMSE[iter_refine,2] = error$xval
        ######################################################
        
        if(verbose>1){
          if(iter_refine>1){
            print(paste("previous training RRMSE =",toString(RRMSE[iter_refine-1,1])))
          }
          print(paste("current  training RRMSE =",toString(RRMSE[iter_refine,1])))
          if(iter_refine>1){
            print(paste("best prv x-valida RRMSE =",toString(min(RRMSE[1:(iter_refine-1),2]))))
          }
          print(paste("current  x-valida RRMSE =",toString(RRMSE[iter_refine,2])))
        }
        if(iter_refine>=nIterRefineMin){
          RRMSEcopy = RRMSE
          RRMSEcopy[1:nIterRefineMin,2] = 1e3
          if((RRMSE[iter_refine,2]-min(RRMSE[(1:(iter_refine-1)),2]))>=tol_RRMSE){
            if(tol_relax<1){
              if(verbose>1){
                print(">>>> RRMSE tolerance reached: break <<<<")
              }
              break
            } else {
              if(verbose>1){
                print(paste(">>>> RRMSE tolerance reached: dropping tol_relax by 1 (",toString(tol_relax),") <<<<",sep=''))
              }
              tol_relax = tol_relax-1
            }
          } else {
            tol_relax = tol_relax # reset counter
          }
        }
        RRMSE_i_test = matrix(Inf,np*ncovar,1)
        pnext_store = matrix(0,np*ncovar,ncovar)
        iter_test = 1
        progress = 1
        if(iter_refine<nIterRefineMax){
          X0test = X0month
          speedtest = speed
          #cat('Countdown: ')
          for(k in 1:np){
            for(l in 1:ncovar){
              #print(paste("> iteration",toString(iter_refine),"trying p-increase",toString(iter_test),"of",toString(np*ncovar),"..."))
              
              ###### KEY OPERATION: generate a candidate ######
              pnext = p[k,,drop=FALSE]
              pnext[,l] = pnext[,l] + 1
              pnext_store[iter_test,] = pnext
              pInvestigate = TRUE
              for(m in 1:np){
                if(sum(abs(p[m,,drop=FALSE]-pnext))==0){
                  pInvestigate = FALSE
                }
                if(sum(pnext>power_limiter)>0){
                  pInvestigate = FALSE
                }
                if(sum(abs(pnext[2:size(pnext)[2]]))>0){
                  if(sum(abs(pnext[1:2]))>0){
                    pInvestigate = FALSE
                  }
                }
              }
              ################################################
              
              if(pInvestigate){
                
                ptest = rbind(p,pnext)
                
                ###### KEY OPERATION: errors for this candidate ######
		if(xvalselect){
                  error = lrRRMSE(X0test,ptest,speedtest,nfold,seed,TRUE)
                  RRMSE_i_test[iter_test,1] = error$xval
                } else {
                  error = lrRRMSE(X0test,ptest,speedtest,NA,NA,FALSE)
                  RRMSE_i_test[iter_test,1] = error$training
                }
                ######################################################
                
              }
              if(ceil(10*iter_test/np/ncovar)>=progress){
                #cat(paste(10-progress,''))
                progress = progress+1
              }
              iter_test = iter_test+1
            }
          }
          #cat('\n')
          if(min(RRMSE_i_test)==Inf){
            if(verbose>1){
              print('>>>>> No valid candidates > break')
            }
            break
          }
          
          ###### KEY OPERATION: greedy selection ########
          if(verbose>1){
            print('>> greedy selection of best new basis vector ...')
          }
          idgreedy = which.min(RRMSE_i_test[,1])
          p = rbind(p,pnext_store[idgreedy,,drop=FALSE])
          if(verbose>1){
            print('>> appended new basis vector')
          }
          ################################################
          
        }
      }
      
      ###### KEY OPERATION: select minimal x-val ########
      if(verbose>1){
        print('>>>> Selecting minimal x-validation ...')
      }
      # find minimal x-validation error
      if(nIterRefineMax>1){
        RRMSEcopy = RRMSE
        RRMSEcopy[1:nIterRefineMin,2] = Inf
        iter_min = which.min(RRMSEcopy[1:iter_refine,2])
      } else {
        iter_min = 1
      }
      # drop part of p that leads to over fitting
      p = p[1:(iter_min+np0-1),,drop=FALSE]
      ###################################################
      
      # collect data
      RRMSE_month[month,1:2] = RRMSE[iter_min,1:2]
      
      # convergence plot
      if(plot){
      png(filename=paste('figures/convergence_month_',toString(year),'_',toString(month),'_member_',toString(iter_ensemble),'.png',sep=''))
      iteration = 1:nIterRefineMax
      plot(iteration,RRMSE[,1],type='b',col=3,ylim=c(0,1))
      lines(iteration,RRMSE[,2],type='b',col=2)
      points(iter_min,RRMSE[iter_min,2],col=1)
      legend(x = "bottomright",          # Position
             legend = c("training", "x-val", "best"),  # Legend texts
             lty = c(1, 1, 1),           # Line types
             col = c(3, 2, 1),           # Line colors
             lwd = 1)                 # Line width
      grid()
      dev.off()
      }
      
      # power matrix plot
      if(plot){
        png(filename=paste('figures/powerMatrix_month_',toString(year),'_',toString(month),'_member_',toString(iter_ensemble),'.png',sep=''))
        pp = matrix(NA,nIterRefineMax+np0-1,size(p)[2])
        pp[1:size(p)[1],] = p
        colnames(pp) = c('lat','lon','era5','alt','slp','tpi','d2c','rgh')
        print(levelplot(t(pp),main=paste('LR power matrix for month',toString(month)),ylab='LR basis vector',xlab='covariate',aspect="fill"))
        dev.off()
      }
      
      # signal
      if(plot){
        png(filename=paste('figures/speedHist_month_',toString(year),'_',toString(month),'_member_',toString(iter_ensemble),'.png',sep=''))
        windspeed = exp(speed)
        hist(windspeed,xlim=c(0,10),breaks=11)
        dev.off()
      }
      
      # residual
      if(plot){
        png(filename=paste('figures/resHist_month_',toString(year),'_',toString(month),'_member_',toString(iter_ensemble),'.png',sep=''))
        residual = exp(linearRegression(X0month,p,speed,X0month))-exp(speed)
        hist(residual,xlim=c(-5,5),breaks=11)
        dev.off()
      }
      
      # background grid
      Nlon = length(grid$lon)
      Nlat = length(grid$lat)
      speedgrid = linearRegression(X0month,p,speed,X0grid)
      grid$LOGBACKGROUND_MEMBER[,month,iter_ensemble] = speedgrid
      if(plot){
        png(filename=paste('figures/grid_logspeed_month_',toString(year),'_',toString(month),'_member_',toString(iter_ensemble),'.png',sep=''))
        filled.contour(grid$lon,grid$lat,t(matrix(grid$LOGBACKGROUND_MEMBER[,month,iter_ensemble],Nlat,Nlon)))
        title(paste('LR background log wind speed for month',toString(month)))
        dev.off()
        png(filename=paste('figures/grid_month_',toString(year),'_',toString(month),'_member_',toString(iter_ensemble),'.png',sep=''))
        filled.contour(grid$lon,grid$lat,t(matrix(exp(grid$LOGBACKGROUND_MEMBER[,month,iter_ensemble]),Nlat,Nlon)))
        title(paste('LR background log wind speed for month',toString(month)))
        dev.off()
      }
    }
  }

  saveRDS(grid, file =paste(outputDir,"backgroundgrid_",toString(year),".rds",sep=""))

}

if(verbose>2){toc()}
