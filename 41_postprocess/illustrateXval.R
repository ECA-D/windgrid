library(pracma)
require(Hmisc)
library(ncdf4)

setwd('~/Jouke/WindGrid/scripts/41_postprocess')

##############
## settings ##
##############
YEAR = 1950:2020
MONTH = 1:12
verbose = 1
windgridInDir = '/data2/Else/Jouke/WindOutput_v2/GPR_gridded/'


########################
## read & loops ##
########################
yearScale = matrix(NA,length(YEAR)*length(MONTH),1)
monthScale = matrix(NA,length(YEAR)*length(MONTH),1)
monthlyRMSExval = matrix(NA,length(YEAR)*length(MONTH),1)
timeCounter = 1
for(year in YEAR){
  
  ystr = toString(year)
  
  for(month in MONTH){
    
    if(month<10){
      mstr = paste("0",toString(month),sep="")
    } else {
      mstr = toString(month)
    }
    
    if(verbose>0){print(paste('== year',ystr,'-> month',mstr,'=='))}
    nDays = monthDays(as.Date(paste(ystr,'-',mstr,'-01',sep='')))
    
    # read windfield statistics
    if(verbose>1){print('Reading windgrid statistics')}
    gridStatistics = readRDS(file = paste(windgridInDir,"/windgrid_ensemble_statistics_",toString(year),"_",mstr,".rds",sep=""))
    
    yearScale[timeCounter] = year + (month-.5)/12
    monthScale[timeCounter] = month
    monthlyRMSExval[timeCounter] = sqrt(mean(gridStatistics$rmse_xval))
    if(verbose>1){print(monthlyRMSExval[timeCounter])}

    timeCounter = timeCounter+1

  }
}

X = cbind(1+0*yearScale,yearScale)
p = qr.solve(X,log(monthlyRMSExval))
yi = exp(X%*%p)

png('rmseXval_vs_year.png')
semilogy(yearScale,monthlyRMSExval,xlab='year (1 Jan)',ylab='x-Validation RMSE [m/s]')
lines(yearScale,yi,col='red')
dev.off()

png('rmseXval_vs_year_noOutliers.png')
semilogy(yearScale,monthlyRMSExval,xlab='year (1 Jan)',ylab='x-Validation RMSE [m/s]',ylim=c(.5,5))
lines(yearScale,yi,col='red')
dev.off()

relRMSE = monthlyRMSExval/yi

png('rmseXval_vs_month.png')
plot(monthScale,relRMSE,xlab='month',ylab='relative x-Validation RMSE [m/s]',ylim=c(0.5,1.5))
dev.off()
