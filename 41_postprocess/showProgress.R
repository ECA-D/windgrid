#dev.off()

library(pracma)
#require(Hmisc)
#library(readr)
#library(matrixStats)
#library(lattice)
#library(ggplot2)
#library(sp)
#library(gstat)


setwd('~/Jouke/WindGrid/scripts/41_postprocess')


png('progress.png')
plot(c(1,12),c(1950,2021),xlab='Month',ylab='Year',col='white')

for(year in c(1950,1981,2001,2011,2021)){
  lines(c(-1,14),c(year-.5,year-.5),col='grey')
}

for(year in 1950:2021){

  for(month in 1:12){

    if(month<10){
      mstr = paste("0",toString(month),sep="")
    } else {
      mstr = toString(month)
    }

    filename = paste("/data2/Else/Jouke/WindOutput_v2/GPR_gridded/windgrid_ensemble_statistics_",toString(year),"_",mstr,".rds",sep="")

    if(file.exists(filename)){
      points(month,year,col='green',pch=16)
    } else {
      points(month,year,col='grey')
    }
  }
}




