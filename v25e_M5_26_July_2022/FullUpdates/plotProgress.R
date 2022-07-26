#dev.off()

library(pracma)
#require(Hmisc)
#library(readr)
#library(matrixStats)
#library(lattice)
#library(ggplot2)
#library(sp)
#library(gstat)
library(yaml)
suppressWarnings(library(randomFunctions)) # Prevent trivial warnings

#yamlfile="/home/besselaa/Jouke/WindGrid/scripts/FullUpdates/eobs_fg_no_homog.yaml"
progdir="/home/besselaa/Jouke/WindGrid/scripts/FullUpdates"
#outputDir="/data3/Else/EOBSv24.0e/Grid_0.1deg/fg/Tempfiles"

#iyaml <- yaml.load_file(yamlfile)

##########################################
#                                        #
#        Get gridding paramters          #
#                                        #
##########################################

print("command line options: yaml=<yaml-file> year=\"<year1>(<;year2>)\" difftime=<hours>")

args <- getArgs()
iyaml <- yaml.load_file(args$yaml)

gen.args <- c("scratch")
postprocess.args <- c("final_dir_10","final_dir_25")

for(x in gen.args){ assign(x,iyaml[[x]]) }
for(x in postprocess.args){ assign(x,iyaml$postprocess[[x]]) }

# Define dir with files to check
if(is.null(args$scratch)){
    outputDir <- scratch
} else {
    outputDir <- args$scratch
}

# Define years
if(is.null(args$year)){
   YEAR <- 1980:2021
   yearstart=min(YEAR)-1
   yearend=max(YEAR)+1
} else {
    yrsplit <- strsplit(args$year,";")[[1]]
    YEAR <- as.numeric(yrsplit)
     if(length(YEAR)==1) {
       yearstart=YEAR-2
       yearend=YEAR+2
     }else{
       yearstart=YEAR[1]-2
       yearend=YEAR[2]+2
     }
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

# Define how long ago files could have been created
if(is.null(args$difftime)) {
    difftime <- 12 # 12 hours
} else {
    difftime <- as.numeric(args$difftime)
}

print(paste0("Checking for files less than ",difftime," hours old"))

#create.time <- file.info("/data3/Else/EOBSv24.0e/Grid_0.1deg/fg/Tempfiles/grid_1950.rds")$ctime
time.today <- Sys.time()


#as.numeric(difftime(time.today,create.time,"%H"))

plot.col1='grey'

png(paste0(progdir,'/ProgressPlots/preprocess.png'))
plot(c(1,1),c(yearstart,format(time.today,"%Y")),ylab='Year',col='white',main='preprocess')
for(year in YEAR) {
    plot.col2 <- 'white'

    filename1 = paste0(outputDir,"/grid_",toString(year),".rds")
    if(file.exists(filename1)){
      create.time <- file.info(filename1)$ctime
      diff.hours <- as.numeric(difftime(time.today,create.time,units='hours'))
      if(diff.hours <= difftime){plot.col2 <- 'green'}else{plot.col2 <- 'white'}
    }
    points(1.,year,col=plot.col2, pch=16)
    points(1.,year,col=plot.col1, pch=1)
}

png(paste0(progdir,'/ProgressPlots/forward_selection.png'))
plot(c(1,1),c(yearstart,format(time.today,"%Y")),ylab='Year',col='white',main='forward selection')
for(year in YEAR){
    plot.col2 <- 'white'

    filename2 = paste0(outputDir,"/backgroundgrid_",toString(year),".rds")
    if(file.exists(filename2)){
      create.time <- file.info(filename2)$ctime
      diff.hours <- as.numeric(difftime(time.today,create.time,units='hours'))
      if(diff.hours <= difftime){plot.col2 <- 'green'}else{plot.col2 <- 'white'}
    }
    points(1.,year,col=plot.col2, pch=16)
    points(1.,year,col=plot.col1, pch=1)
}    


png(paste0(progdir,'/ProgressPlots/gpr_anomaly.png'))
plot(c(1,12),c(yearstart,yearend),xlab='Month',ylab='Year',col='white',main='GPR anomaly')
for(year in YEAR){
  for(month in MONTH){
    if(month<10){
      mstr = paste("0",toString(month),sep="")
    } else {
      mstr = toString(month)
    }
    plot.col2 <- 'white'

    filename3 = paste0(outputDir,"/windgrid_ensemble_statistics_",toString(year),"_",mstr,".rds")
    if(file.exists(filename3)){
      create.time <- file.info(filename3)$ctime
      diff.hours <- as.numeric(difftime(time.today,create.time,units='hours'))
      if(diff.hours <= difftime){plot.col2 <- 'green'}else{plot.col2 <- 'white'}
    }
    points(month,year,col=plot.col2, pch=16)
    points(month,year,col=plot.col1, pch=1)
  }
}   





