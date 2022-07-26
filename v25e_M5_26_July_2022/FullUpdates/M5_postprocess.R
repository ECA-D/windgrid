library(pracma)
require(Hmisc)
library(ncdf4)
library(lubridate)
library(yaml)
suppressWarnings(library(randomFunctions)) # Prevent trivial warnings

##########################################
#                                        #
#        Get gridding paramters          #
#                                        #
##########################################

args <- getArgs()
iyaml <- yaml.load_file(args$yaml)

Nmember = 20
outlierTreshold = 35

########################### Set Arguments ############################
gen.args <- c("lonmin","lonmax","latmin","latmax","resolution",
	      "nEnsemble","verbose","txtload_source","NstationMax",
	      "degTreshold")

postprocess.args <- c("final_dir_10","final_dir_25","plots","version",
		     "scalefactor")

for(x in gen.args){ assign(x,iyaml[[x]]) }
for(x in postprocess.args){ assign(x,iyaml$postprocess[[x]]) }

## Define outputDir directory for 0.1deg
if(is.null(args$final_dir)){
    outputDir_10 <- iyaml$postprocess$final_dir_10
} else {
    outputDir_10 <- args$final_dir_10
}

## Define outputDir directory for 0.25deg
if(is.null(args$final_dir)){
    outputDir_25 <- iyaml$postprocess$final_dir_25
} else {
    outputDir_25 <- args$final_dir_25
}

## Define inputDir directory
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

iyear <- YEAR[1]

########################
## read & write loops ##
########################
print('Starting year:')
print(min(YEAR))
print('Ending year:')
print(max(YEAR))
print('Input dir:')
print(inputDir)
print('Output dir 10 km:')
print(outputDir_10)
print('Output dir 25 km:')
print(outputDir_25)

for(year in YEAR){

  timeCounter = 1

  # get coordinates
  grid0 = readRDS(file = paste0(inputDir,"/windgrid_ensemble_statistics_",year,"_01.rds",sep=""))
  lonvals10km = grid0$lon
  latvals10km = grid0$lat
  rm(grid0)
  #lonvals25km = seq(min(lonvals10km-0.05+0.125),max(lonvals10km+0.05-0.125),.25) # to check -> wrong output coordinates compared to Tmin
  lonvals25km = seq(min(lonvals10km+0.05+0.125),max(lonvals10km+0.05-0.125),.25) # to check -> OK
  latvals25km = seq(min(latvals10km-0.05+0.125),max(latvals10km+0.05-0.125),.25) # to check -> OK

  # define spatial dims
  lon10km <- ncdim_def('longitude','degrees_east', lonvals10km, unlim=FALSE,create_dimvar=TRUE,calendar=NA,longname='Longitude values')
  lat10km <- ncdim_def('latitude','degrees_north', latvals10km, unlim=FALSE,create_dimvar=TRUE,calendar=NA,longname='Latitude values')
  lon25km <- ncdim_def('longitude','degrees_east', lonvals25km, unlim=FALSE,create_dimvar=TRUE,calendar=NA,longname='Longitude values')
  lat25km <- ncdim_def('latitude','degrees_north', latvals25km, unlim=FALSE,create_dimvar=TRUE,calendar=NA,longname='Latitude values')

  # define time dim
  origin = as.Date('1950-01-01')
  day1 = as.Date(paste(toString(year),'-',toString(min(MONTH)),'-01',sep=''))
  nDays = monthDays(as.Date(paste(toString(year),'-',toString(max(MONTH)),'-01',sep='')))
  day2 = as.Date(paste(toString(year),'-',toString(max(MONTH)),'-',toString(nDays),sep=''))
  time1 = difftime(day1,origin, units = "days")
  time2 = difftime(day2,origin, units = "days")
  time <- ncdim_def('time','days since 1950-01-01 00:00:00',time1:time2, unlim=TRUE,create_dimvar=TRUE,calendar="standard",longname='Time in days')
  member <- ncdim_def('member','ensemble member',1:Nmember,unlim=FALSE,create_dimvar=TRUE,calendar=NA,longname='Ensemble member')

  # define vars
  meanSpeed10km <- ncvar_def('fg','m/s',dim=list(lon10km,lat10km,time),missval=-9999.,longname='Ensemble mean wind speed',prec="short",compression=1)
  spreadSpeed10km <- ncvar_def('fg','m/s',dim=list(lon10km,lat10km,time),missval=-9999.,longname='Ensemble spread wind speed',prec="short",compression=1)
  ensSpeed10km <- ncvar_def('fg','m/s',dim=list(lon10km,lat10km,member,time),missval=-9999.,longname='Ensemble members wind speed',prec="short",compression=1)
  meanSpeed25km <- ncvar_def('fg','m/s',dim=list(lon25km,lat25km,time),missval=-9999.,longname='Ensemble mean wind speed',prec="short",compression=1)
  spreadSpeed25km <- ncvar_def('fg','m/s',dim=list(lon25km,lat25km,time),missval=-9999.,longname='Ensemble spread wind speed',prec="short",compression=1)
  ensSpeed25km <- ncvar_def('fg','m/s',dim=list(lon25km,lat25km,member,time),missval=-9999.,longname='Ensemble members wind speed',prec="short",compression=1)

  # open files for writing
  ncnewMean10km <- nc_create(paste(outputDir_10,"/fg_ens_mean_0.1deg_reg_",year,".nc",sep=""),meanSpeed10km)
  ncnewSpread10km <- nc_create(paste(outputDir_10,"/fg_ens_spread_0.1deg_reg_",year,".nc",sep=""),spreadSpeed10km)
  ncnewEnsMem10km <- nc_create(paste(outputDir_10,"/fg_ens_members_0.1deg_reg_",year,".nc",sep=""),ensSpeed10km)
  ncnewMean25km <- nc_create(paste(outputDir_25,"/fg_ens_mean_0.25deg_reg_",year,".nc",sep=""),meanSpeed25km)
  ncnewSpread25km <- nc_create(paste(outputDir_25,"/fg_ens_spread_0.25deg_reg_",year,".nc",sep=""),spreadSpeed25km)
  ncnewEnsMem25km <- nc_create(paste(outputDir_25,"/fg_ens_members_0.25deg_reg_",year,".nc",sep=""),ensSpeed25km)

  # add standard names
  ncatt_put(ncnewMean10km,'longitude',"standard_name",'longitude')
  ncatt_put(ncnewSpread10km,'longitude',"standard_name",'longitude')
  ncatt_put(ncnewEnsMem10km,'longitude',"standard_name",'longitude')
  ncatt_put(ncnewMean25km,'longitude',"standard_name",'longitude')
  ncatt_put(ncnewSpread25km,'longitude',"standard_name",'longitude')
  ncatt_put(ncnewEnsMem25km,'longitude',"standard_name",'longitude')

  ncatt_put(ncnewMean10km,'latitude',"standard_name",'latitude')
  ncatt_put(ncnewSpread10km,'latitude',"standard_name",'latitude')
  ncatt_put(ncnewEnsMem10km,'latitude',"standard_name",'latitude')
  ncatt_put(ncnewMean25km,'latitude',"standard_name",'latitude')
  ncatt_put(ncnewSpread25km,'latitude',"standard_name",'latitude')
  ncatt_put(ncnewEnsMem25km,'latitude',"standard_name",'latitude')

  ncatt_put(ncnewMean10km,'longitude',"axis",'X')
  ncatt_put(ncnewSpread10km,'longitude',"axis",'X')
  ncatt_put(ncnewEnsMem10km,'longitude',"axis",'X')
  ncatt_put(ncnewMean25km,'longitude',"axis",'X')
  ncatt_put(ncnewSpread25km,'longitude',"axis",'X')
  ncatt_put(ncnewEnsMem25km,'longitude',"axis",'X')

  ncatt_put(ncnewMean10km,'latitude',"axis",'Y')
  ncatt_put(ncnewSpread10km,'latitude',"axis",'Y')
  ncatt_put(ncnewEnsMem10km,'latitude',"axis",'Y')
  ncatt_put(ncnewMean25km,'latitude',"axis",'Y')
  ncatt_put(ncnewSpread25km,'latitude',"axis",'Y')
  ncatt_put(ncnewEnsMem25km,'latitude',"axis",'Y')

  ncatt_put(ncnewMean10km,'fg',"standard_name",'wind_speed')
  ncatt_put(ncnewSpread10km,'fg',"standard_name",'wind_speed')
  ncatt_put(ncnewEnsMem10km,'fg',"standard_name",'wind_speed')
  ncatt_put(ncnewMean25km,'fg',"standard_name",'wind_speed')
  ncatt_put(ncnewSpread25km,'fg',"standard_name",'wind_speed')
  ncatt_put(ncnewEnsMem25km,'fg',"standard_name",'wind_speed')

  ncatt_put(ncnewMean10km,'time',"standard_name",'time')
  ncatt_put(ncnewSpread10km,'time',"standard_name",'time')
  ncatt_put(ncnewEnsMem10km,'time',"standard_name",'time')
  ncatt_put(ncnewMean25km,'time',"standard_name",'time')
  ncatt_put(ncnewSpread25km,'time',"standard_name",'time')
  ncatt_put(ncnewEnsMem25km,'time',"standard_name",'time')

  ncatt_put(ncnewMean10km,'fg',"missing_value",'-9999s')
  ncatt_put(ncnewSpread10km,'fg',"missing_value",'-9999s')
  ncatt_put(ncnewEnsMem10km,'fg',"missing_value",'-9999s')
  ncatt_put(ncnewMean25km,'fg',"missing_value",'-9999s')
  ncatt_put(ncnewSpread25km,'fg',"missing_value",'-9999s')
  ncatt_put(ncnewEnsMem25km,'fg',"missing_value",'-9999s')

  # add attributes
  ncatt_put(ncnewMean10km,'fg',"scale_factor",scalefactor,prec="float")
  ncatt_put(ncnewMean10km,'fg',"add_offset",0.0,prec="float")
  ncatt_put(ncnewSpread10km,'fg',"scale_factor",scalefactor,prec="float")
  ncatt_put(ncnewSpread10km,'fg',"add_offset",0.0,prec="float")
  ncatt_put(ncnewEnsMem10km,'fg',"scale_factor",scalefactor,prec="float")
  ncatt_put(ncnewEnsMem10km,'fg',"add_offset",0.0,prec="float")
  ncatt_put(ncnewMean25km,'fg',"scale_factor",scalefactor,prec="float")
  ncatt_put(ncnewMean25km,'fg',"add_offset",0.0,prec="float")
  ncatt_put(ncnewSpread25km,'fg',"scale_factor",scalefactor,prec="float")
  ncatt_put(ncnewSpread25km,'fg',"add_offset",0.0,prec="float")
  ncatt_put(ncnewEnsMem25km,'fg',"scale_factor",scalefactor,prec="float")
  ncatt_put(ncnewEnsMem25km,'fg',"add_offset",0.0,prec="float")

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
    gridStatistics = readRDS(file = paste(inputDir,"/windgrid_ensemble_statistics_",toString(year),"_",mstr,".rds",sep=""))

    # loop over days
    if(verbose>1){print('Writing windgrid NC files')}
    for(day in 1:nDays){
      
      if(day<10){
        dstr = paste("0",toString(day),sep="")
      } else {
        dstr = toString(day)
      }
      
      if(verbose>2){print(paste(ystr,mstr,dstr))}
      
      # todays grid
      lon = gridStatistics$lon
      lat = gridStatistics$lat
      Nlon = length(lon)
      Nlat = length(lat)
      meanToday10km = t(matrix(gridStatistics$WINDSPEED_ENSEMBLE_MEAN[,day],Nlat,Nlon))
      spreadToday10km = t(matrix(gridStatistics$WINDSPEED_ENSEMBLE_TUNED_SPREAD_QUANT[,day],Nlat,Nlon))
      ensembleToday10km = array(0,dim=c(Nlon,Nlat,Nmember))
      for(member in 1:Nmember){
        ensembleToday10km[,,member] = t(matrix(gridStatistics$WINDSPEED_ENSEMBLE_TUNED[,member,day],Nlat,Nlon))
      }
      
      # filter grid (outlier detection)
      idout = (meanToday10km<=0)
      meanToday10km[idout] = NA
      spreadToday10km[idout] = NA
      for(member in 1:Nmember){
        ensTemp = ensembleToday10km[,,member]
        ensTemp[idout] = NA
        ensembleToday10km[,,member] = ensTemp
      }

      idout = (meanToday10km>outlierTreshold)
      meanToday10km[idout] = NA
      spreadToday10km[idout] = NA
      for(member in 1:Nmember){
        ensTemp = ensembleToday10km[,,member]
        ensTemp[idout] = NA
        ensembleToday10km[,,member] = ensTemp
      }
      
      idout = (spreadToday10km<=0)
      meanToday10km[idout] = NA
      spreadToday10km[idout] = NA
      for(member in 1:Nmember){
        ensTemp = ensembleToday10km[,,member]
        ensTemp[idout] = NA
        ensembleToday10km[,,member] = ensTemp
      }
      
      idout = (spreadToday10km>outlierTreshold)
      meanToday10km[idout] = NA
      spreadToday10km[idout] = NA
      for(member in 1:Nmember){
        ensTemp = ensembleToday10km[,,member]
        ensTemp[idout] = NA
        ensembleToday10km[,,member] = ensTemp
      }

      ensMin = apply(ensembleToday10km,c(1,2),min)
      idout = (ensMin<=0)
      meanToday10km[idout] = NA
      spreadToday10km[idout] = NA
      for(member in 1:Nmember){
        ensTemp = ensembleToday10km[,,member]
        ensTemp[idout] = NA
        ensembleToday10km[,,member] = ensTemp
      }

      ensMax = apply(ensembleToday10km,c(1,2),max)
      idout = (ensMax>outlierTreshold)
      meanToday10km[idout] = NA
      spreadToday10km[idout] = NA
      for(member in 1:Nmember){
        ensTemp = ensembleToday10km[,,member]
        ensTemp[idout] = NA
        ensembleToday10km[,,member] = ensTemp
      }

      # interpolate to 25 km
      COOR = meshgrid(lonvals25km,latvals25km)
      LON = COOR$X
      LAT = COOR$Y
      meanToday25km = t(matrix(interp2(latvals10km,lonvals10km,meanToday10km,LAT,LON),length(latvals25km),length(lonvals25km)))
      spreadToday25km = t(matrix(interp2(latvals10km,lonvals10km,spreadToday10km,LAT,LON),length(latvals25km),length(lonvals25km)))
      ensembleToday25km = array(0,dim=c(length(lonvals25km),length(latvals25km),Nmember))
      for(member in 1:Nmember){
        ensembleToday25km[,,member] = t(matrix(interp2(latvals10km,lonvals10km,ensembleToday10km[,,member],LAT,LON),length(latvals25km),length(lonvals25km)))
      }
      
      # append to 10 km netcdf files
      ncvar_put(ncnewMean10km,meanSpeed10km,meanToday10km/scalefactor,start=c(1,1,timeCounter),count=c(length(lonvals10km),length(latvals10km),1))
      ncvar_put(ncnewSpread10km,spreadSpeed10km,spreadToday10km/scalefactor,start=c(1,1,timeCounter),count=c(length(lonvals10km),length(latvals10km),1))
      for(member in 1:Nmember){
        ncvar_put(ncnewEnsMem10km,ensSpeed10km,ensembleToday10km[,,member]/scalefactor,start=c(1,1,member,timeCounter),count=c(length(lonvals10km),length(latvals10km),1,1))
      }
      
      # append to 25 km netcdf files
      ncvar_put(ncnewMean25km,meanSpeed25km,meanToday25km/scalefactor,start=c(1,1,timeCounter),count=c(length(lonvals25km),length(latvals25km),1))
      ncvar_put(ncnewSpread25km,spreadSpeed25km,spreadToday25km/scalefactor,start=c(1,1,timeCounter),count=c(length(lonvals25km),length(latvals25km),1))
      for(member in 1:Nmember){
        ncvar_put(ncnewEnsMem25km,ensSpeed25km,ensembleToday25km[,,member]/scalefactor,start=c(1,1,member,timeCounter),count=c(length(lonvals25km),length(latvals25km),1,1))
      }
      
      timeCounter = timeCounter + 1
      
      if(plots){
        png(filename=paste('figures/10km_final_grid_mean_',toString(year),'_',toString(month),'_',toString(day),'.png',sep=''))
        filled.contour(lonvals10km,latvals10km,meanToday10km)
        dev.off()
        png(filename=paste('figures/10km_final_grid_spread_',toString(year),'_',toString(month),'_',toString(day),'.png',sep=''))
        filled.contour(lonvals10km,latvals10km,spreadToday10km)
        dev.off()
        png(filename=paste('figures/25km_final_grid_mean_',toString(year),'_',toString(month),'_',toString(day),'.png',sep=''))
        filled.contour(lonvals25km,latvals25km,meanToday25km)
        dev.off()
        png(filename=paste('figures/25km_final_grid_spread_',toString(year),'_',toString(month),'_',toString(day),'.png',sep=''))
        filled.contour(lonvals25km,latvals25km,spreadToday25km)
        dev.off()
      }


    }
  }


  # close nc files
  nc_close(ncnewMean10km)
  nc_close(ncnewSpread10km)
  nc_close(ncnewEnsMem10km)
  nc_close(ncnewMean25km)
  nc_close(ncnewSpread25km)
  nc_close(ncnewEnsMem25km)

  # add global attributes
  print("adding global attributes")

  filename <- paste0(outputDir_10,"/fg_ens_mean_0.1deg_reg_",year,".nc")
  system(sprintf("ncatted -h -O -a history,global,d,, %s",filename))
  system(sprintf("ncatted -O -h -a E-OBS_version,global,a,c,%s %s",version,filename))
  system(sprintf("ncatted -O -h -a References,global,a,c,'http://surfobs.climate.copernicus.eu/dataaccess/access_eobs.php' %s",filename))

  filename <- paste0(outputDir_25,"/fg_ens_mean_0.25deg_reg_",year,".nc")
  system(sprintf("ncatted -h -O -a history,global,d,, %s",filename))
  system(sprintf("ncatted -O -h -a E-OBS_version,global,a,c,%s %s",version,filename))
  system(sprintf("ncatted -O -h -a References,global,a,c,'http://surfobs.climate.copernicus.eu/dataaccess/access_eobs.php' %s",filename))

  filename <- paste0(outputDir_10,"/fg_ens_spread_0.1deg_reg_",year,".nc")
  system(sprintf("ncatted -h -O -a history,global,d,, %s",filename))
  system(sprintf("ncatted -O -h -a E-OBS_version,global,a,c,%s %s",version,filename))
  system(sprintf("ncatted -O -h -a References,global,a,c,'http://surfobs.climate.copernicus.eu/dataaccess/access_eobs.php' %s",filename))

  filename <- paste0(outputDir_25,"/fg_ens_spread_0.25deg_reg_",year,".nc")
  system(sprintf("ncatted -h -O -a history,global,d,, %s",filename))
  system(sprintf("ncatted -O -h -a E-OBS_version,global,a,c,%s %s",version,filename))
  system(sprintf("ncatted -O -h -a References,global,a,c,'http://surfobs.climate.copernicus.eu/dataaccess/access_eobs.php' %s",filename))

  filename <- paste0(outputDir_10,"/fg_ens_members_0.1deg_reg_",year,".nc")
  system(sprintf("ncatted -h -O -a history,global,d,, %s",filename))
  system(sprintf("ncatted -O -h -a E-OBS_version,global,a,c,%s %s",version,filename))
  system(sprintf("ncatted -O -h -a References,global,a,c,'http://surfobs.climate.copernicus.eu/dataaccess/access_eobs.php' %s",filename))

  filename <- paste0(outputDir_25,"/fg_ens_members_0.25deg_reg_",year,".nc")
  system(sprintf("ncatted -h -O -a history,global,d,, %s",filename))
  system(sprintf("ncatted -O -h -a E-OBS_version,global,a,c,%s %s",version,filename))
  system(sprintf("ncatted -O -h -a References,global,a,c,'http://surfobs.climate.copernicus.eu/dataaccess/access_eobs.php' %s",filename))

  print("year completed!")

}
