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
plots = FALSE
windgridInDir = '/data2/Else/Jouke/WindOutput_v2/GPR_gridded/'
windgridOutDir = '/data2/Else/Jouke/WindOutput_v2/Final_NC/'
version = 'v23.1e'
scalefactor = 0.01

####################
## create nc file ##
####################

# get coordinates
grid0 = readRDS(file = paste(windgridInDir,"/windgrid_ensemble_statistics_1980_01.rds",sep=""))
lonvals10km = grid0$lon
latvals10km = grid0$lat
rm(grid0)
lonvals25km = seq(min(lonvals10km-0.05+0.125),max(lonvals10km+0.05-0.125),.25) # to check
latvals25km = seq(min(latvals10km-0.05+0.125),max(latvals10km+0.05-0.125),.25) # to check

# define spatial dims
lon10km <- ncdim_def('longitude','degrees_east', lonvals10km, unlim=FALSE,create_dimvar=TRUE,calendar=NA,longname='Longitude values')
lat10km <- ncdim_def('latitude','degrees_north', latvals10km, unlim=FALSE,create_dimvar=TRUE,calendar=NA,longname='Latitude values')
lon25km <- ncdim_def('longitude','degrees_east', lonvals25km, unlim=FALSE,create_dimvar=TRUE,calendar=NA,longname='Longitude values')
lat25km <- ncdim_def('latitude','degrees_north', latvals25km, unlim=FALSE,create_dimvar=TRUE,calendar=NA,longname='Latitude values')

# define time dim
origin = as.Date('1950-01-01')
day1 = as.Date(paste(toString(min(YEAR)),'-',toString(min(MONTH)),'-01',sep=''))
nDays = monthDays(as.Date(paste(toString(max(YEAR)),'-',toString(max(MONTH)),'-01',sep='')))
day2 = as.Date(paste(toString(max(YEAR)),'-',toString(max(MONTH)),'-',toString(nDays),sep=''))
time1 = difftime(day1,origin, units = "days")
time2 = difftime(day2,origin, units = "days")
time <- ncdim_def('time','days since 1950-01-01 00:00:00',time1:time2, unlim=TRUE,create_dimvar=TRUE,calendar="standard",longname='Time in days')

# define vars
meanSpeed10km <- ncvar_def('fg','m/s',dim=list(lon10km,lat10km,time),missval=-9999.,longname='Ensemble mean wind speed',prec="short",compression=1)
spreadSpeed10km <- ncvar_def('fg','m/s',dim=list(lon10km,lat10km,time),missval=-9999.,longname='Ensemble spread wind speed',prec="short",compression=1)
meanSpeed25km <- ncvar_def('fg','m/s',dim=list(lon25km,lat25km,time),missval=-9999.,longname='Ensemble mean wind speed',prec="short",compression=1)
spreadSpeed25km <- ncvar_def('fg','m/s',dim=list(lon25km,lat25km,time),missval=-9999.,longname='Ensemble spread wind speed',prec="short",compression=1 )

# open files for writing
ncnewMean10km <- nc_create(paste(windgridOutDir,"/fg_ens_mean_0.1deg_reg_",toString(min(YEAR)),"-",toString(max(YEAR)),"_",version,".nc",sep=""),meanSpeed10km)
ncnewSpread10km <- nc_create(paste(windgridOutDir,"/fg_ens_spread_0.1deg_reg_",toString(min(YEAR)),"-",toString(max(YEAR)),"_",version,".nc",sep=""),spreadSpeed10km)
ncnewMean25km <- nc_create(paste(windgridOutDir,"/fg_ens_mean_0.25deg_reg_",toString(min(YEAR)),"-",toString(max(YEAR)),"_",version,".nc",sep=""),meanSpeed25km)
ncnewSpread25km <- nc_create(paste(windgridOutDir,"/fg_ens_spread_0.25deg_reg_",toString(min(YEAR)),"-",toString(max(YEAR)),"_",version,".nc",sep=""),spreadSpeed25km)

# add standard names
ncatt_put(ncnewMean10km,'longitude',"standard_name",'longitude')
ncatt_put(ncnewSpread10km,'longitude',"standard_name",'longitude')
ncatt_put(ncnewMean25km,'longitude',"standard_name",'longitude')
ncatt_put(ncnewSpread25km,'longitude',"standard_name",'longitude')

ncatt_put(ncnewMean10km,'latitude',"standard_name",'latitude')
ncatt_put(ncnewSpread10km,'latitude',"standard_name",'latitude')
ncatt_put(ncnewMean25km,'latitude',"standard_name",'latitude')
ncatt_put(ncnewSpread25km,'latitude',"standard_name",'latitude')

ncatt_put(ncnewMean10km,'longitude',"axis",'X')
ncatt_put(ncnewSpread10km,'longitude',"axis",'X')
ncatt_put(ncnewMean25km,'longitude',"axis",'X')
ncatt_put(ncnewSpread25km,'longitude',"axis",'X')

ncatt_put(ncnewMean10km,'latitude',"axis",'Y')
ncatt_put(ncnewSpread10km,'latitude',"axis",'Y')
ncatt_put(ncnewMean25km,'latitude',"axis",'Y')
ncatt_put(ncnewSpread25km,'latitude',"axis",'Y')

ncatt_put(ncnewMean10km,'fg',"standard_name",'wind_speed')
ncatt_put(ncnewSpread10km,'fg',"standard_name",'wind_speed')
ncatt_put(ncnewMean25km,'fg',"standard_name",'wind_speed')
ncatt_put(ncnewSpread25km,'fg',"standard_name",'wind_speed')

ncatt_put(ncnewMean10km,'time',"standard_name",'time')
ncatt_put(ncnewSpread10km,'time',"standard_name",'time')
ncatt_put(ncnewMean25km,'time',"standard_name",'time')
ncatt_put(ncnewSpread25km,'time',"standard_name",'time')

ncatt_put(ncnewMean10km,'fg',"missing_value",'-9999s')
ncatt_put(ncnewSpread10km,'fg',"missing_value",'-9999s')
ncatt_put(ncnewMean25km,'fg',"missing_value",'-9999s')
ncatt_put(ncnewSpread25km,'fg',"missing_value",'-9999s')

# add attributes
ncatt_put(ncnewMean10km,'fg',"scale_factor",scalefactor,prec="float")
ncatt_put(ncnewMean10km,'fg',"add_offset",0.0,prec="float")
ncatt_put(ncnewSpread10km,'fg',"scale_factor",scalefactor,prec="float")
ncatt_put(ncnewSpread10km,'fg',"add_offset",0.0,prec="float")
ncatt_put(ncnewMean25km,'fg',"scale_factor",scalefactor,prec="float")
ncatt_put(ncnewMean25km,'fg',"add_offset",0.0,prec="float")
ncatt_put(ncnewSpread25km,'fg',"scale_factor",scalefactor,prec="float")
ncatt_put(ncnewSpread25km,'fg',"add_offset",0.0,prec="float")


 

########################
## read & write loops ##
########################
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
      spreadToday10km = t(matrix(gridStatistics$WINDSPEED_ENSEMBLE_RAW_SPREAD_05_95[,day],Nlat,Nlon))
      
      # filter grid (outlier detection)
      idout = (meanToday10km<=0)
      meanToday10km[idout] = NA
      spreadToday10km[idout] = NA

      idout = (meanToday10km>35)
      meanToday10km[idout] = NA
      spreadToday10km[idout] = NA
      
      idout = (spreadToday10km<=0)
      meanToday10km[idout] = NA
      spreadToday10km[idout] = NA
      
      idout = (spreadToday10km>35)
      meanToday10km[idout] = NA
      spreadToday10km[idout] = NA
      
      # interpolate to 25 km
      COOR = meshgrid(lonvals25km,latvals25km)
      LON = COOR$X
      LAT = COOR$Y
      meanToday25km = t(matrix(interp2(latvals10km,lonvals10km,meanToday10km,LAT,LON),length(latvals25km),length(lonvals25km)))
      spreadToday25km = t(matrix(interp2(latvals10km,lonvals10km,spreadToday10km,LAT,LON),length(latvals25km),length(lonvals25km)))
      
      # append to 10 km netcdf files
      ncvar_put(ncnewMean10km,meanSpeed10km,meanToday10km/scalefactor,start=c(1,1,timeCounter),count=c(length(lonvals10km),length(latvals10km),1))
      ncvar_put(ncnewSpread10km,spreadSpeed10km,spreadToday10km/scalefactor,start=c(1,1,timeCounter),count=c(length(lonvals10km),length(latvals10km),1))
      
      # append to 25 km netcdf files
      ncvar_put(ncnewMean25km,meanSpeed25km,meanToday25km/scalefactor,start=c(1,1,timeCounter),count=c(length(lonvals25km),length(latvals25km),1))
      ncvar_put(ncnewSpread25km,spreadSpeed25km,spreadToday25km/scalefactor,start=c(1,1,timeCounter),count=c(length(lonvals25km),length(latvals25km),1))
      
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
}

# close nc files
nc_close(ncnewMean10km)
nc_close(ncnewSpread10km)
nc_close(ncnewMean25km)
nc_close(ncnewSpread25km)

# add global attributes
message("adding global attributes")

filename <- paste0(windgridOutDir,"/fg_ens_mean_0.1deg_reg_",toString(min(YEAR)),"-",toString(max(YEAR)),"_",version,".nc")
system(sprintf("ncatted -h -O -a history,global,d,, %s",filename))
system(sprintf("ncatted -O -h -a E-OBS_version,global,a,c,%s %s",version,filename))
system(sprintf("ncatted -O -h -a References,global,a,c,'http://surfobs.climate.copernicus.eu/dataaccess/access_eobs.php' %s",filename))

filename <- paste0(windgridOutDir,"/fg_ens_mean_0.25deg_reg_",toString(min(YEAR)),"-",toString(max(YEAR)),"_",version,".nc")
system(sprintf("ncatted -h -O -a history,global,d,, %s",filename))
system(sprintf("ncatted -O -h -a E-OBS_version,global,a,c,%s %s",version,filename))
system(sprintf("ncatted -O -h -a References,global,a,c,'http://surfobs.climate.copernicus.eu/dataaccess/access_eobs.php' %s",filename))

filename <- paste0(windgridOutDir,"/fg_ens_spread_0.1deg_reg_",toString(min(YEAR)),"-",toString(max(YEAR)),"_",version,".nc")
system(sprintf("ncatted -h -O -a history,global,d,, %s",filename))
system(sprintf("ncatted -O -h -a E-OBS_version,global,a,c,%s %s",version,filename))
system(sprintf("ncatted -O -h -a References,global,a,c,'http://surfobs.climate.copernicus.eu/dataaccess/access_eobs.php' %s",filename))

filename <- paste0(windgridOutDir,"/fg_ens_spread_0.25deg_reg_",toString(min(YEAR)),"-",toString(max(YEAR)),"_",version,".nc")
system(sprintf("ncatted -h -O -a history,global,d,, %s",filename))
system(sprintf("ncatted -O -h -a E-OBS_version,global,a,c,%s %s",version,filename))
system(sprintf("ncatted -O -h -a References,global,a,c,'http://surfobs.climate.copernicus.eu/dataaccess/access_eobs.php' %s",filename))
