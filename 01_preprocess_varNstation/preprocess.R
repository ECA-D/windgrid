#rm()
#dev.off()

###############
## libraries ##
###############
library(pracma)
library(sp)
#library(ggplot2)
library(ncdf4)
library(raster)
library(readr)
require(Hmisc)
library(matrixStats)
#library(plot3D)
#library(MASS)
#library(lattice)

setwd('~/Jouke/WindGrid/scripts/01_preprocess_varNstation/')
source('../91_utilities/textload.R')

##############
## settings ##
##############
YEAR = 2021
MNTH = 1:5
lonmin = -11.05 # -11.05
lonmax = 39.95 # 39.95
latmin = 35.05 # 35.05
latmax = 69.95 # 69.95
resolution = 0.1 # .1
NstationMax = 5000 # 5000
nanValue = -999.9 # currently checking: <0
blur = TRUE # TRUE
minRoughnessYear = 1980
maxRoughnessYear = 2018
verbose = 1
degTreshold = 0.001

outputDir = '/data2/Else/Jouke/WindOutput_v2/Preprocess/'

#########
## DEM ##
#########
gridname = '/data2/Richard/eobs/gtopo30_gis_1km.gri'
rawDEM <- stack(gridname)

# extract DEM for grid
grid = list(lon=t(seq(lonmin,lonmax,resolution)),lat=t(seq(latmin,latmax,resolution)))
Nlon = length(grid$lon)
Nlat = length(grid$lat)
coord = meshgrid(grid$lon,grid$lat)
grid$LON = matrix(coord$X,Nlon*Nlat,1)
grid$LAT = matrix(coord$Y,Nlon*Nlat,1)
if(blur){
  grid$ALT = matrix(extract(rawDEM$alt_blur,cbind(grid$LON,grid$LAT)),Nlon*Nlat,1)
  grid$SLOPE = matrix(extract(rawDEM$slope_blur,cbind(grid$LON,grid$LAT)),Nlon*Nlat,1)
  grid$TPI = matrix(extract(rawDEM$tpi_blur,cbind(grid$LON,grid$LAT)),Nlon*Nlat,1)
  grid$D2C = matrix(extract(rawDEM$dist2coast_blur,cbind(grid$LON,grid$LAT)),Nlon*Nlat,1)
} else {
  grid$ALT = matrix(extract(rawDEM$alt,cbind(grid$LON,grid$LAT)),Nlon*Nlat,1)
  grid$SLOPE = matrix(extract(rawDEM$slope,cbind(grid$LON,grid$LAT)),Nlon*Nlat,1)
  grid$TPI = matrix(extract(rawDEM$tpi,cbind(grid$LON,grid$LAT)),Nlon*Nlat,1)
  grid$D2C = matrix(extract(rawDEM$dist2coast,cbind(grid$LON,grid$LAT)),Nlon*Nlat,1)
}

##########################
## ERA5 monthly average ##
##########################
# read ERA5 file
era5name = '/data2/Else/Jouke/WindInput/ERA5Wind/WDSP_monmean_800hPa_1980-2020_EOBSgrid.nc'
ncin <- nc_open(era5name)
era5lon = ncvar_get(ncin,"lon") # list in documentation
era5lat = ncvar_get(ncin,"lat") # list in documentation
era5speed = ncvar_get(ncin,"u") # sqrt(u^2+v^2) list in documentation

#####################
## Loop over years ##
#####################
for(iterYear in 1:length(YEAR)){
  year = YEAR[iterYear]
  
  # read roughness file
  ncfroughnessname = paste('/data2/Else/Jouke/WindInput/Roughness/',toString(max(min(year,maxRoughnessYear),minRoughnessYear)),'.nc',sep='')
  ncin <- nc_open(ncfroughnessname)
  roughlon = ncvar_get(ncin,"lon") # document
  roughlat = ncvar_get(ncin,"lat") # document
  roughsr = ncvar_get(ncin,"sr") # document
  nRoughDays = size(roughsr)[3]

  # initialize data for this year
  yearDay = 1
  Nmonth = length(MNTH)
  station = list(lon=matrix(NA,NstationMax,Nmonth),lat=matrix(NA,NstationMax,Nmonth))
  station$stationid = matrix(NA,NstationMax,Nmonth)
  station$alt = matrix(NA,NstationMax,Nmonth)
  station$slope = matrix(NA,NstationMax,Nmonth)
  station$tpi = matrix(NA,NstationMax,Nmonth)
  station$d2c = matrix(NA,NstationMax,Nmonth)
  station$avg_windspeed = matrix(NA,NstationMax,Nmonth)
  station$n_in_avg = matrix(NA,NstationMax,Nmonth)
  station$avg_roughness = matrix(NA,NstationMax,Nmonth)
  station$era5_avg_windspeed = matrix(NA,NstationMax,Nmonth)
  grid$ROUGH = matrix(NA,Nlon*Nlat,Nmonth)
  grid$ERA5SPEED = matrix(NA,Nlon*Nlat,Nmonth)
  
  ######################
  ## Loop over months ##
  ######################
  for(monthIter in 1:Nmonth){
    month = MNTH[monthIter]
    print(paste('Year:',toString(year),'-> Month:',toString(MNTH[monthIter]),'...'))
    ystr = toString(year)
    if(MNTH[monthIter]<10){
      mstr = paste("0",toString(MNTH[monthIter]),sep="")
    } else {
      mstr = toString(MNTH[monthIter])
    }

    # get number of days in month
    nDays = monthDays(as.Date(paste(ystr,'-',mstr,'-01',sep='')))

    # read station data
    stationData = readStationData(year,month)
    Nstation = size(stationData$speed)[1]
    station$stationid[1:Nstation,monthIter] = stationData$id[,1]
    station$lon[1:Nstation,monthIter] = stationData$lon[,1]
    station$lat[1:Nstation,monthIter] = stationData$lat[,1]

    ## crop stations
    idcrop = (station$lon[1:Nstation,monthIter]<lonmin) |
             (station$lon[1:Nstation,monthIter]>lonmax) |
             (station$lat[1:Nstation,monthIter]<latmin) |
             (station$lat[1:Nstation,monthIter]>latmax)
    stationData$speed[idcrop,] = NA
    
    # extract DEM for stations
    if(blur){
      station$alt[1:Nstation,monthIter] = matrix(extract(rawDEM$alt_blur,cbind(station$lon[1:Nstation,monthIter,drop=FALSE],station$lat[1:Nstation,monthIter,drop=FALSE])),Nstation,1)
      station$slope[1:Nstation,monthIter] = matrix(extract(rawDEM$slope_blur,cbind(station$lon[1:Nstation,monthIter,drop=FALSE],station$lat[1:Nstation,monthIter,drop=FALSE])),Nstation,1)
      station$tpi[1:Nstation,monthIter] = matrix(extract(rawDEM$tpi_blur,cbind(station$lon[1:Nstation,monthIter,drop=FALSE],station$lat[1:Nstation,monthIter,drop=FALSE])),Nstation,1)
      station$d2c[1:Nstation,monthIter] = matrix(extract(rawDEM$dist2coast_blur,cbind(station$lon[1:Nstation,monthIter,drop=FALSE],station$lat[1:Nstation,monthIter,drop=FALSE])),Nstation,1)
    } else {
      station$alt[1:Nstation,monthIter] = matrix(extract(rawDEM$alt,cbind(station$lon[1:Nstation,monthIter,drop=FALSE],station$lat[1:Nstation,monthIter,drop=FALSE])),Nstation,1)
      station$slope[1:Nstation,monthIter] = matrix(extract(rawDEM$slope,cbind(station$lon[1:Nstation,monthIter,drop=FALSE],station$lat[1:Nstation,monthIter,drop=FALSE])),Nstation,1)
      station$tpi[1:Nstation,monthIter] = matrix(extract(rawDEM$tpi,cbind(station$lon[1:Nstation,monthIter,drop=FALSE],station$lat[1:Nstation,monthIter,drop=FALSE])),Nstation,1)
      station$d2c[1:Nstation,monthIter] = matrix(extract(rawDEM$dist2coast,cbind(station$lon[1:Nstation,monthIter,drop=FALSE],station$lat[1:Nstation,monthIter,drop=FALSE])),Nstation,1)
    }
    
    # replace nan by 0 in DEM for stations
    id = is.na(station$alt[1:Nstation,monthIter])
    station$alt[id,monthIter] = 0
    id = is.na(station$slope[1:Nstation,monthIter])
    station$slope[id,monthIter] = 0
    id = is.na(station$tpi[1:Nstation,monthIter])
    station$tpi[id,monthIter] = 0
    id = is.na(station$d2c[1:Nstation,monthIter])
    station$d2c[id,monthIter] = 0
    
    # era5 windspeed
    station$era5_avg_windspeed[1:Nstation,monthIter] = interp2(era5lat,era5lon,era5speed[,,monthIter],station$lat[1:Nstation,monthIter],station$lon[1:Nstation,monthIter])
    grid$ERA5SPEED[,monthIter] = interp2(era5lat,era5lon,era5speed[,,monthIter],grid$LAT,grid$LON)
    
    ###################################################
    # get daily wind speed and roughness for stations #
    ###################################################
    todaySpeed = matrix(NA,Nstation,nDays)
    todayStationRough = matrix(NA,Nstation,nDays)
    todayGridRough = matrix(NA,Nlon*Nlat,nDays)
    for(dayIter in 1:nDays){
      if(dayIter<10){
        dstr = paste("0",toString(dayIter),sep="")
      } else {
        dstr = toString(dayIter)
      }
      #print(paste(ystr,mstr,dstr))
      
      todaySpeed[1:Nstation,dayIter] = stationData$speed[1:Nstation,dayIter]
      
      # roughness
      todayStationRough[,dayIter] = interp2(roughlat,roughlon,roughsr[,,yearDay],station$lat[1:Nstation,monthIter],station$lon[1:Nstation,monthIter])
      todayGridRough[,dayIter] = interp2(roughlat,roughlon,roughsr[,,yearDay],grid$LAT,grid$LON)
      if(yearDay<nRoughDays){
        yearDay = yearDay + 1
      }
        
    }
    
    # compute monthly-average windspeed and roughness
    station$avg_windspeed[1:Nstation,monthIter] = rowMeans(todaySpeed,na.rm=TRUE)
    station$n_in_avg[1:Nstation,monthIter] = rowSums(is.finite(todaySpeed))
    station$avg_roughness[1:Nstation,monthIter] = rowMeans(todayStationRough,na.rm=TRUE)
    grid$ROUGH[,monthIter] = rowMeans(todayGridRough,na.rm=TRUE)
    
    # remove average if n_in_avg < 70%
    idout = (station$n_in_avg[,monthIter]/nDays)<0.70
    station$avg_windspeed[idout,monthIter] = NA
    station$n_in_avg[idout,monthIter] = 0
    
  }
  
  ## write data for current year to disk
  saveRDS(grid, file =paste(outputDir,"grid_",toString(year),".rds",sep=""))
  saveRDS(station, file =paste(outputDir,"station_",toString(year),".rds",sep=""))
  
}

