readStationData <- function(year,month)
{
  # local settings
  NstationMax = 5000

  # get number of days in month
  ystr = toString(year)
  if(month<10){
    mstr = paste("0",toString(month),sep="")
  } else {
    mstr = toString(month)
  }
  nDays = monthDays(as.Date(paste(ystr,'-',mstr,'-01',sep='')))

  # read station data
  stationID = matrix(NA,NstationMax,nDays)
  stationLon = matrix(NA,NstationMax,nDays)
  stationLat = matrix(NA,NstationMax,nDays)
  stationSpeed = matrix(NA,NstationMax,nDays)
  NNstation = 0
  for(day in 1:nDays){
    
    if(day<10){
      dstr = paste("0",toString(day),sep="")
    } else {
      dstr = toString(day)
    }

    tempstation0 = read_table(paste("/data2/Else/Jouke/WindInput/Stations/fg_ECAD_",ystr,"/fg_ECAD_",ystr,mstr,dstr,".txt",sep=""),col_names = FALSE,col_types = cols())
    tempstation = list(stationid=tempstation0$X1,longitude=tempstation0$X2,latitude=tempstation0$X3,speed=tempstation0$X4)
    tempstation$speed[tempstation$speed<0] = NA
    Nstation = length(tempstation$stationid)
    NNstation = max(Nstation,NNstation)
    if(Nstation>NstationMax){stop('preallocate more stations in utilities/textload.R -> readStationData()')}
    stationID[1:Nstation,day] = tempstation$stationid
    stationLon[1:Nstation,day] = tempstation$longitude
    stationLat[1:Nstation,day] = tempstation$latitude
    stationSpeed[1:Nstation,day] = tempstation$speed
  }
    
  # tailor station matrices
  id = (rowSums(is.finite(stationSpeed))>0)
  stationID = stationID[id,,drop=FALSE]
  stationLon = stationLon[id,,drop=FALSE]
  stationLat = stationLat[id,,drop=FALSE]
  stationSpeed = stationSpeed[id,,drop=FALSE]
  NNstation = size(stationSpeed)[1]
    
  # store original station data for return
  stationData = list(idRaw = stationID, lonRaw = stationLon, latRaw = stationLat, speedRaw = stationSpeed)

  # find unique station IDs
  uniqueID = unique(matrix(stationID,size(stationID)[1]*size(stationID)[2],1))
  uniqueID = uniqueID[is.finite(uniqueID)]
  nunique = length(uniqueID)
  
  # search for associated unique station coordinates
  uniqueLon = NA*uniqueID
  uniqueLat = NA*uniqueID
  idColumn = matrix(stationID,size(stationID)[1]*size(stationID)[2],1)
  lonColumn = matrix(stationLon,size(stationID)[1]*size(stationID)[2],1)
  latColumn = matrix(stationLat,size(stationID)[1]*size(stationID)[2],1)
  for(sid in 1:nunique){
    for(sidSearch in 1:length(idColumn)){
      if(is.finite(idColumn[sidSearch])){
        if(uniqueID[sid] == idColumn[sidSearch]){
          uniqueLon[sid] = lonColumn[sidSearch]
          uniqueLat[sid] = latColumn[sidSearch]
          break
        }
      }
    }
  }
  
  # search for associated daily averaged speed
  IDfixed = matrix(NA,nunique,nDays)
  Lonfixed = matrix(NA,nunique,nDays)
  Latfixed = matrix(NA,nunique,nDays)
  Speedfixed = matrix(NA,nunique,nDays)
  for(dy in 1:nDays){
    for(sid in 1:nunique){
      IDfixed[sid,dy] = uniqueID[sid]
      Lonfixed[sid,dy] = uniqueLon[sid]
      Latfixed[sid,dy] = uniqueLat[sid]
      for(sidSearch in 1:size(stationID)[1]){
        if(is.finite(stationID[sidSearch,dy])){
          if(IDfixed[sid,dy] == stationID[sidSearch,dy]){
            Speedfixed[sid,dy] = stationSpeed[sidSearch,dy]
          }
        }
      }
    }
  }

  # store final station data for return
  stationData$id = IDfixed
  stationData$lon = Lonfixed
  stationData$lat = Latfixed
  stationData$speed = Speedfixed

  return(stationData)

}

