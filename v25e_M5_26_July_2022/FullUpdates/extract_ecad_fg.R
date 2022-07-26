#!/usr/bin/Rscript

library(eobs)
library(parallel)
library(methods)
suppressWarnings(library(randomFunctions)) #Prevent trivial warnings

fetch.data.def <- function(obj,ivar,dbin,homog,stn,synop){
    require(RMySQL)
    require(data.table)

    if(synop){
        datatable <- sprintf("synops_%s",ivar)
        statement <- paste0("SELECT DATE_FORMAT(syn_date,'%%Y') as year, ",
                            "DATE_FORMAT(syn_date,'%%m') as month, ",
                            "DATE_FORMAT(syn_date,'%%d') as day, ",
                            "%s as var, qc, syn_id as id, syn_id as blend_id ",
                            "FROM %s ",
                            "WHERE syn_id=%d")
        statement <- sprintf(statement, ivar, datatable, stn)
    } else {
        itable <- ifelse(homog,"homblend","blended")
        datatable <- sprintf("series_%s_%s_mixed",ivar,itable)
        statement <- paste0("SELECT DATE_FORMAT(ser_date,'%%Y') as year, ",
                            "DATE_FORMAT(ser_date,'%%m') as month, ",
                            "DATE_FORMAT(ser_date,'%%d') as day, ",
                            "%s as var, qc, %s.sta_id as id, blend_id ",
                            "FROM %s ",
                            "WHERE %s.sta_id=%d ")
        statement <- sprintf(statement, ivar, datatable, datatable, datatable, stn)
    }
    
    rs <- dbSendQuery(dbin, statement) 
    dat <- fetch(rs, n=-1)
    dat <- as.data.table(dat)          #Construct a data.table from dat

    ## QC
    dat[qc>0,var:=NA]

    ## Synop codes are in a different part of the database.
    ## Hence, we need to split this ele_name construction
    ids <- unique(dat$blend_id)
    ids.synop <- ids[ids>=900000]
    ids.normal <- ids[ids<900000]

    if(length(ids.normal)>0){
        statement.ids <- paste0("SELECT series.ser_id as blend_id,ele_name ",
                                "FROM series,elements ",
                                "WHERE series.ele_id=elements.ele_id ",
                                "AND series.ser_id in (%s)")
        ids.normal <- paste(ids.normal,collapse = ",")
        statement.ids <- sprintf(statement.ids,ids.normal)
        rs <- dbSendQuery(dbin, statement.ids) 
        ele.names <- fetch(rs, n=-1)
    }

    if(length(ids.synop)>0){
        statement.ids <- paste0("SELECT syn_id as blend_id,ele_name ",
                                "FROM synops_derived WHERE ele_name like '%%%s%%' ",
                                "AND syn_id in (%s)")
        ids.synop <- paste(ids.synop,collapse = ",")
        statement.ids <- sprintf(statement.ids, ivar, ids.synop)
        rs <- dbSendQuery(dbin, statement.ids) 
        ele.synop <- fetch(rs, n=-1)
        ele.names <- tryCatch(rbind(ele.names,ele.synop), error=function(e){ele.synop})
    }
    
    ele.names <- data.table(ele.names)
    setnames(ele.names, c("blend_id","ele_name"))
    dat <- merge(dat, ele.names, by="blend_id", all.x = TRUE, all.y = FALSE)
    
    ## Prepare the data
    dat[var==-9999, var:=NA]           #Change missing values to NA
    if(!synop) dat[qc>0,var:=NA]                  #Set suspects to NA if required

    # wind speed: m/s
    dat[,var:=var/10.]                 #Convert to decimal units
    date.vars <- c("year","month","day")
    dat[,(date.vars):=lapply(.SD,as.numeric), .SDcols=date.vars]
    
    if(nrow(dat)==0) stop() # If there are no data raise an error

    list("dat"=dat)
}

fetch.data.dtr <- function(obj,ivar,dbin,homog,stn,synop){
    dat.tx <- fetch.data.def(obj,"tx",dbin,homog,stn,synop)
    dat.tn <- fetch.data.def(obj,"tn",dbin,homog,stn,synop)
    list("tmax"=dat.tx$dat,"tmin"=dat.tn$dat)
}

setGeneric("fetch.data", function(obj, ...) { standardGeneric("fetch.data") })
setMethod("fetch.data", "eobs", fetch.data.def)
setMethod("fetch.data", "eobs_rain", fetch.data.def)
setMethod("fetch.data", "eobs_dtr", fetch.data.dtr)
setMethod("fetch.data", "eobs_tmn", fetch.data.dtr)

fetch.stns.def <- function(obj,ivar, dbin, homog, synop){
    require(RMySQL)

    if(!synop){
        metatable <- "series_blended_mixed_derived"
        stn.statement <- paste0("SELECT sta_id as id ",
                                "FROM %s ",
                                "WHERE %s.syn_kind='%s'")
        stn.statement <- sprintf(stn.statement, metatable, metatable, ivar)
    } else {
        stn.statement <- "SELECT syn_id as id FROM wmostations"
    }
    
    stn.rs <- dbSendQuery(dbin, stn.statement) 
    stns <- fetch(stn.rs, n=-1)

    stns
}

fetch.stns.dtr <- function(obj,ivar,dbin,homog){
    stns.tx <- fetch.stns.def(obj,"tx", dbin, homog)
    stns.tn <- fetch.stns.def(obj,"tn", dbin, homog)
    stns <- rbind(stns.tx,stns.tn)
    icnt <- as.data.frame(table(stns$id))
    id <- as.numeric(icnt$Var1[icnt$Freq==2])
    stns <- data.frame(id=id)
    stns
}

setGeneric("fetch.stns", function(obj,...) standardGeneric("fetch.stns"))
setMethod("fetch.stns", "eobs", fetch.stns.def)
setMethod("fetch.stns", "eobs_rain", fetch.stns.def)
setMethod("fetch.stns", "eobs_dtr", fetch.stns.dtr)
setMethod("fetch.stns", "eobs_tmn", fetch.stns.dtr)

run.call <- function(obj, homog=FALSE, mon.update=FALSE, synop=FALSE){
    require(RMySQL)
    require(RSQLite)
    require(data.table)
    require(eobs)
    require(lubridate)
    
    ivar <- obj@var
    wgs.84 <- "+proj=longlat +datum=WGS84"
    laa.proj <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"
    filetag <- ifelse(homog,"homog","no_homog")
    
    if(mon.update){
        ## Select the data for the current year Up to end of last
        ## month. Unless we're in January, in which case select all of
        ## last year
        if(month(Sys.Date())==1){
            start <- make_datetime(year(Sys.Date())-1)
            end <- make_datetime(year(Sys.Date())-1,12,31)
        } else {
            end <- as.Date(cut(Sys.Date(), "month"))
            start <- as.Date(cut(Sys.Date(), "year"))
            end <- seq(end,length=2, by="-1 day")[2]
        }
        ofile <-  sprintf("eobs_%s_%s_monupdate.sqlite", ivar,filetag)
    } else {
        ofile <- sprintf("eobs_%s_%s.sqlite", ivar,filetag)
    }
    
    if(synop) ofile <- paste0("synop_", ofile)

    ## The output sqlite database (delete previous version)    
    if(file.exists(ofile)){
        warning(paste("Removing file", ofile))
        iremove <- file.remove(ofile)
        stopifnot(iremove)
    }

    ## Open the ouput database
    dbout <- dbConnect(SQLite(), dbname=ofile)

    ## Open the input MySQL database
#    dbin <- dbConnect(MySQL(), user='rccuser', password='rcc123', dbname='ECAD')
    dbin <- dbConnect(MySQL(), host='bhle4m.knmi.nl', user='ecadro', password='R0nly_Ekat', dbname='ECAD')

    ## Read stations
    stns <- fetch.stns(obj,ivar,dbin,homog,synop)

    ## Loop through each station
    for(stn in stns$id){
    
        FD <- tryCatch({
            fetch.data(obj, ivar, dbin, homog, stn, synop)
        }, error = function(e) {
            message("Missing data. Skipping...")
            NULL
        })
        
        if(is.null(FD)) next            #Iterate if error is raised

        prep.args <- list("obj"=obj, id=stn, "baseThresh"=0.1, "gammaFunc"="GammaFunc.SCI")
        prep.args <- c(prep.args, FD)
        DT <- do.call("prep_data", prep.args)

        if(any(c(nrow(DT$day),nrow(DT$month))==0)) next #Iterate if no data

        if(mon.update){
            day.dat <- DT$day
            day.dat[,date:=as.Date(paste(year,month,day,sep="-"),tz="GMT")]
            day.dat <- day.dat[date>=start&date<=end]
            day.dat[,date:=NULL]
            
            mon.dat <- DT$month
            mon.dat[,date:=as.Date(paste(year,period,1,sep="-"),tz="GMT")]
            mon.dat <- mon.dat[date>=start&date<=end]
            mon.dat[,date:=NULL]

            DT$day <- day.dat
            DT$month <- mon.dat
        }
        
        iwrite <- sapply(names(DT), function(x){
            dbWriteTable(dbout, x, DT[[x]], append=TRUE, row.names = FALSE)
        })
        
    }

    ## Read the meta-data
    if(synop){
        ## Get the variable id from data in the dbout database
        var.stns <- dbGetQuery(dbout,"SELECT id from day UNION SELECT id from month UNION SELECT id from climatol")
        var.stns <- paste(var.stns$id,collapse=",")
        
        meta.statement <- paste0("SELECT wmostations.syn_id AS id, ",
                                 "wmostations.name AS station, country.name AS country, ",
                                 "wmostations.lat/3600 AS lat, wmostations.lon/3600 AS lon, ",
                                 "elev/10.0 AS elev ",
                                 "FROM wmostations,country ",
                                 "WHERE wmostations.syn_id in (%s) ",
                                 "AND wmostations.coun_id=country.coun_id")
        meta.statement <- sprintf(meta.statement,var.stns)
        
    } else {
        meta.statement <- paste0("SELECT %s.sta_id AS id, ",
                                 "stations.name AS station, country.name AS country, ",
                                 "stations.lat/3600 AS lat, stations.lon/3600 AS lon, ",
                                 "elev/10.0 AS elev ",
                                 "FROM %s,stations,country ",
                                 "WHERE %s.syn_kind='%s' ",
                                 "AND %s.sta_id=stations.sta_id ",
                                 "AND stations.coun_id=country.coun_id")
        
        meta.var <- ifelse(ivar%in%c("dtr","tmn"),"tx",ivar) #Just use tx for dtr/tmn metadata
        metatable <- "series_blended_mixed_derived"
        meta.statement <- sprintf(meta.statement, metatable,metatable,metatable,
                                  meta.var, metatable)
    }
    
    meta.rs <- dbSendQuery(dbin, meta.statement) 
    meta.dat <- fetch(meta.rs, n=-1)
    metawrite.out <- dbWriteTable(dbout, "meta", meta.dat, row.names=FALSE)

    ## Create indices for quicker time slicing
    dbSendQuery(dbout,"CREATE INDEX day_idx ON day(year,period)")
    dbSendQuery(dbout,"CREATE INDEX mon_idx ON month(year,period)")
    dbSendQuery(dbout,"CREATE INDEX climate_idx ON climatol(period)")
    
    closedb.out <- dbDisconnect(dbout)  
    closedb.in <- dbDisconnect(dbin)  
}

setwd("/data2/Else/Jouke/WindInput/Sqlite/")

args <- getArgs()
mon.update <- ifelse(is.null(args$mon), FALSE, TRUE)
synop <- ifelse(is.null(args$synop), FALSE, TRUE)
homog <- ifelse(is.null(args$homog), FALSE, TRUE)

cl <- makeCluster(6)
clusterExport(cl, list("run.call","fetch.data","fetch.stns","fetch.data.def","fetch.stns.def","mon.update","homog","synop"))

## ivars <- c("tx","tn","tg","rr","dtr","tmn","pp")
#ivars <- c("tx","tn","tg","rr","pp")
ivars <- c("fg")

## The years required can be specified below. This is set to a large
## future range but only data that are actually present are
## returned. In the case of monthly updates this is also cut to the
## relevant month.
names(ivars) <- ivars
obj.list <- lapply(ivars, function(v){
    obj <- startEOBS(c(0,0),c(0,0),c(1961,1990),c(1900,2099),1:12, var=v)
})

idone <- parLapply(cl, obj.list, function(x) run.call(x, homog=homog, mon.update=mon.update, synop=synop))

stopCluster(cl)
