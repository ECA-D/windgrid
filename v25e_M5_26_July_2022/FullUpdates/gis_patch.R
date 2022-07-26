#!/usr/bin/Rscript

run.call <- function(ifile, gis){
    require(RMySQL)
    require(RSQLite)
    require(rgeos)
    require(FNN)

    ## coast.file <- file.path("/nobackup/users/cornes/gridding/gis/ne_110m/",
    ##                         "ne_110m_coastline.shp")
    coast.file <- file.path("/data2/Richard/DEM/",
                            "ne_110m_coastline.shp")
    wgs.84 <- "+proj=longlat +datum=WGS84"

    ## Open the output sqlite database
    dbout <- dbConnect(SQLite(), dbname=ifile)
    meta.dat <- dbGetQuery(dbout, "select * from meta")

    ## Add the GIS metadata
    coordinates(meta.dat) <- ~lon+lat
    proj4string(meta.dat) <- CRS(wgs.84)
    meta.dat.proj <- spTransform(meta.dat, CRS(proj4string(gis)))

    ## Find nearest neighbours
    nn <- get.knnx(coordinates(gis.proj), coordinates(meta.dat.proj),1)
    stn.gis <- as.data.frame(gis.proj)[nn$nn.index[,1],]
    colnames(stn.gis)[which(colnames(stn.gis)=="dist2coast")] <- "dist2coast_dem"
    stn.gis <- stn.gis[,!colnames(stn.gis)%in%c("x","y")]

    ## Distance to coast per point
    coasts  <- shapefile(coast.file)
    coasts.proj <- spTransform(coasts, CRS(laa.proj))
    dist2coast <- gDistance(meta.dat.proj,coasts.proj,byid = TRUE)
    dist2coast <- apply(dist2coast,2,min)

    meta.out <- cbind(as.data.frame(meta.dat), stn.gis,
                      dist2coast=dist2coast)

    ## Replace missing elevation values with DEM-derived values
    miss.elev <- which(meta.out$elev<=-999.9)
    meta.out$elev[miss.elev] <- meta.out$alt[miss.elev]
    meta.out$fill.alt <- 0
    meta.out$fill.alt[miss.elev] <- 1   #Flag filled altitude values

    ## Change the alitude names
    colnames(meta.out)[which(colnames(meta.out)=="alt")] <- "alt_dem"
    colnames(meta.out)[which(colnames(meta.out)=="elev")] <- "alt"

    metawrite.out <- dbWriteTable(dbout, "meta", meta.out,
                                  overwrite = TRUE, row.names=FALSE)

    closedb.out <- dbDisconnect(dbout)  
}

library(raster)

setwd("/home/besselaa/Jouke/WindGrid/scripts/Monupdates/")
file.list <- list.files("/data2/Else/Jouke/WindInput/Sqlite/",
                        "eobs_fg_no_homog.sqlite", full.names = TRUE)
gis.file <- file.path("/data2/Richard/eobs/","gtopo30_gis_1km")

laa.proj <- paste("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000",
                  "+ellps=GRS80 +units=m +no_defs", sep=" ")

gis.raster <- brick(gis.file)
gis.proj <- rasterToPoints(gis.raster,spatial=TRUE)
gis.proj <- spTransform(gis.proj, CRS(laa.proj))

## Run the patch
idone <- lapply(file.list, run.call, gis=gis.proj)
