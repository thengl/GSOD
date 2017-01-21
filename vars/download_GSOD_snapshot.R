## Import, formating and binding of the Daily Global Weather Measurements, 1929-2009 (NCDC, GSOD) yrs = 1901:2016
## Description of files available at [http://www7.ncdc.noaa.gov/CDO/GSOD_DESC.txt]
# by: Tom.Hengl@isric.org

library(RCurl)
library(rgdal)
library(readr)
library(parallel)
library(devtools)
#devtools::install_github('adamhsparks/GSODR')
library("GSODR")
library(snowfall)
library(data.table)
library(spacetime)
library(plotKML)
library(snowfall)
library(raster)
library(tools)

## Save objects in parallel:
source("/data/models/saveRDS_functions.R")

## Download all tar files from server:
setwd("/data/GSOD")
ftp.GSOD <- "ftp://ftp.ncdc.noaa.gov/pub/data/gsod/"
system(paste("wget -A tar -m -p -E -k -K -np", ftp.GSOD))
## unzip files:
tar.lst <- list.files(pattern=glob2rx("*.tar"), full.names = TRUE, recursive = TRUE)
## 116
## TAKES ABOUT 8 hrs TO EXTRACT IN PARALLEL
x = parallel::mclapply(tar.lst, FUN=function(x){system(paste("7za e -ttar", x, "-y"))}, mc.cores=48)
GSOD.list <- list.files(pattern=glob2rx("*.gz"), full.names = TRUE, path = "/data/GSOD")
## 436,576 files! 3.1GB of data
setwd("./tmp")
saveRDS.gz(GSOD.list, "GSOD.list.rds")

## we use only few functions from the GSODR package:
source_https <- function(url, ...) {
  require(RCurl)
  cat(getURL(url, followlocation = TRUE, cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl")), file = basename(url))
  source(basename(url))
}
source_https("https://raw.githubusercontent.com/adamhsparks/GSODR/master/R/get_GSOD.R")

## Get station data:
stations = .fetch_stations()
stations$STNID <- as.factor(paste(stations$USAF, stations$WBAN, sep="-"))  # WMO/DATSAV3 number
stations.sp = as.data.frame(stations)[,c("STNID","LON","LAT")]
stations.sp = stations.sp[!is.na(stations.sp$LAT)&!is.na(stations.sp$LON),]
coordinates(stations.sp) <- c("LON", "LAT")
proj4string(stations.sp) <- CRS("+proj=longlat +datum=WGS84") 

Read_gz <- function(x){
  out = data.table(.read_gz(x))
  out$STNID = paste(out$STN, out$WBAN, sep="-")
  #out$DATE = as.Date(paste(out$YEAR, substr(out$MODA, 1, 2), substr(out$MODA, 3, 4), sep="-"))
  return(out)
}

v1 = .read_gz(GSOD.list[5])
## Read all data in parallel (TAKES ca 1 HR):
sfInit(parallel=TRUE, cpus=48)
sfExport("GSOD.list", "stations",".read_gz")
sfLibrary(readr)
sfLibrary(data.table)
GSOD.tbl <- sfClusterApplyLB(GSOD.list, function(x){try( Read_gz )})
sfStop()
## 31GB
sel.r = unlist(parallel::mclapply(GSOD.tbl, function(x){all(names(x) %in% names(v1))}, mc.cores = 48))
summary(sel.r)
## Efficient bind:
#system.time( GSOD.tbl <- rbindlist(GSOD.tbl[sel.r]) )
system.time( GSOD.tbl <- do.call(rbind, GSOD.tbl[sel.r]) )
## Save 32GB object in parallel:
saveRDS.gz(GSOD.tbl, paste0("GSOD_snapshot_", format(Sys.Date(), "%Y_%b_%d"), ".rds"))
dim(GSOD.tbl)
## 124,709,417 records 

## Subset based on variable of interest
#GSOD.tbl$STNID <- as.factor(paste(GSOD.tbl$STN, GSOD.tbl$WBAN, sep="-")) 
t.vars = c("TEMP","DEWP","SLP","STP","VISIB","WDSP","MXSPD","MAX","MIN","PRCP","SNDP","I_")
for(j in 1:length(t.vars)){
  out.rds = paste0("GSOD_snapshot_", t.vars[j], "_", format(Sys.Date(), "%Y_%b_%d"), ".rds")
  if(!file.exists(out.rds)){
    cnames = c("STNID", "YEAR", "MODA", names(GSOD.tbl)[grep(pattern=t.vars[j], names(GSOD.tbl))])
    values = subset(GSOD.tbl, select=cnames)
    values = subset(values, complete.cases(values))
    x = list(stations=stations, values=values)
    saveRDS.gz(x, out.rds)
    rm(x)
    gc()
  }
}

## examples with temperatures:
tmp.f = readRDS.gz(list.files(pattern="TEMP")) ## 4GB
tmp.f$values = subset(tmp.f$values, tmp.f$values$YEAR=="2015")
tmp.f$values$DATE <- as.Date(paste(tmp.f$values$YEAR, substr(tmp.f$values$MODA, 1, 2), substr(tmp.f$values$MODA, 3, 4), sep="-"))
tmp.f$values = subset(tmp.f$values, tmp.f$values$DATE>"2015-05-01" & tmp.f$values$DATE<"2015-05-10")

summary(tmp.f$values$TEMP)
## convert to Celsius
tmp.f$values$TEMPC <- round((tmp.f$values$TEMP-32)*5/9, 1)
rn.t = quantile(tmp.f$values$TEMPC, c(0.1,0.9))
tmp.sp = plyr::join(tmp.f$stations, tmp.f$values)
str(tmp.sp)
tmp.sp = tmp.sp[!is.na(tmp.sp$TEMP),]
coordinates(tmp.sp) <- c("LON", "LAT")
proj4string(tmp.sp) <- CRS("+proj=longlat +datum=WGS84") 

tmp_ST <- STIDF(sp=as(tmp.sp, "SpatialPoints"), time=tmp.sp$DATE-0.5, data=tmp.sp@data[,c("TEMPC","STNID")], endTime=tmp.sp$DATE+0.5)
shape = "http://maps.google.com/mapfiles/kml/pal2/icon18.png"
kml(tmp_ST, dtime = 24*3600, colour = TEMPC, shape = shape, labels = TEMPC, file.name="Temperatures_daily_2015.kml", folder.name="TEMP", z.lim=rn.t)
system("zip -m Temperatures_daily_2015.kmz Temperatures_daily_2015.kml")

## Overlay GSOD stations and CHELSA climate / MODIS snow probs / EarthEnv climatic layers
clim.tifs = c(list.files(path = "/data/Climate", pattern = glob2rx("*.tif$"), full.names = TRUE), list.files(path = "/data/EarthEnv", pattern = glob2rx("MODCF_*.tif$"), full.names = TRUE), list.files(path = "/data/ESA_global", pattern = glob2rx("ESACCI_snow_prob_*.tif$"), full.names = TRUE))
## 71 in total
## overlay function
extract.tif = function(x, y){
  r = raster(x)
  if(!is.na(proj4string(r))){
    y = spTransform(y, proj4string(r))
  }  
  out = raster::extract(y=y, x=r)
  return(out)
}
## overlay in parallel TAKES ca 10 MINS:
sfInit(parallel=TRUE, cpus=48)
sfExport("stations.sp", "clim.tifs", "extract.tif")
sfLibrary(rgdal)
sfLibrary(raster)
clim.out <- data.frame(sfClusterApplyLB(clim.tifs, function(i){try( extract.tif(i, stations.sp) )}))
sfStop()
names(clim.out) = file_path_sans_ext(basename(clim.tifs))
gsod.clim = cbind(as.data.frame(stations.sp), clim.out)
str(gsod.clim)
sel.rm = rowSums(gsod.clim[,names(clim.out)])
GSOD_clim = gsod.clim[!is.na(sel.rm),]
#names(GSOD_clim)[which(names(GSOD_clim)=="CHELSA_temp_5_1979-2013_V1_1")] = "CHELSA_temp_5_1979-2013"
## 23,927 stations
save(GSOD_clim, file="GSOD_clim.rda", compress="xz")
write.csv(data.frame(Names=names(GSOD_clim)), "GSOD_clim.names.csv")

## ----------------------------------


fmd0$Day = as.Date(as.POSIXct(fmd0$data_rapporto, format="%d/%m/%Y"))
fmd0$sex = as.factor(fmd0$sesso)
fmd0$age = as.numeric(gsub(",", "\\.", fmd0$eta))

proj4string(fmd0) <- CRS("+proj=longlat +datum=WGS84")
writeOGR(fmd0, "emergency_calls_pnts.shp", "emergency_calls_pnts", "ESRI Shapefile")

tmp.f$PREC <- round(ifelse(tmp.f$PREC==99.9, NA, tmp.f$PREC*25.4), 2)  # convert to mm
## plot in Google Earth:
coordinates(GSOD.2008.XY) <- ~LON+LAT
proj4string(GSOD.2008.XY) <- CRS("+proj=longlat +datum=WGS84")

GSOD.20080501.XY <- subset(GSOD.2008.XY, GSOD.2008.XY$DATE==as.Date("2008-05-01"))
# bubble(GSOD.20080501.XY, "TEMPC")
writeOGR(GSOD.20080501.XY[,c("STNID", "TEMPC", "PREC")], "GSOD20080501_XY.kml", "GSOD20080501_XY", "KML")
writeOGR(GSOD.20080501.XY[,c("STNID", "YEARMODA", "TEMPC", "TEMP.count", "PREC", "PREC.flag")], "GSOD20080501_XY.shp", "GSOD20080501_XY", "ESRI Shapefile")



system(paste("7za e -tgzip", i, " *.op -aos"), show.output.on.console=FALSE)
tmp <- readLines(paste(in.dir, "/", strsplit(i, ".gz")[[1]][1], sep=""))
tmp <- gsub(pattern='[[:space:]]+', replacement=',', tmp)
tmp.txt <- tempfile(,fileext = ".txt")
write.table(t(matrix(unlist(strsplit(tmp, ",")[-1]), nrow=22)), tmp.txt, col.names=FALSE, row.names=FALSE)
tmp.f <- read.table(tmp.txt, col.names=c("STN", "WBAN", "YEARMODA", "TEMP", "TEMP.count", "DEWP", "DEWP.count", "SLP", "SLP.count", "STP", "STP.count", "VISIB", "VISIB.count", "WDSP", "WDSP.count", "MXSPD", "GUST", "MAX", "MIN", "PRCP", "SNDP", "FRSHTT"), stringsAsFactors=FALSE, as.is=TRUE)

## tidy up:
tmp.f$TEMPC <- round(ifelse(tmp.f$TEMP==9999.9, NA, (tmp.f$TEMP-32)*5/9), 1)  # convert to Celsius
tmp.f$PREC <- as.numeric(substr(tmp.f$PRCP, 1, nchar(tmp.f$PRCP)-1))
tmp.f$PREC <- round(ifelse(tmp.f$PREC==99.9, NA, tmp.f$PREC*25.4), 2)  # convert to mm
tmp.f$PREC.flag <- substr(tmp.f$PRCP, nchar(tmp.f$PRCP), nchar(tmp.f$PRCP))
tmp.f$PREC.flag <- ifelse(tmp.f$PREC.flag=="I"|tmp.f$PREC.flag=="A"|tmp.f$PREC.flag=="B"|tmp.f$PREC.flag=="C"|tmp.f$PREC.flag=="D"|tmp.f$PREC.flag=="E"|tmp.f$PREC.flag=="F"|tmp.f$PREC.flag=="G"|tmp.f$PREC.flag=="H", tmp.f$PREC.flag, "no flag") 
return(tmp.f)
unlink(tmp.txt)
}

str(GSOD.2008)  # 3,414,183 records!
GSOD.2008$DATE <- as.Date(paste(substr(GSOD.2008$YEARMODA, 1, 4), substr(GSOD.2008$YEARMODA, 5, 6), substr(GSOD.2008$YEARMODA, 7, 8), sep="-"))
GSOD.2008$STNID <- as.factor(paste(GSOD.2008$STN, GSOD.2008$WBAN, sep="-"))  # 9938 stations

## Attach coordinates and export to a table/shapefile
# read the location coordinates of stations:
download.file(paste(ftp.GSOD, "ish-history.csv", sep=""), destfile=paste(getwd(), "/ish-history.csv", sep=""), mode='wb', method='wget')
stations <- read.csv("ish-history.csv")
str(stations)
stations$STNID <- as.factor(paste(stations$USAF, stations$WBAN, sep="-"))  # WMO/DATSAV3 number
## format coordinates to decimal degrees:
stations$LAT <- ifelse(stations$LAT > 90.0*1000|stations$LAT < -90.0*1000, NA, stations$LAT/1000)
stations$LON <- ifelse(stations$LON > 180*1000|stations$LON < -180*1000, NA, stations$LON/1000)
stations$ELEV..1M. <- ifelse(stations$ELEV..1M.==-99999|stations$ELEV..1M.==-999.999, NA, stations$ELEV..1M./10)
plot(stations$LON, stations$LAT)

## merge two tables and write to a shape file:
GSOD.2008.XY <- merge(GSOD.2008, stations[,c("STNID", "STATION.NAME", "LAT", "LON", "ELEV..1M.")], by=c("STNID"))
# this takes cca 10 minutes!
str(GSOD.2008.XY)
GSOD.2008.XY <- subset(GSOD.2008.XY, !is.na(GSOD.2008.XY$LAT)&!is.na(GSOD.2008.XY$LON))
coordinates(GSOD.2008.XY) <- ~LON+LAT
proj4string(GSOD.2008.XY) <- CRS("+proj=longlat +datum=WGS84")

GSOD.20080501.XY <- subset(GSOD.2008.XY, GSOD.2008.XY$DATE==as.Date("2008-05-01"))
# bubble(GSOD.20080501.XY, "TEMPC")
writeOGR(GSOD.20080501.XY[,c("STNID", "TEMPC", "PREC")], "GSOD20080501_XY.kml", "GSOD20080501_XY", "KML")
writeOGR(GSOD.20080501.XY[,c("STNID", "YEARMODA", "TEMPC", "TEMP.count", "PREC", "PREC.flag")], "GSOD20080501_XY.shp", "GSOD20080501_XY", "ESRI Shapefile")
# unlink("GSOD_20080501_XY.shp")

# subset Bilogora case study:
bilogora.bbox <- data.frame(x=c(466500, 874500), y=c(4878500, 5286500))
coordinates(bilogora.bbox) <- ~x+y
proj4string(bilogora.bbox) <- CRS("+proj=utm +zone=33 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
bilogora.ll <- spTransform(bilogora.bbox, CRS("+proj=longlat +datum=WGS84"))
GSOD.2008.bilogora <- subset(GSOD.2008.XY, GSOD.2008.XY@coords[,1]>floor(bilogora.ll@coords[1,"x"])&GSOD.2008.XY@coords[,1]<ceiling(bilogora.ll@coords[2,"x"])&GSOD.2008.XY@coords[,2]>floor(bilogora.ll@coords[1,"y"])&GSOD.2008.XY@coords[,2]<ceiling(bilogora.ll@coords[2,"y"]))
# 38889 records
bubble(subset(GSOD.2008.bilogora, GSOD.2008.bilogora$DATE==as.Date("2008-05-01")&!is.na(GSOD.2008.bilogora$PREC)), "PREC")

# export as csv file:
write.table(GSOD.2008.bilogora@data, "GSOD_2008_bilogora.csv", row.names=FALSE, sep=";", quote=FALSE) 

# to read data to R:
# GSOD.2008.bilogora <- read.table("GSOD_2008_bilogora.csv", sep=";", header=TRUE)
