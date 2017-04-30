# IMAGE CRITICAL INFO TO FEED IN SEQUENCE OF R SCRIPTS

# INSTALL AND LOAD R PACKAGES
################################################################
ipak <- function(pkg){ # check to see if packages are installed. Install them if they are not, then load them into the R session.
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, repos="http://cran.rstudio.com/", dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
packages <- c("XML", "proj4", "rgdal", "raster")
ipak(packages)

# IMAGE INPUT
#################################################################
args <- commandArgs(trailingOnly = TRUE)
input <- Sys.getenv('input') # input="053734892180_01"
image <- paste(input, "_P001_MUL", ".tif", sep = "")
imagebase <- unlist(strsplit(image, "[.]"))[1]

# DEFINE WORKING DIRECTORIES
#################################################################
root <- Sys.getenv('baseDG')
#root <- "/home/stratouliasd/stars/derived/DG_v8" # FROM linux352.itc.utwente.nl
pathin <- paste(root, "0_categ", input, imagebase, sep = "/")
pathin_readme <- paste(root, "0_categ", input, sep = "/")
pathin_shp <- paste(root, "0_shapefiles", sep = "/")

# DEFINE PROJECTION SYSTEMS
#################################################################
lat_long <- "+proj=longlat +ellps=WGS84"
#UTM30N <- "+proj=utm +zone=30 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

# SHAPEFILES OF STUDY AREAS - CONVERT TO LAT/LONG
#################################################################

shp_latlong <- function(shp){ 
  setwd(pathin_shp)
  readshp <- readOGR(dsn = path.expand(pathin_shp), layer = shp)
  readshp@proj4string@projargs
  transf_latlong <- spTransform(readshp, CRS(lat_long))
}

# READ SHAPEFILES
MALI_latlong <- shp_latlong("MALI10_10_NEW_BOX")
NIGERIA_latlong <- shp_latlong("KKMa3_Kofa_box")
TANZANIA_NJOMBE_latlong <- shp_latlong("TZ_Njombe20x20km")
TANZANIA_KILOSA_latlong <- shp_latlong("TZ_Kilosa20x20km")
TANZANIA_SAME_latlong <- shp_latlong("TZ_Same20x20km")
UGANDA_latlong <- shp_latlong("UG_Moroto20x20km")
BANGLADESH_latlong <- shp_latlong("BD_6_IrMASaT_Sites_DG_20150128_10by10km")
MALI_BAMAKO_CAMPUS_latlong <- shp_latlong("ML_Bamako_campus")
TANZANIA_MOROGORO_CAMPUS_latlong <- shp_latlong("TZ_Morogoro_campus")

# GET EXTENT OF SHAPEFILES IN LAT LONG
NIGERIA_latlong_extent <- extent(NIGERIA_latlong@bbox[1,1], NIGERIA_latlong@bbox[1,2], NIGERIA_latlong@bbox[2,1], NIGERIA_latlong@bbox[2,2]) 
MALI_latlong_extent <- extent(MALI_latlong@bbox[1,1], MALI_latlong@bbox[1,2], MALI_latlong@bbox[2,1], MALI_latlong@bbox[2,2])
BANGLADESH_latlong_extent <- extent(BANGLADESH_latlong@bbox[1,1], BANGLADESH_latlong@bbox[1,2], BANGLADESH_latlong@bbox[2,1], BANGLADESH_latlong@bbox[2,2]) 
TANZANIA_KILOSA_latlong_extent <- extent(TANZANIA_KILOSA_latlong@bbox[1,1], TANZANIA_KILOSA_latlong@bbox[1,2], TANZANIA_KILOSA_latlong@bbox[2,1], TANZANIA_KILOSA_latlong@bbox[2,2]) 
TANZANIA_NJOMBE_latlong_extent <- extent(TANZANIA_NJOMBE_latlong@bbox[1,1], TANZANIA_NJOMBE_latlong@bbox[1,2], TANZANIA_NJOMBE_latlong@bbox[2,1], TANZANIA_NJOMBE_latlong@bbox[2,2]) 
TANZANIA_SAME_latlong_extent <- extent(TANZANIA_SAME_latlong@bbox[1,1], TANZANIA_SAME_latlong@bbox[1,2], TANZANIA_SAME_latlong@bbox[2,1], TANZANIA_SAME_latlong@bbox[2,2]) 
UGANDA_latlong_extent <- extent(UGANDA_latlong@bbox[1,1], UGANDA_latlong@bbox[1,2], UGANDA_latlong@bbox[2,1], UGANDA_latlong@bbox[2,2]) 
MALI_BAMAKO_CAMPUS_latlong_extent <- extent(MALI_BAMAKO_CAMPUS_latlong@bbox[1,1], MALI_BAMAKO_CAMPUS_latlong@bbox[1,2], MALI_BAMAKO_CAMPUS_latlong@bbox[2,1], MALI_BAMAKO_CAMPUS_latlong@bbox[2,2]) 
TANZANIA_MOROGORO_CAMPUS_latlong_extent <- extent(TANZANIA_MOROGORO_CAMPUS_latlong@bbox[1,1], TANZANIA_MOROGORO_CAMPUS_latlong@bbox[1,2], TANZANIA_MOROGORO_CAMPUS_latlong@bbox[2,1], TANZANIA_MOROGORO_CAMPUS_latlong@bbox[2,2]) 

# Plot maps
#windows()
#map('worldHires')
#plot(BANGLADESH_latlong, add=TRUE, col = "red", lwd=20)
#plot(NIGERIA_latlong, add=TRUE, col="red", lwd=20)
#plot(MALI_latlong, add=TRUE, col="red", lwd=20)

# DG IMAGES GET BAND NUMBER AND STUDY AREA VIA POINT IN POLYGON OPERATION
#################################################################

setwd(pathin)
file <- list.files(pattern = '*TIF$', , recursive = TRUE)[1]
image.info <- GDALinfo(file, returnStats=FALSE)
projection <- attributes(image.info)$projection
bandnu <- image.info[["bands"]]
  
 
setwd(pathin_readme)  
readmexml <- list.files(pattern = '*README.XML', recursive = F)
listoffiles <- lapply(readmexml, read.table, header=TRUE, colClasses="character", sep="\t")

xmlfiles <- xmlTreeParse(readmexml)
xmltop = xmlRoot(xmlfiles)
xmlvalues <- xmlSApply(xmltop, function(x) xmlSApply(x, xmlValue))
df <- data.frame(t(xmlvalues),row.names=NULL)
  
collection <- as.character(df$COLLECTIONSTART) 
date <- substring(collection, 1, 10)
substring(date, 5) <- "Y"
substring(date, 8) <- "M"
date <- paste(date, "D", sep = "")

timestart <- substring(collection, 12, 19)
timestart <- paste(substr(timestart, 1, 2), "hh", substr(timestart, 4, 5), "mm", substr(timestart, 7, 8), "ss", sep = "")
#timestart <- gsub(":", "_", timestart_raw)

# FIND EXTENT OF THE RASTER (DELIVERY) IN LAT LONG
NWLAT <- as.numeric(df$NWLAT)
NWLONG <- as.numeric(df$NWLONG)
SELAT <- as.numeric(df$SELAT)
SELONG <- as.numeric(df$SELONG)
raster_latlong_extent <- extent(NWLONG, SELONG, SELAT, NWLAT)


# FIND WHICH SHAPEFILE INTERSECTS WITH THE EXTENT OF THE RASTER LAYER
 if(!is.null(intersect(BANGLADESH_latlong_extent, raster_latlong_extent))){
  study_area <<- "BG"
  study_area_oid <<- 0000
  } else if(!is.null(intersect(NIGERIA_latlong_extent, raster_latlong_extent))){
  study_area <<- "NG_Kofa"
  study_area_oid <<- 1005
  } else if(!is.null(intersect(MALI_latlong_extent, raster_latlong_extent))){
  study_area <<- "ML_Sukumba"
  study_area_oid <<- 1000
  } else if(!is.null(intersect(TANZANIA_KILOSA_latlong_extent, raster_latlong_extent))){
  study_area <<- "TZ_Kilosa"
  study_area_oid <<- 1015 
  } else if(!is.null(intersect(TANZANIA_NJOMBE_latlong_extent, raster_latlong_extent))){
  study_area <<- "TZ_Njombe"
  study_area_oid <<- 1010 
  } else if(!is.null(intersect(TANZANIA_SAME_latlong_extent, raster_latlong_extent))){
  study_area <<- "TZ_Same"
  study_area_oid <<- 1020 
  } else if(!is.null(intersect(UGANDA_latlong_extent, raster_latlong_extent))){
  study_area <<- "UG_Moroto"
  study_area_oid <<- 1025 
  } else if(!is.null(intersect(MALI_BAMAKO_CAMPUS_latlong_extent, raster_latlong_extent))){
  study_area <<- "ML_Bamako_campus"
  study_area_oid <<- 1045 
  } else if(!is.null(intersect(TANZANIA_MOROGORO_CAMPUS_latlong_extent, raster_latlong_extent))){
  study_area <<- "TZ_Morogoro_campus"
  study_area_oid <<- 1050 
  } else {
  study_area <<- "Unknown"
  study_area_oid <<- 9999
  }
  
setwd(pathin) 
readimagexml <- list.files(pattern = '*P001.XML', recursive = T)
readimagexml <- Filter(function(x) grepl("-M", x), readimagexml)
listoffiles <- lapply(readimagexml, read.table, header=TRUE, colClasses="character", comment.char = "\t", sep="\t")

  xmlfiles <- xmlTreeParse(readimagexml)
  xmltop = xmlRoot(xmlfiles)
  xmlvalues <- xmlSApply(xmltop, function(x) xmlSApply(x, xmlValue))
  df <- data.frame(t(xmlvalues),row.names=NULL)
  
  sat <- as.character(xmlValue(xmlRoot(xmlfiles)[["IMD"]][["IMAGE"]][["SATID"]]))
    if (sat == "WV03"){
  sat <- "WorldView-3"
  } else{
	if (sat == "WV02"){
	sat <- "WorldView-2"
	} else{
	if (sat == "WV01"){
	sat <- "WorldView-1"
	} else{
	if (sat == "QB02") {
	sat <- "QuickBird"
	} else{
	if (sat == "GE01") {
	sat <- "GeoEye-1"
	} else{sat <- "Unknown"}
}}}}

  sunaz <- as.character(xmlValue(xmlRoot(xmlfiles)[["IMD"]][["IMAGE"]][["MEANSUNAZ"]]))
  sunel <- as.character(xmlValue(xmlRoot(xmlfiles)[["IMD"]][["IMAGE"]][["MEANSUNEL"]]))
  sataz <- as.character(xmlValue(xmlRoot(xmlfiles)[["IMD"]][["IMAGE"]][["MEANSATAZ"]]))
  satel <- as.character(xmlValue(xmlRoot(xmlfiles)[["IMD"]][["IMAGE"]][["MEANSATEL"]]))
  nadir <- as.character(xmlValue(xmlRoot(xmlfiles)[["IMD"]][["IMAGE"]][["MEANOFFNADIRVIEWANGLE"]]))
  cloudcover <- as.character(xmlValue(xmlRoot(xmlfiles)[["IMD"]][["IMAGE"]][["CLOUDCOVER"]]))
  
  #level <- as.character(xmlvalues[["IMD"]][["PRODUCTLEVEL"]])
  level <- as.character(xmlvalues[["IMD"]][["IMAGEDESCRIPTOR"]])
  
data <- rbind(input, sat, bandnu, sunaz, sunel, sataz, satel, nadir, study_area, level, date, timestart, cloudcover, study_area_oid)
rownames(data) <- c("Folder", "Satellite", "Number_of_bands", "Sun_azimuth", "Sun_elevation", "Satellite_azimuth", "Satellite_elevation", "off-nadir_angle", "Study_area", "Processing_level", "Acquisition_date", "Acquisition_time", "Cloud_cover", "Study_area_OID")
 
write.table(data, file = paste("metadata_", input, ".txt", sep = ""), sep="\t", row.names = TRUE, col.names = FALSE, quote = FALSE)
