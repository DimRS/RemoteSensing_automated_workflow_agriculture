# CLASSIFY THE FOLDERS ACCORDING TO THE AREA OF STUDY

# REMOVE ALL THE VARIABLES FROM THE WORKSPACE
#################################################################
rm(list=ls(all=TRUE))

Server <- TRUE

# INSTALL AND LOAD R PACKAGES
#################################################################
require("rgdal")
require("proj4")
require("gdalUtils")

# IMAGE INPUT
#################################################################
args <- commandArgs(trailingOnly = TRUE)
# input <- args[1]
input <- Sys.getenv('input') # input="053734892180_01"
image <- paste(input, "_P001_MUL", ".tif", sep = "")
imagebase <- unlist(strsplit(image, "[.]"))[1]

# DEFINE WORKING DIRECTORIES
#################################################################
root <- Sys.getenv('baseDG')
#root <- "/home/stratouliasd/stars/derived/DG_v9" # FROM linux352.itc.utwente.nl
pathmeta <- paste(root, "0_categ", input, imagebase, sep = "/")
pathin <- paste(root, "3_orthorectification_otb", input, sep="/")
pathin_alt <- paste(root, "2_mosaic_r", input, sep="/")
pathout <- paste(root, "4_categ", sep = "/")

# GET METADATA
#################################################################
setwd(pathmeta)
metadata <- read.table(paste("metadata_", input, ".txt", sep = ""))
satellite <- as.character(metadata[2,2])
study_area <- as.character(metadata[9,2])
proc_level <- as.character(metadata[10,2])
acq_date <- as.character(metadata[11,2])
Band_nr <- as.character(metadata[3,2])
pathout_area <- paste(pathout, "/", study_area, "/", sep = "")

if(Band_nr==4)
  {nR_false <- 4
   nG_false <- 3
   nB_false <- 2
   nR_true <- 3
   nG_true <- 2
   nB_true <- 1
  } else if(Band_nr==8) 
    {nR_false <- 7
     nG_false <- 5
     nB_false <- 3
	 nR_true <- 5
	 nG_true <- 3
	 nB_true <- 2	
    } else if(Band_nr==5) 
      {nR_false <- 5
       nG_false <- 3
       nB_false <- 2
	   nR_true <- 3
	   nG_true <- 2
	   nB_true <- 1
      } else 
        {nR_false <- "Unknown"
         nG_false <- "Unknown"
         nB_false <- "Unknown"
		 nR_true <-  "Unknown"
	     nG_true <-  "Unknown"
	     nB_true <-  "Unknown"
         }

filename_out <- paste(input, study_area, acq_date, satellite, proc_level, sep = "_")

setwd(pathin)
listtif <- list.files(pattern = '*.tif$', recursive = T, ignore.case = TRUE)

# CATEGORIZE ACCORDING TO STUDY AREA
#################################################################
if (length(listtif) == 0) {
	setwd(pathin_alt)
	listtif_alt <- list.files(pattern = '*.tif', recursive = T, ignore.case = TRUE)
	sapply(listtif_alt, file.copy, to = pathout_area, recursive = TRUE, overwrite = FALSE) 
	} else {
		setwd(pathin)
		file.copy(from = listtif, to = pathout_area)
		setwd(pathout_area)
		file.rename(listtif, paste(filename_out, ".tif", sep = ""))
		#setwd(pathout_area)
		#gdal_translate(paste(pathin, listtif, sep="/"), paste(filename_out, ".tif", sep = ""), of = "GTiff", ot = "UInt16")
	  
# # GENERATE THUMBNAILS
# #################################################################

# # ATTEMPT 1: KLM export
# # Linear stretch of a multi-band raster	 
# linstretch<-function(img,minmax=NA){
  # if(is.na(minmax)) minmax<-c(min(getValues(img),na.rm=T),max(getValues(img),na.rm=T))
  # temp<-calc(img,fun=function(x) (255*(x-minmax[1]))/(minmax[2]-minmax[1]))
  # #set all values above or below minmax to 0 or 255
  # temp[temp<0]<-0;temp[temp>255]<-255;
  # return(temp)
# }
# # Histogram equalization stretch of a multi-band raster
# eqstretch<-function(img){
  # unique <- na.omit(getValues(img))
  # if (length(unique>0)){
    # ecdf<-ecdf(unique)
    # out <- calc(img,fun=function(x) ecdf(x)*255)
  # }
  # return(out)
# }
	  
	  # projection(stacked)
	  # p <- projectRaster(rasterized, crs="+proj=longlat +datum=WGS84", method='bilinear') # transform to longitude/latitude
	  # p <- linstretch(p)
	  # KML(p, file='test.kml')

# ATTEMPT 2: readGDAL and image 
	  histstretch<-function(data) {
	    cur.lim<-quantile(data,c(0.025,0.975),na.rm=TRUE)
	    data<-pmax(cur.lim[1],pmin(cur.lim[2],data))
	    data<-floor(255*(data-cur.lim[1])/(cur.lim[2]-cur.lim[1]))
	    data[is.na(data)]<-0
	    data
	  }
	 
	  #READ FALSE-COLOUR COMPOSITE 	 
	  false_colour <- readGDAL(paste(pathin, listtif, sep="/"))
	  false_colour@data <- data.frame(false_colour@data[,c(nR_false,nG_false,nB_false)])
	  for(nb in 1:3){
	    false_colour@data[,nb] <- histstretch(false_colour@data[,nb])
	  }
	  
	  #READ TRUE-COLOUR COMPOSITE 	 
	  true_colour <- readGDAL(paste(pathin, listtif, sep="/"))
	  true_colour@data <- data.frame(true_colour@data[,c(nR_true,nG_true,nB_true)])
	  for(nb in 1:3){
	    true_colour@data[,nb] <- histstretch(true_colour@data[,nb])
	  }
	  
	  # GENERATE FALSE-COLOUR COMPOSITE PREVIEW IMAGE
	  #jpeg(filename = paste("preview_", filename_out, ".jpeg", sep = ""), width = 1000, height = 1000, units = "px", quality = 65, type = "cairo")
	  png(filename = paste("preview_", "falsecolour_", filename_out, ".png", sep = ""), width = 500, height = 500, units = "px", bg = "white", type = "cairo")
	  par(mar=c(0,0,0,0), mai=c(0,0,0,0))
	  image(false_colour, red=1, green=2, blue=3, axes=FALSE, box=FALSE)
	  #plotRGB(false_colour, 1, 2, 3, stretch='lin')
	  #title(main=filename_out)
	  #dev.copy(jpeg,"thumbnail.jpg")
	  dev.off()
	  # GENERATE TRUE-COLOUR COMPOSITE PREVIEW IMAGE
	  png(filename = paste("preview_", "truecolour_", filename_out, ".png", sep = ""), width = 500, height = 500, units = "px", bg = "white", type = "cairo")
	  par(mar=c(0,0,0,0), mai=c(0,0,0,0))
	  image(true_colour, red=1, green=2, blue=3, axes=FALSE, box=FALSE)
	  dev.off()
}
