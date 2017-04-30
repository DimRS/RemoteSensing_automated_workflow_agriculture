# TILES STITCHING IN R 

# REMOVE ALL THE VARIABLES FROM THE WORKSPACE
#################################################################
rm(list=ls(all=TRUE))
WriteLog <- TRUE

# IMAGE INPUT
#################################################################
args <- commandArgs(trailingOnly = TRUE)
input <- Sys.getenv('input') # input="053734892180_01"
image <- paste(input, "_P001_MUL", ".tif", sep = "")
imagebase <- unlist(strsplit(image, "[.]"))[1]

# INSTALL AND LOAD R PACKAGES
#################################################################
require("rgdal")
require("sp")
require("spatial.tools")
require("proj4")
require("gdalUtils")

# DEFINE WORKING DIRECTORIES
#################################################################
root <- Sys.getenv('baseDG')
#root <- "/home/stratouliasd/stars/derived/DG_v8" # FROM linux352.itc.utwente.nl
pathin <- paste(root, "1_atcor_6s", input, sep = "/")
pathout <- paste(root, "2_mosaic_r", input, sep = "/")

# Open logfile
if(WriteLog){
	setwd(pathout)
	sink(file = paste(input, "_mosaic",".log", sep=""), split = TRUE)
	cat(paste("started at",Sys.time(),"\n",sep=" "))
}

setwd(pathin)

list_files <- list.files(pattern = "\\.tif$", recursive = FALSE, ignore.case = TRUE)

# set filename
filename <- paste(strsplit(list_files[[1]], "[_]")[[1]][1], strsplit(list_files[[1]], "[_]")[[1]][2], strsplit(list_files[[1]], "[_]")[[1]][3], strsplit(list_files[[1]], "[_]")[[1]][4], sep = "_")

if(length(list_files) == 0) { stop("No .tif files found")}
if(WriteLog) {cat(paste("The following raster tiles are mosaiced: ", "\n",sep = ""))}
if(WriteLog) {cat(paste(writeLines(list_files), "\n", sep = ""))}

if (length(list_files)==1){
	file.copy(list_files[1], to = pathout, recursive = TRUE)
} else{
	x <- gdalwarp(srcfile = list_files, dstfile = paste(pathout, "/" , filename, ".tif", sep = ""), of = "GTiff", type = "UInt16")
} 

# Close logfile
if(WriteLog){
	cat(paste("ended at", Sys.time(),sep=" "))
	sink()
} 
# 