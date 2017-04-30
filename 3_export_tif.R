# EXPORT FROM HDF TO TIF AND APPEND METADATA
rm(list=ls(all=TRUE))

# INSTALL AND LOAD R PACKAGES
#################################################################
require("rgdal")
require("gdalUtils")

#drivers <- gdalDrivers() 
#any(drivers$name %in% c("HDF4", "HDF"))

args <- commandArgs(trailingOnly = TRUE)
input <- Sys.getenv('input') # input="053734892180_01"
image <- paste(input, "_P001_MUL", ".tif", sep = "")
imagebase <- unlist(strsplit(image, "[.]"))[1]

# DEFINE WORKING DIRECTORIES
#################################################################
root <- Sys.getenv('baseDG')
#root <- "/home/stratouliasd/stars/derived/DG_v8" # FROM linux352.itc.utwente.nl
pathin <- paste(root, "1_atcor_6s", input, sep = "/")
pathin_ref <- paste(root, "0_categ", input, imagebase, sep = "/")

setwd(pathin)
list_hdf <- list.files(pattern = "\\.hdf$")

if(length(list_hdf) == 0) { stop("No .HDF files found")}

# REPLACE "-"WITH "_" IN ALL .HDF4 FILENAMES
for (i in 1:length(list_hdf)) {
hdf <- list_hdf[i]
hdfinfo <- gdalinfo(list_hdf[i])
file.rename(list_hdf[i], gsub("[-]", "_", list_hdf[i]))
}

#RE-READ RENAMED LIST OF .HDF4 FILES
list_hdf <- list.files(pattern = "\\.hdf$")

for (i in 1:length(list_hdf)) {
	x <- list_hdf[[i]]
	info <- gdalinfo(x)
	data <- grep("HDF4_SDS", info, value = TRUE)
	datasub <- gsub("SUBDATASET_\\d+_NAME=|\\s+", "", data)
	
	filename_split <- unlist(strsplit(x[[1]], "[.]"))
	name_split <- unlist(strsplit(filename_split[[2]], "[_]"))
	length(name_split)

	# Check if there is only one tile or several in the delivery - i.e. if there is R1C1 etc......
	if (length(name_split) == 11) {
		pat <- paste(name_split[[7]], name_split[[8]], sep = "_")
		filename_out <- paste(name_split[[9]], name_split[[10]], name_split[[11]], sep = "_")
	}	else if (length(name_split) == 10) {
		pat <- name_split[[7]]
		filename_out <- paste(name_split[[8]], name_split[[9]], name_split[[10]], sep = "_")
		}	else {
			pat <- "Unknown naming convention"
			filename_out <- "Unknown naming convention"
			}

	# GET METADATA FROM RAW IMAGE
	setwd(pathin_ref)

	ref_tile_name <- list.files(pattern = paste(pat, ".*TIF$", sep = ""), recursive = FALSE)
	ref_tile <- readGDAL(ref_tile_name)
	nb <- length(ref_tile@data)
	bounds_ullr <- c(ref_tile@bbox[1,1], ref_tile@bbox[2,2], ref_tile@bbox[1,2], ref_tile@bbox[2,1])
	
	my_data <- data.frame(array(0,c(nrow(ref_tile),nb)))
	
	for(j in 1:nb)names(my_data)[j] <- paste("band",j,sep="")
	A <- SpatialGridDataFrame(ref_tile@grid, my_data, ref_tile@proj4string)
	# CONTINUE WITH CONVERSION
	setwd(pathin)

	# EXPORT INDIVIDUAL HDF4 FILES TO TIF
	temp <- paste(pathin, "/temp_", i, "_", pat, sep = "")
	dir.create(temp, showWarnings = TRUE, recursive = FALSE, mode = "0777")
	for (j in 1:length(datasub)) {
		hdftotif <- gdal_translate(src_dataset = x, dst_dataset = paste(temp, "/band0", j, "_surface_reflectance", ".tif", sep = ""), of = "GTiff", ot = "UInt16", output_Raster = TRUE, sd_index=j, a_srs = ref_tile@proj4string, a_ullr = bounds_ullr)
	}
	
	# STACK INDIVIDUAL BANDS FROM HDF
	for (j in 1:nb) {
		y <- readGDAL(datasub[j])
		tmp_band <- y@data$band1	
		tmp_band[tmp_band<0] <- 2^16-1	
		A@data[,j] <- tmp_band
	}

	writeGDAL(dataset = A, fname = paste(filename_out, "_", pat, ".tif", sep = ""), drivername = "GTiff", type = "UInt16", mvFlag = 2^16-1)
	#x <- readGDAL(fname = paste(filename_out, "_", pat, ".tif", sep = ""))
	
	#unlink(temp, recursive = TRUE, force = FALSE) # delete content of the directory!!
}
