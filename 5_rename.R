# COPY AUX FILES AND FILE RENAME 

# REMOVE ALL THE VARIABLES FROM THE WORKSPACE
#################################################################
rm(list=ls(all=TRUE))

# IMAGE INPUT
#################################################################
args <- commandArgs(trailingOnly = TRUE)
input <- Sys.getenv('input') # input="053734892180_01"
image <- paste(input, "_P001_MUL", ".tif", sep = "")
imagebase <- unlist(strsplit(image, "[.]"))[1]

# INSTALL AND LOAD R PACKAGES
#################################################################
require("rgdal")
require("raster")

# DEFINE WORKING DIRECTORIES
#################################################################
root <- Sys.getenv('baseDG')
#root <- "/home/stratouliasd/stars/derived/DG_v8" # FROM linux352.itc.utwente.nl
pathin <- paste(root, "0_categ", input, imagebase, sep = "/")
pathout <- paste(root, "2_mosaic_r", input, sep = "/")
pathout_otb <- paste(root, "3_orthorectification_otb", input, sep = "/")
dempath <- paste(dirname(root), "DEM", "DEM_SRTM", sep = "/")

# COPY .RPB FILE
setwd(pathin)
listRPBraw <- list.files(pattern = '*RPB$')
invisible(sapply(listRPBraw, file.copy, to = pathout, recursive = TRUE, overwrite = FALSE))
setwd(pathout)
listRPB <- list.files(pattern = '*RPB$')
file.rename(listRPB, gsub("[-]", "_", listRPB))
listRPB <- list.files(pattern = '*RPB$')
listTIF <- list.files(pattern = '*tif$', ignore.case = TRUE)

# SET FILENAME
filename_split <- unlist(strsplit(listRPB[[1]], "[.]"))
name_split <- unlist(strsplit(filename_split[[1]], "[_]"))
length(name_split)

if (length(name_split) == 5) {
	#pat <- paste(name_split[[3]], name_split[[4]], sep = "_")
	filename_out <- paste(name_split[[3]], name_split[[4]], name_split[[5]], name_split[[2]], sep = "_")
}	else {
		pat <- "Unknown naming convention"
		filename_out <- "Unknown naming convention"
		}	
	
filename_rpb <- paste(filename_out, ".RPB", sep = "")
filename_tif <- paste(filename_out, ".tif", sep = "")

# # Check if there is only one tile or several in the delivery - i.e. if there is R1C1 etc......
# if (length(name_split) == 7) {
	# pat <- paste(name_split[[3]], name_split[[4]], sep = "_")
	# filename_out <- paste(name_split[[5]], name_split[[6]], name_split[[7]], sep = "_")
# }	else if (length(name_split) == 6) {
	# pat <- name_split[[3]]
	# filename_out <- paste(name_split[[4]], name_split[[5]], name_split[[6]], sep = "_")
	# }	else {
		# pat <- "Unknown naming convention"
		# filename_out <- "Unknown naming convention"
		# }	
		
# RENAME .RPB and .tif FILES
file.rename(from = listRPB, to = filename_rpb)
file.rename(from = listTIF, to = filename_tif)

# BUILD OTB BATCH FILE FOR THE ORTHORECTIFICATION
#################################################################

# Header
header <- capture.output(cat(
  "# ORTHORECTIFICATION IN OTB 
\n# !! skip cartographic information for ortho-ready products (?&skipcarto=true)
# !! use SRTM DEM downloaded from OTB
# !! change working directory to:
#D:
#cd S:\\derived
#for Windows: type the name of the file followed by .bat in the OSGeo4W shell
\n\n"))

# Individual images
command_image <- function(){ # c
  paste("call otbcli_OrthoRectification", " ",
        "-io.in", " ", "\"", pathout, "/", imagebase, ".tif", "?&skipcarto=true", "\"", " ",
        "-io.out"," ", pathout_otb, "/", input, "_", "otb_cli_skipcarto_DEM", "_", ".tif",  " ",
        "-elev.dem", " ", dempath, " ",
        "-interpolator nn",
        sep = "")}
commands <- command_image()

# Combine header and commands
batch_file <- c(header, commands)

setwd(pathout_otb)
write.table(batch_file, paste(imagebase, ".bat", sep = ""), row.names = FALSE, col.names = FALSE, eol = "\n", quote = FALSE)