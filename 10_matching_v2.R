#======================================================================================
# Variable definitions, data import, preparation
#======================================================================================
rm(list=ls(all=TRUE))

graphics.off()

require(rgdal)
require(Rcpp)

# IMAGE INPUT
#################################################################
args <- commandArgs(trailingOnly = TRUE)
input <- Sys.getenv('input') # input="053734892180_01"
image <- paste(input, "_P001_MUL", ".tif", sep = "")
imagebase <- unlist(strsplit(image, "[.]"))[1]

Server <- TRUE
WriteLog <- TRUE
Show_graph <- FALSE

if(Server){
	root <- Sys.getenv('baseDG')
	#root <- "/home/stratouliasd/stars/derived/DG_v8" # FROM linux352.itc.utwente.nl
	Path_lib <- Sys.getenv('basescript')
	pathmeta <- paste(root, "0_categ", input, imagebase, sep = "/")
	Path_in_s <- paste(root, "5_blob_detection", input,sep="/")
	Path_out <- paste(root, "6_image_matching",input,sep="/")
} else{
	root <- "S:/derived/DG_v9"
	Path_lib <- "S:/derived/scripts/v9"
	pathmeta <- paste("S:/derived/DG_v9/0_categ", input, imagebase, sep = "/")
}

setwd(Path_lib)
source("scale_space_lib_v7_1.R")
sourceCpp("matcher_v5.cpp")

# Define area
setwd(pathmeta)
metadata <- read.table(paste("metadata_", input, ".txt", sep = ""), stringsAsFactors=FALSE)
study_area <- as.character(unlist(metadata[9,2]))

# Define master image as reference
  if(study_area == "BG"){
  master_id <- "xxxxxxxxxxxxx"
  } else if(study_area == "NG_Kofa"){
  master_id <- "053734892380_01"
  } else if(study_area == "ML_Sukumba"){
  master_id <- "054112895030_01"
  } else if(study_area == "TZ_Kilosa"){
  master_id <- "xxxxxxxxxxxxx"
  } else if(study_area == "TZ_Njombe"){
  master_id <- "xxxxxxxxxxxxx"
  } else if(study_area == "TZ_Same"){
  master_id <- "xxxxxxxxxxxxx"
  } else if(study_area == "UG_Moroto"){
  master_id <- "xxxxxxxxxxxxx"
  } else {
  master_id <- "Unknown"
  }

imagebase <- unlist(strsplit(paste(master_id, "_P001_MUL", ".tif", sep = ""), "[.]"))[1]
pathmeta_master <- paste(root, "0_categ", master_id, imagebase, sep = "/") 
Path_images <- paste(root,"4_categ", study_area ,sep="/") # needed in subroutine
Path_in_m <- paste(root, "5_tree_measurement",master_id,sep="/")

setwd(Path_images)

# # exclude master image from processing and copy from 4_categ to 6_image_matching
# if(master_id == input){
	# if(WriteLog){
		# setwd(Path_out)
		# sink(file=paste("This is the master image, no reason of image registration",".txt",sep=""))
		# cat("This is the master image, no reason of image registration")
		# sink()}
		# setwd(Path_images)
		# master_to_copy <- list.files(".",pattern = paste(master_id, ".*TIF$", sep = ""), ignore.case = TRUE)
		# file.copy(from = master_to_copy, to = Path_out)
	# quit(save = "no")}
	
# # quit processing if there is no master image processed
# if(master_id == "xxxxxxxxxxxxx"){
	# if(WriteLog){
		# setwd(Path_out)
		# sink(file=paste("There is no master image processed for this area yet",".txt",sep=""))
		# cat("There is no master image processed for this area yet")
		# sink()}
	# quit(save = "no")}
	
setwd(Path_images)

# Read master and slave raster datasets
master_file <- list.files(".",pattern = paste(master_id, ".*TIF$", sep = ""), ignore.case = TRUE)
ms.imagefn_m <- master_file
ms.imageinfo_m <- GDALinfo(paste(Path_images,ms.imagefn_m,sep="/"))
ps.ms_m <- c(ms.imageinfo_m[["res.x"]],ms.imageinfo_m[["res.x"]])

slave_file <- list.files(".",pattern = paste(input, ".*TIF$", sep = ""), ignore.case = TRUE)
ms.imagefn_s <- slave_file
ms.imageinfo_s <- GDALinfo(paste(Path_images,"/",ms.imagefn_s,sep=""))
ps.ms_s <- c(ms.imageinfo_s[["res.x"]],ms.imageinfo_s[["res.x"]])

#===========================================================================
	# Read image geometry: sun and satellite position (determined from metadata)
	#===========================================================================

	# GET METADATA
	setwd(pathmeta_master)
	metadata <- read.table(paste("metadata_", master_id, ".txt", sep = ""), stringsAsFactors=FALSE)

	sun_sat <- as.numeric(unlist(metadata[4:8,2]))

	sun_az_m <- sun_sat[1] * pi/180
	sat_az_m <- sun_sat[3] * pi/180
	theta_sun_m <- (90-sun_sat[2]) * pi/180
	theta_sat_m <- (90-sun_sat[4]) * pi/180

	alpha_sun_m <- pi/2 - sun_az_m
	alpha_sat_m <- pi/2 - sat_az_m

	psi_m <- atan2(tan(theta_sat_m)*sin(alpha_sat_m)-tan(theta_sun_m)*sin(alpha_sun_m),tan(theta_sat_m)*cos(alpha_sat_m)-tan(theta_sun_m)*cos(alpha_sun_m))
	corr_factor_m <- cos(psi_m) / (tan(theta_sat_m)*cos(alpha_sat_m)-tan(theta_sun_m)*cos(alpha_sun_m))
	#===========================================================================

	# GET METADATA
	setwd(pathmeta)
	metadata <- read.table(paste("metadata_", input, ".txt", sep = ""), stringsAsFactors=FALSE)

	sun_sat <- as.numeric(unlist(metadata[4:8,2]))

	sun_az_s <- sun_sat[1] * pi/180
	sat_az_s <- sun_sat[3] * pi/180
	theta_sun_s <- (90-sun_sat[2]) * pi/180
	theta_sat_s <- (90-sun_sat[4]) * pi/180

	alpha_sun_s <- pi/2 - sun_az_s
	alpha_sat_s <- pi/2 - sat_az_s

	psi_s <- atan2(tan(theta_sat_s)*sin(alpha_sat_s)-tan(theta_sun_s)*sin(alpha_sun_s),tan(theta_sat_s)*cos(alpha_sat_s)-tan(theta_sun_s)*cos(alpha_sun_s))
	corr_factor_s <- cos(psi_s) / (tan(theta_sat_s)*cos(alpha_sat_s)-tan(theta_sun_s)*cos(alpha_sun_s))
	#===========================================================================
	# Define tiles
	Mtile <- 1000
	Ntile <- 1000
	# overlap of tiles; prevents loosing points at the margins
	tile_over <- 10

	# Open logfile
	if(WriteLog){
		setwd(Path_out)
		sink(file=paste("matching_log_master=",master_id,"_slave=",input,".txt",sep=""))
		cat(paste(Sys.time(),"starting","\n",sep=" "))
	}

	# Get image names from the input dir: filename contains tif but not aux
	setwd(Path_images)
	Files <- list.files(".",pattern=master_id,ignore.case = TRUE)

	tif_files <- grep(".tif",Files,ignore.case = TRUE)
	Files <- Files[tif_files]

	aux_files <- grep(".aux",Files,ignore.case = TRUE)
	if(length(aux_files)>0) Files <- Files[-aux_files]

	im_master <- 1

	ms.imagefn_m <- Files[im_master]

	ms.imageinfo_m <- GDALinfo(paste(Path_images,ms.imagefn_m,sep="/"))
	ps.ms_m <- c(ms.imageinfo_m[["res.x"]],ms.imageinfo_m[["res.x"]])

	Files <- list.files(".",pattern=input,ignore.case = TRUE)
	tif_files <- grep(".tif",Files,ignore.case = TRUE)
	Files <- Files[tif_files]

	aux_files <- grep(".aux",Files,ignore.case = TRUE)
	if(length(aux_files)>0) Files <- Files[-aux_files]

	im_slave <- 1

	ms.imagefn_s <- Files[im_slave]
	ms.imageinfo_s <- GDALinfo(paste(Path_images,"/",ms.imagefn_s,sep=""))
	ps.ms_s <- c(ms.imageinfo_s[["res.x"]],ms.imageinfo_s[["res.x"]])

	# Define tiles on the basis of master image
	N0.ms <- ms.imageinfo_m[["rows"]]
	M0.ms <- ms.imageinfo_m[["columns"]]

	# Read blobs from master image (already corrected for topographic shift)
	setwd(Path_in_m)
	datafile <- paste("corrected_trees_size_shape_top_100_perc_v6.RData",sep="")

	if(!file.exists(datafile)){
		if(WriteLog) cat(as.character(Sys.time()),"Error: blobs file for master image not found","\n")
		next
	}

	load(file=datafile)
	#Master_blobs <- Blobs_all
	Blobs_A <- Master_blobs

	# Read blobs from slave image
	setwd(Path_in_s)
	datafile <- paste("all_blobs_image_",input,".RData",sep="")

	if(!file.exists(datafile)){
		if(WriteLog) cat(as.character(Sys.time()),"Error: blobs file for slave image not found","\n")
		next
	}

	load(file=datafile)
	Blobs_B <- Blobs[,]

	Nb <- ms.imageinfo_m[["bands"]]

	if(Nb==8){
		# Set RGB composition
		nR <- 7
		nG <- 5
		nB <- 3
	}else{
		if(Nb==4){
			nR <- 4
			nG <- 3
			nB <- 2
		}else{
			nR <- 1
			nG <- 1
			nB <- 1
		}
	}

	N0.ms_m <- ms.imageinfo_m[["rows"]]
	M0.ms_m <- ms.imageinfo_m[["columns"]]

	N0.ms_s <- ms.imageinfo_s[["rows"]]
	M0.ms_s <- ms.imageinfo_s[["columns"]]

	Nb_s <- ms.imageinfo_s[["bands"]]
	if(Nb_s==8){
		# Set RGB composition
		nR_s <- 7
		nG_s <- 5
		nB_s <- 3
	}else{
		if(Nb_s==4){
			nR_s <- 4
			nG_s <- 3
			nB_s <- 2
		}else{
			nR_s <- 1
			nG_s <- 1
			nB_s <- 1
		}
	}


	ntx <- ceiling(M0.ms/Mtile)
	nty <- ceiling(N0.ms/Ntile)

	ix_arr <- 1:ntx
	iy_arr <- 1:nty

	if(WriteLog) cat(paste(Sys.time()," Tile size ",Mtile," by ",Ntile,"; tile overlap ",tile_over,"\n",sep=""))
	if(WriteLog) cat(paste(Sys.time()," Processing tiles ",min(ix_arr),":",max(ix_arr)," by ",min(iy_arr),":",max(iy_arr),"\n",sep=""))

	summary_matcher <- data.frame(array(0,c(0,6)))

	row_counter <- 0

	Matched_Blobs <- data.frame(array(0,c(0,13)))

	# process a single tile
	#ix <- 2
	#iy <- 4
	
	# Shift the master blobs according to slave image geometry
	h  <- Master_blobs$h
	hf <- Master_blobs$hf
	n  <- Master_blobs$n

	x <- Master_blobs$x
	y <- Master_blobs$y
	R <- 0.5 * Master_blobs$d

	h1 <- h*hf
	h2 <- h-h1

	start_time <- Sys.time()
	topo_shift_s <- project_pollock_quantitative_matrixC(h1, h2, R, theta_sat_s, n)
	Sys.time() - start_time

	x_s <- x - topo_shift_s * cos(alpha_sat_s)
	y_s <- y - topo_shift_s * sin(alpha_sat_s)

	Blobs_A$x <- x_s
	Blobs_A$y <- y_s

	for(ix in ix_arr)
	for(iy in iy_arr){
		i1 <- max((ix-1)*Mtile + 1 - tile_over,1)
		i2 <- min(ix*Mtile + tile_over,M0.ms_m) 
		j1 <- max((iy-1)*Ntile + 1 - tile_over,1)
		j2 <- min(iy*Ntile + tile_over,N0.ms_m)

		# read master image subset 
		ijr <- c(i1,i2,j1,j2)

		Path_in <- Path_images
		Msub <- read_subset(ms.imagefn_m,ijr[1],ijr[2],ijr[3],ijr[4])
		if(Show_graph){
			Mdisp <- Msub
			for(k in 1:Nb)Mdisp@data[,k] <- histstretch(Msub@data[,k])
		}

		if(FALSE){
			bb <- bbox(Msub)
			xrl <- bb[1,]
			yrl <- bb[2,]

			# read slave image subset 
			ijr <- xy_to_rowcol(cbind(xrl,yrl),ms.imageinfo_s)

			if(ijr[1]<0) ijr[1] <-1
			
			
			if(!((min(ijr)>0)&(ijr[2]<M0.ms_s)&(ijr[4]<N0.ms_s))){
				#cat(as.character(Sys.time()),"pol",i.pol,"Error: polygon is outside the slave image",im_slave,"\n")
				next
			}

			Ssub <- read_subset(ms.imagefn_s,ijr[1],ijr[2],ijr[3],ijr[4])

			if(Show_graph){
				Sdisp <- Ssub
				for(k in 1:Nb_s)Sdisp@data[,k] <- histstretch(Ssub@data[,k])
			}
		}


		r_search <- 40.0
		
		bb <- bbox(Msub)
		xr <- bb[1,]
		yr <- bb[2,]
		# Subset master blobs
		ind <- which((Blobs_A$x>xr[1])&(Blobs_A$x<xr[2])&(Blobs_A$y>yr[1])&(Blobs_A$y<yr[2]))
		Blobs_tile_A <- Blobs_A[ind,]
		# Subset slave blobs
		ind <- which((Blobs_B$x>xr[1]-r_search)&(Blobs_B$x<xr[2]+r_search)&(Blobs_B$y>yr[1]-r_search)&(Blobs_B$y<yr[2]+r_search))
		Blobs_tile_B <- Blobs_B[ind,]
		
		#  Clean the very long vector
		#ind <- which(Blobs_tile_B$ndvi>=0.6)
		#Blobs_tile_B <- Blobs_tile_B[ind,]
		
		Blobs_tile_A <- unique(Blobs_tile_A)
		Blobs_tile_B <- unique(Blobs_tile_B)
		
		consensus_set <- blob_match_v2(Blobs_tile_A,Blobs_tile_B,r_search,3.0,500)
		n_cons <- nrow(consensus_set)
		#n_cons

		if(length(consensus_set)==1) if(consensus_set==-1) next

		Blobs_A1 <- Blobs_tile_A[consensus_set[,1],]
		Blobs_B1 <- Blobs_tile_B[consensus_set[,2],]

		xydiff <- Blobs_A1[,1:2]-Blobs_B1[,1:2]
		dxy <- colMeans(xydiff)

		dx <- dxy[1]
		dy <- dxy[2]

		rmse <- sqrt(mean((xydiff[,1]-dx)^2 + (xydiff[,2]-dy)^2))

		# Apply transform and show results
		Btransf <- geo_transform(Blobs_B1,dx,dy)
		
		if(Show_graph){
			#windows()
			par(mfrow=c(1,2))
			image(Mdisp,red=nR,green=nG,blue=nB,axes=TRUE)
			title(main=paste("master image, tile x=",ix,", y=",iy,sep=""))
			#display_all_blobs(Blobs_m,"green")
			display_all_blobs(Blobs_A1,"green")
			#if(nrow(Blobs_A1)>0)text(Blobs_A1$x,Blobs_A1$y,labels=1:nrow(Blobs_A1),pos=4)

			#image(Sdisp,red=nR_s,green=nG_s,blue=nB_s,axes=TRUE)
			image(Mdisp,red=nR,green=nG,blue=nB,axes=TRUE)
			title(main=paste("slave image ",sep=""))
			display_all_blobs(Blobs_B1,"blue")
			#if(nrow(Blobs_B1)>0)text(Blobs_B1$x,Blobs_B1$y,labels=1:nrow(Blobs_B1),pos=4)

			display_all_blobs(Btransf,"white")
		}

		newdf <- Blobs_A1
		for(i in 1:5)names(newdf)[i] <- paste(names(newdf)[i],"_m",sep="") 

		newdf <- cbind(newdf,Blobs_B1)
		for(i in (ncol(Blobs_tile_A)+1):(ncol(Blobs_tile_A)+ncol(Blobs_B1)))names(newdf)[i] <- paste(names(newdf)[i],"_s",sep="") 

		row_first <- row_counter+1
		row_last  <- row_counter+nrow(newdf)
		rownames(newdf) <- row_first:row_last
		row_counter <- row_last

		if((ix==1)&(iy==1)){
			Matched_Blobs <- newdf
		}else{
			Matched_Blobs <- rbind(Matched_Blobs,newdf)
		}
		
		if(WriteLog) cat(as.character(Sys.time())," tile ix=",ix," iy=",iy,"success", dx,dy,n_cons,rmse,row_first,row_last,"\n")
		summary_matcher <- rbind(summary_matcher,c(dx,dy,n_cons,rmse,row_first,row_last))
	}	

	names(summary_matcher) <- c("dx","dy","n_cons","rmse","row_start","row_end")

	setwd(Path_out)
	save(Matched_Blobs,file=paste("matched_blobs_master_",im_master,"_slave_",im_slave,".RData",sep=""))
	save(summary_matcher,file=paste("summary_of_matching",".RData",sep=""))
	

	setwd(Path_out)
	load(file=paste("matched_blobs_master_",im_master,"_slave_",im_slave,".RData",sep=""))
	load(file=paste("summary_of_matching",".RData",sep=""))

	#str(Matched_Blobs)
	#summary_matcher

	# Apply topographic shift of master GCPs
	dx_m <- 0
	dy_m <- 0
	# Shift master to fit the FMU polygons
	#dx_m <- -3.5*ps.ms_m[1]
	#dy_m <- +2.0*ps.ms_m[2]

	Blobs_A <- data.frame(cbind(Matched_Blobs$x_m,Matched_Blobs$y_m,Matched_Blobs$d_m))
	names(Blobs_A) <- c("x","y","d")

	Blobs_B <- data.frame(cbind(Matched_Blobs$x_s,Matched_Blobs$y_s,Matched_Blobs$d_s))
	names(Blobs_A) <- c("x","y","d")

	dx_all <- Matched_Blobs$x_m - Matched_Blobs$x_s
	dy_all <- Matched_Blobs$y_m - Matched_Blobs$y_s

	dx <- median(dx_all)
	dy <- median(dy_all)

	if(WriteLog) cat(paste(Sys.time()," Determined overall transformation parameters: dx=",dx," dy=",dy,"\n",sep=""))


	#		Blobs_A$xtrue <- Blobs_A$x + dx_m
	#		Blobs_A$ytrue <- Blobs_A$y + dy_m

	#r_search <- 40
	#consensus_set <- blob_match_v2(Blobs_A[1:15000,],Blobs_B[1:15000,],r_search,3.0,500)
	#n_cons <- nrow(consensus_set)
	#====================================================================
	# Transform entire scenes
	#====================================================================
	if(FALSE){
		Msub <- readGDAL(paste(Path_images,"/",ms.imagefn_m,sep=""))

		# transform the master image
		cellcentre.offset <- Msub@grid@cellcentre.offset + c(dx_m,dy_m)
		m_grid <- GridTopology(cellcentre.offset, Msub@grid@cellsize, Msub@grid@cells.dim)
		Mtransf <- SpatialGridDataFrame(m_grid, Msub@data, proj4string = Msub@proj4string)		

		# Save transformed image
		
		text_str <- unlist(strsplit(ms.imagefn_m, split = "[.]"))[[1]]
		text_str <- paste(text_str,"_transformed_topocorr",".tif",sep="")

		setwd(Path_out)
		OUT <- Mtransf
		OUT.tif<-create2GDAL(OUT,drivername="GTiff",type="UInt16")
		saveDataset(OUT.tif,text_str)
		GDAL.close(OUT.tif)
	}

	if(WriteLog) cat(paste(Sys.time()," Opening slave image","\n",sep=""))

	# transform the slave image
	Ssub <- readGDAL(paste(Path_images,"/",ms.imagefn_s,sep=""))

	if(WriteLog) cat(paste(Sys.time()," Applying transformation to slave image","\n",sep=""))

	dx_s <- dx + dx_m
	dy_s <- dy + dy_m

	# transform the slave image
	cellcentre.offset <- Ssub@grid@cellcentre.offset + c(dx_s,dy_s)
	s_grid <- GridTopology(cellcentre.offset, Ssub@grid@cellsize, Ssub@grid@cells.dim)
	Stransf <- SpatialGridDataFrame(s_grid, Ssub@data, proj4string = Ssub@proj4string)		

	# Save transformed image

	text_str <- unlist(strsplit(ms.imagefn_s, split = "[.]"))[[1]]
	text_str <- paste(text_str,"_transformed_topocorr",".tif",sep="")

	if(WriteLog) cat(paste(Sys.time()," Saving corrected slave image","\n",sep=""))

	setwd(Path_out)
	OUT <- Stransf
	OUT.tif<-create2GDAL(OUT,drivername="GTiff",type="UInt16")
	saveDataset(OUT.tif,text_str)
	GDAL.close(OUT.tif)

	# Close logfile
	if(WriteLog){
		cat(paste(Sys.time(),"Process ended",sep=" "))
		sink()
	}

# }
#==================================================================================
# The End
#==================================================================================
