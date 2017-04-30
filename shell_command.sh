#!/bin/bash
#bash
#nohup ./shell_command.sh &>/dev/null &
#cd /home/stratouliasd/stars/derived/scripts/v9

# SET PATHNAME 
#This is non-permanent change, LD_LIBRARY_PATH reset to default value at next session)
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib # extends the session's library-path
# export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu:/usr/local/lib:/usr/lib
export base=~/stars/derived
export baseDG=~/stars/derived/DG_v9
export basescript=~/stars/derived/scripts/v9
export R_LIBS=~/stars/rlibs

# SET IMAGE FILENAME
export input=$1
# export input="054330851010_01" # BG ortho shift 
# export input="054393058010_01" # BG ortho shift
export input="054294362010_01" # Nigeria cssl issue
# export input="053613698020_01" # Mali
echo $input

# COPY DATA FROM ARCHIVE TO PROCESSING NODE
# Rscript $basescript/0_copy_data.R $input

# DERIVE METADATA
# Rscript $basescript/1_metadata.R $input

# ASSIGN VARIABLE TO METADATA ELEMENTS
META=($input)
i=0
while IFS=$'\t' read -r -a myArray
do
 META[++i]=${myArray[1]}
 echo "${myArray[1]}"
done < $baseDG/0_categ/$input/''$input''_P001_MUL/metadata_$input.txt

# PROCESS BY AREA
AREA="NG_Kofa" # ML_Sukumba   |   NG_Kofa
if [ "${META[9]}" = "${AREA}" ]; then
  echo "Study area is ${AREA}"
else
  echo "Study area is not ${AREA}"
  exit
fi 

# CREATE IMAGE-SPECIFIC DIRECTORIES
# mkdir -p $baseDG/1_atcor_6s/"$input"
# mkdir -p $baseDG/2_mosaic_r/"$input"
# mkdir -p $baseDG/3_orthorectification_otb/"$input"
# mkdir -p $baseDG/5_blob_detection/"$input"
# mkdir -p $baseDG/6_image_matching/"$input"
# mkdir -p $baseDG/7_cloud_mask/"$input"
# mkdir -p $baseDG/8_tree_mask/"$input"

# # ATMOSPHERIC CORRECTION
# cd $baseDG/1_atcor_6s_source/PythonWrapper
# python $baseDG/1_atcor_6s_source/PythonWrapper/pythonWrV8.py --band4n8alldir=''$baseDG'/0_categ/'$input'/'$input'_P001_MUL' --outputdir=''$baseDG'/1_atcor_6s/'$input'/'

# CONVERT TO .TIF
# Rscript $basescript/3_export_tif.R $input

# STITCH TILES
# Rscript $basescript/4_mosaic.R $input

# PREPARE AND RENAME FILES
# Rscript $basescript/5_rename.R $input

# ORTHO-RECTIFICATION 
if [ ${META[10]} = "OrthoRectified3" ]
then
  echo "OrthoRectified3"
  # /usr/local/bin/otbcli_OrthoRectification -io.in ""$baseDG"/2_mosaic_r/"$input"/"$input"_P001_M2AS.tif?&skipcarto=true" -io.out ""$baseDG"/3_orthorectification_otb/"$input"/"$input"_otb_skipcarto_DEM_nn.tif" uint16 -elev.dem ""$base"/DEM/DEM_SRTM" -interpolator nn
elif [ ${META[10]} = "Standard2A" ]
then
  echo "Standard2A"
  # /usr/local/bin/otbcli_OrthoRectification -io.in ""$baseDG"/2_mosaic_r/"$input"/"$input"_P001_M2AS.tif?&skipcarto=true" -io.out ""$baseDG"/3_orthorectification_otb/"$input"/"$input"_otb_skipcarto_DEM_nn.tif" uint16 -elev.dem ""$base"/DEM/DEM_SRTM" -interpolator nn
else [ ${META[10]} = "ORStandard2A" ]
  echo "ORStandard2A"
  # /usr/local/bin/otbcli_OrthoRectification -io.in ""$baseDG"/2_mosaic_r/"$input"/"$input"_P001_M2AS.tif?&skipcarto=true" -io.out ""$baseDG"/3_orthorectification_otb/"$input"/"$input"_otb_skipcarto_DEM_nn.tif" uint16 -elev.dem ""$base"/DEM/DEM_SRTM" -interpolator nn
fi

# RENAME AND CATEGORIZE ACORDING TO AREA OF STUDY
# Rscript $basescript/7_categorize_area.R $input

# BLOB DETECTION
# Rscript $basescript/8_blob_detector_server_clean_v8.R $input 
 
##### THE MASTER IMAGE FOR EACH SITE IS PROCESSSED AT 5_TREE_MEASUREMENT

# IMAGE MATCHING
# Rscript $basescript/10_matching_v2.R $input

# CLOUD MASKING
# Rscript $basescript/11_cloud_masking_v3.r $input

# TREE MASKING
# Rscript $basescript/12_tree_masking.R $input

# EXTRACT POLYGON-BASED SPECTRAL STATISTICS FOR THE CSSL
# Rscript $basescript/13_cssl_statistics_spectral_v7.R $input

# EXTRACT POLYGON-BASED TEXTURE STATISTICS FOR THE CSSL
# Rscript $basescript/13_cssl_statistics_textural_v4.R $input

# # # #################################
# # # #report  processing time
# # # # proc.time()

# # # #ls -LR | grep .tif$ >> list_tif.txt # list all tif files
# # # #find . -maxdepth 1 -type d >> subdirectories.txt # list subdirectories
# # # # find . -maxdepth 3 -iname "metadata*" > 0_metadata.txt
# # # # find . -iname "*.tif" | cut -d'/' -f 3-4 > list_tif.txt
# # # #ls -l --block-size=MB -LR | grep .tif$ >> tif_size.log # size of .tif files
# # # #du -h --max-depth=1 | sort -hr >> dir_size.txt # size of directories
