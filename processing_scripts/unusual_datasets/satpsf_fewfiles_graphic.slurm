#!/bin/bash
#SBATCH -p daceg -o my.stdout --mail-user=anthony.cheetham@unige.ch --mail-type=ALL --error=error_output_%j --ntasks=10 --mincpus=10

# module add openmpi/gcc/1.8.4
# module add python/2.7.8

#Steps for analysis:
# 1. rm_overscan to make images square
# 2. register to measure the psf centre etc.
# 3. NaNmedian on sky frames
# 4. sky_sub to subtract median sky from data
# 5. bad_pix 
# 6. medianfilter to smooth the images
# 6b. (Optional) Bin the data
# 7. psf_sub to subtract the psf (LOCI)
# 8. deromed to derotate and combine the images.

# Steps 1-5 need to be repeated for the flux frames.
# If the SNR is low, steps 2-4 may need to be repeated. (i.e. 1234 234 5678)

# This is for targets that don't have many frames. It uses simple sky subtraction (closest in time) and skips frame selection

# GRAPHIC_VERSION="3.3"

# GRAPHIC_DIR="/obs/sphere/GRAPHIC/production/graphic_git/"$GRAPHIC_VERSION"/"
# GRAPHIC_DIR="/home/spectro/cheetham/code/graphic_git/"$GRAPHIC_VERSION"/"
# GRAPHIC_DIR="/Users/cheetham/code/graphic_git/"$GRAPHIC_VERSION"/"

# GRAPHIC_OPTIONS="-mca btl ^openib"
GRAPHIC_OPTIONS=""

GRAPHIC_N_CORES="" #for running on a cluster
GRAPHIC_N_CORES="-n 6" # for running on a laptop. Make sure n_cores is less than or equal to the number of files!

GRAPHIC_N_SKY="3" # the number of cubes to bin together to make a sky

GRAPHIC_N_TARG="$(ls -l Targ/NACO*.fits | wc -l)" # Number of target cubes to bin together.
L_MAX="" #to cut the images down

# The name of the file used to save the output of each process
OUTPUT_FILE="../satpsf_output.txt"

#################
# # NOW WORK ON THE DATA
#################
cd Targ

# Clean up the old txt file
rm $OUTPUT_FILE

# #  Remove the 2 empty pixel rows at the top to have square frames
mpirun $GRAPHIC_N_CORES $GRAPHIC_OPTIONS python $GRAPHIC_DIR"GRAPHIC_rm_overscan_"$GRAPHIC_VERSION".py" --pattern NACO. --fix_naco_columns >> $OUTPUT_FILE

# # Register it with graphic but don't calculate the centre yet. Then make some rough skies without blocking the star out and subtract that from the data.
mpirun $GRAPHIC_N_CORES $GRAPHIC_OPTIONS python $GRAPHIC_DIR"GRAPHIC_naco_register_"$GRAPHIC_VERSION".py" --pattern o_NACO -no_psf -naco >> $OUTPUT_FILE
mpirun $GRAPHIC_N_CORES $GRAPHIC_OPTIONS python $GRAPHIC_DIR"GRAPHIC_n_NaNmedian_"$GRAPHIC_VERSION".py" --pattern o_NACO --num $GRAPHIC_N_SKY --win 0 --info_pattern all_info_ -nofit >> $OUTPUT_FILE
mpirun $GRAPHIC_N_CORES $GRAPHIC_OPTIONS python $GRAPHIC_DIR"GRAPHIC_sky_sub_"$GRAPHIC_VERSION".py" --pattern o_NACO --sky_pattern med --sky_dir sky-num --info_dir cube-info --info_pattern all_info_ -nofit >> $OUTPUT_FILE
# # mpirun $GRAPHIC_N_CORES $GRAPHIC_OPTIONS python $GRAPHIC_DIR"GRAPHIC_pca_sky_sub_"$GRAPHIC_VERSION".py" --pattern o_NACO --sky_pattern med --sky_dir sky-num --info_dir cube-info --info_pattern all_info_*o_NACO -nofit --pca_modes 1 >> $OUTPUT_FILE

# # Delete the rough sky frames (since we don't need them anymore)
rm sky-num/med*.fits

# # Now calculate the actual centre using the rough sky subtracted frames, make a clean sky and subtract that.
mpirun $GRAPHIC_N_CORES $GRAPHIC_OPTIONS python $GRAPHIC_DIR"GRAPHIC_naco_register_"$GRAPHIC_VERSION".py" --pattern nomed_ -naco -saturated >> $OUTPUT_FILE
mpirun $GRAPHIC_N_CORES $GRAPHIC_OPTIONS python $GRAPHIC_DIR"GRAPHIC_n_NaNmedian_"$GRAPHIC_VERSION".py" --pattern o_NACO --num $GRAPHIC_N_SKY --win 40 --info_pattern all_info_*nomed >> $OUTPUT_FILE
rm nomed*.fits # Delete the rough sky subtracted frames from before
mpirun $GRAPHIC_N_CORES $GRAPHIC_OPTIONS python $GRAPHIC_DIR"GRAPHIC_sky_sub_"$GRAPHIC_VERSION".py" --pattern o_NACO --sky_pattern med --sky_dir sky-num --info_dir cube-info --info_pattern all_info_ >> $OUTPUT_FILE
# # mpirun $GRAPHIC_N_CORES $GRAPHIC_OPTIONS python $GRAPHIC_DIR"GRAPHIC_pca_sky_sub_"$GRAPHIC_VERSION".py" --pattern o_NACO --sky_pattern med --sky_dir sky-num --info_dir cube-info --info_pattern all_info_*nomed -nofit --pca_modes 5 >> $OUTPUT_FILE

# # # Maybe we also want to do the star centroid again if it didn't work well on the rough sky subtracted data
rm cube-info/all_info*nomed*.rdb
mpirun $GRAPHIC_N_CORES $GRAPHIC_OPTIONS python $GRAPHIC_DIR"GRAPHIC_naco_register_"$GRAPHIC_VERSION".py" --pattern nomed_ -naco -saturated >> $OUTPUT_FILE

# # # # # Perform frame selection:
# mpirun $GRAPHIC_OPTIONS python $GRAPHIC_DIR"GRAPHIC_frame_selection_"$GRAPHIC_VERSION".py" --pattern nomed_ --info_pattern all_info_*nomed -dithered --centering_nsigma 5 --flux_nsigma 5 --psf_width_nsigma 5 >> $OUTPUT_FILE

# # # # Bin the data into 1 file
mpirun $GRAPHIC_N_CORES $GRAPHIC_OPTIONS python $GRAPHIC_DIR"GRAPHIC_recenter_cubes_"$GRAPHIC_VERSION".py" --pattern nomed --info_pattern "all_info_*nomed*" --naxis3 $GRAPHIC_N_TARG --lmax 512 >> $OUTPUT_FILE

# # # Convert the output to the right format for the GRAPHIC PCA and PHISHI
mpirun $GRAPHIC_N_CORES $GRAPHIC_OPTIONS python $GRAPHIC_DIR"GRAPHIC_convert_for_pca.py" --pattern rc --output_dir "../ADI/" >> $OUTPUT_FILE


cd ..
