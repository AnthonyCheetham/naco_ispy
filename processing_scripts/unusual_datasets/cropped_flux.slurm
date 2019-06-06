#!/bin/bash
#SBATCH -p daceg -o my.stdout --mail-user=anthony.cheetham@unige.ch --mail-type=ALL --error=error_output_%j

# module add openmpi/gcc/1.8.4
# module add python/2.7.8

# This script is intended to process the flux frames only, and is common for the AGPM and saturated psf observations

GRAPHIC_VERSION="3.3"

# GRAPHIC_DIR="/obs/sphere/GRAPHIC/production/graphic_git/"$GRAPHIC_VERSION"/"
# GRAPHIC_DIR="/home/spectro/cheetham/code/graphic_git/"$GRAPHIC_VERSION"/"
# GRAPHIC_DIR="/Users/cheetham/code/graphic_git/"$GRAPHIC_VERSION"/"

# In case we need to add any command line options to mpirun
GRAPHIC_OPTIONS=""

# Number of processors used to run GRAPHIC
GRAPHIC_N_CORES_FLUX="-n 3"
GRAPHIC_N_SKY_FLUX="3" # Number of cubes to bin together to make a sky

# Number of flux cubes to bin together to make the final cube (i.e. all of them)
GRAPHIC_N_FLUX="$(ls -l Flux/NACO*.fits | wc -l)" 

# The name of the file used to save the output of each process
OUTPUT_FILE="../flux_output.txt"

#################
## Prepare the FLUX frames
#################

cd Flux

# Clean up the old txt file
rm $OUTPUT_FILE

# # # Register it with graphic but don't calculate the centre yet. Then make some rough skies without blocking the star out and subtract that from the data.
mpirun $GRAPHIC_N_CORES_FLUX $GRAPHIC_OPTIONS python $GRAPHIC_DIR"GRAPHIC_register_"$GRAPHIC_VERSION".py" --pattern NACO -naco -no_psf >> $OUTPUT_FILE
mpirun $GRAPHIC_N_CORES_FLUX $GRAPHIC_OPTIONS python $GRAPHIC_DIR"GRAPHIC_n_NaNmedian_"$GRAPHIC_VERSION".py" --pattern NACO --num $GRAPHIC_N_SKY_FLUX --win 0 --info_pattern all_info_ -nofit >> $OUTPUT_FILE
mpirun $GRAPHIC_N_CORES_FLUX $GRAPHIC_OPTIONS python $GRAPHIC_DIR"GRAPHIC_sky_sub_"$GRAPHIC_VERSION".py" --pattern NACO --sky_pattern med --sky_dir sky-num --info_dir cube-info --info_pattern all_info_ -nofit >> $OUTPUT_FILE

# # # Delete the rough sky frames (since we don't need them anymore)
rm sky-num/med*.fits
rm cube-info/all_info_*.rdb

# # Now calculate the actual centre using the rough sky subtracted frames, make a clean sky and subtract that.
mpirun $GRAPHIC_N_CORES_FLUX $GRAPHIC_OPTIONS python $GRAPHIC_DIR"GRAPHIC_register_"$GRAPHIC_VERSION".py" --pattern nomed -naco --psf_width 3 --smooth_width 3  >> $OUTPUT_FILE
mpirun $GRAPHIC_N_CORES_FLUX $GRAPHIC_OPTIONS python $GRAPHIC_DIR"GRAPHIC_n_NaNmedian_"$GRAPHIC_VERSION".py" --pattern NACO --num $GRAPHIC_N_SKY_FLUX --win 10 --info_pattern all_info_*nomed >> $OUTPUT_FILE
rm nomed*.fits # Delete the rough sky frames from before
mpirun $GRAPHIC_N_CORES_FLUX $GRAPHIC_OPTIONS python $GRAPHIC_DIR"GRAPHIC_sky_sub_"$GRAPHIC_VERSION".py" --pattern NACO --sky_pattern med --sky_dir sky-num --info_dir cube-info --info_pattern all_info_ -nofit >> $OUTPUT_FILE

# # Clean bad pixels using the sky frames as reference
mpirun $GRAPHIC_N_CORES_FLUX $GRAPHIC_OPTIONS python $GRAPHIC_DIR"GRAPHIC_bad_pix_"$GRAPHIC_VERSION".py" --pattern nomed --dark_pattern med --dark_dir sky-num --coef 6 >> $OUTPUT_FILE

# # Maybe we also want to do the star centroid again if it didn't work well on the rough sky subtracted data
rm cube-info/*nomed*.rdb
mpirun $GRAPHIC_N_CORES_FLUX $GRAPHIC_OPTIONS python $GRAPHIC_DIR"GRAPHIC_register_"$GRAPHIC_VERSION".py" --pattern nomed -naco --psf_width 3 --smooth_width 3 >> $OUTPUT_FILE

# Perform frame selection:
python $GRAPHIC_DIR"GRAPHIC_frame_selection_"$GRAPHIC_VERSION".py" --pattern nomed --info_pattern all_info_*nomed --centering_nsigma 1e5 --flux_nsigma 5 --psf_width_nsigma 5 --debug 3 -dithered >> $OUTPUT_FILE

# # Bin the data into 1 file
# mpirun $GRAPHIC_N_CORES_FLUX $GRAPHIC_OPTIONS python $GRAPHIC_DIR"GRAPHIC_recenter_cubes_"$GRAPHIC_VERSION".py" --pattern nomed --info_pattern all_info_framesel*nomed --naxis3 $GRAPHIC_N_FLUX -combine_frames --lmax 40 >> $OUTPUT_FILE
mpirun $GRAPHIC_N_CORES_FLUX $GRAPHIC_OPTIONS python $GRAPHIC_DIR"GRAPHIC_recenter_cubes_"$GRAPHIC_VERSION".py" --pattern nomed --info_pattern all_info_framesel*nomed --naxis3 $GRAPHIC_N_FLUX --lmax 40 >> $OUTPUT_FILE

# # # Convert the output to the right format for the GRAPHIC PCA and PHISHI
python $GRAPHIC_DIR"GRAPHIC_convert_for_pca.py" --pattern rc -collapse_cube --output_dir "../ADI/" -skip_parang --output_file 'flux.fits' >> $OUTPUT_FILE
python $GRAPHIC_DIR"GRAPHIC_convert_for_pca.py" --pattern rc --output_dir "../ADI/" -skip_parang --output_file 'flux_cube.fits' >> $OUTPUT_FILE

cd ..

# # clean up the un-needed files
cd Flux/
rm cl_*.fits
rm nomed_*.fits
cd ..
