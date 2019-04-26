#!/bin/bash

# This script is for AGPM targets with only a few sky frames
# It is different from the normal reduction in the following ways:
# - Regular sky sub instead of PCA due to low number of sky frames

GRAPHIC_VERSION="3.3"

# Number of processors used to run GRAPHIC
GRAPHIC_N_CORES="" #for running on a cluster
GRAPHIC_N_CORES="-n 6" # for running on a laptop. Make sure n_cores is less than or equal to the number of files!

# Number of sky files
GRAPHIC_N_SKY="$(ls -l Sky/NACO*.fits | wc -l)"
GRAPHIC_N_CORES_SKY="-n $(($((GRAPHIC_N_SKY+1))<${GRAPHIC_N_CORES: -1}?$((GRAPHIC_N_SKY+1)):${GRAPHIC_N_CORES: -1}))" # take the minimum of n_cores and (n_sky+1)
#some of the processes have the rank=0 process doing nothing, so we need n_sky+1 processes to do them all at the same time.

# Number of target cubes to bin together to make the final cube (i.e. all of them)
GRAPHIC_N_TARG="$(ls -l Targ/NACO*.fits | wc -l)"

# The name of the file used to save the output of each process
OUTPUT_FILE="../agpm_output.txt"

##################
# # Make the SKIES
#################
# Create sky frames from the sky directory
cd Sky

# # Clean up the old txt file
rm $OUTPUT_FILE

# # #  Remove the region outside the AGPM 
mpirun $GRAPHIC_N_CORES_SKY python $GRAPHIC_DIR"GRAPHIC_rm_overscan_"$GRAPHIC_VERSION".py" --pattern NACO --l_max 600 --centre_offset 391 425 --fix_naco_columns >> $OUTPUT_FILE

# # #  # Calculate parallactic angle, but don't look for the star centre
mpirun $GRAPHIC_N_CORES_SKY python $GRAPHIC_DIR"GRAPHIC_register_"$GRAPHIC_VERSION".py" --pattern o_NACO -no_psf -naco >> $OUTPUT_FILE

# # # # Calculate the sky frames 
mpirun $GRAPHIC_N_CORES_SKY python $GRAPHIC_DIR"GRAPHIC_n_NaNmedian_"$GRAPHIC_VERSION".py" --pattern o_NACO --num 1 --win 0 --info_pattern all_info_ -nofit >> $OUTPUT_FILE

# Get the AGPM centre from the sky frames
mpirun $GRAPHIC_N_CORES_SKY python $GRAPHIC_DIR"GRAPHIC_naco_agpm_offset_"$GRAPHIC_VERSION".py" --pattern sky-num/med1  >> $OUTPUT_FILE

cd ..

# # clean up the un-needed files
cd Sky/
# rm o_*.fits
cd ..

#################
# # NOW WORK ON THE DATA
#################
cd Targ

#  Remove the region outside the AGPM 
mpirun $GRAPHIC_N_CORES python $GRAPHIC_DIR"GRAPHIC_rm_overscan_"$GRAPHIC_VERSION".py" --pattern NACO --l_max 600 --centre_offset 391 425 --fix_naco_columns >> $OUTPUT_FILE

# # # Register it with graphic but don't calculate the centre yet. Then sky subtract and remove bad pixels
mpirun $GRAPHIC_N_CORES python $GRAPHIC_DIR"GRAPHIC_register_"$GRAPHIC_VERSION".py" --pattern o_NACO -naco -no_psf >> $OUTPUT_FILE
mpirun $GRAPHIC_N_CORES python $GRAPHIC_DIR"GRAPHIC_sky_sub_"$GRAPHIC_VERSION".py" --pattern o_NACO --sky_pattern med --sky_dir ../Sky/sky-num --info_dir cube-info --info_pattern all_info_ -nofit -fix_naco_bad_columns >> $OUTPUT_FILE
# mpirun $GRAPHIC_N_CORES $GRAPHIC_OPTIONS python $GRAPHIC_DIR"GRAPHIC_pca_sky_sub_"$GRAPHIC_VERSION".py" --pattern o_NACO --sky_pattern med --sky_dir ../Sky/sky-num --info_dir cube-info --info_pattern all_info_*o_NACO -nofit --pca_modes 5 -fix_naco_bad_columns >> $OUTPUT_FILE

# # # Now calculate the actual centre using the sky subtracted frames
# # OPTION 1: Fit to the frames before sky subtraction, so you can get the combined star+agpm centre:
# rm cube-info/all_info_3_o_*.rdb # delete the o_ info files and create new ones
# mpirun $GRAPHIC_N_CORES $GRAPHIC_OPTIONS python $GRAPHIC_DIR"GRAPHIC_register_"$GRAPHIC_VERSION".py" --pattern o_ -naco -saturated -remove_striping --search_region 60 >> $OUTPUT_FILE
# # Option 2: Fit to the AGPM centre
# # mpirun $GRAPHIC_N_CORES $GRAPHIC_OPTIONS python $GRAPHIC_DIR"GRAPHIC_register_"$GRAPHIC_VERSION".py" --pattern nomed_ -naco -saturated -agpm_centre --search_region 60 >> $OUTPUT_FILE
# # Option 3: Nico's AGPM fit:
mpirun $GRAPHIC_N_CORES python $GRAPHIC_DIR"GRAPHIC_naco_agpm_register_"$GRAPHIC_VERSION".py" --pattern nomed  >> $OUTPUT_FILE

# # Run frame selection on the data
python $GRAPHIC_DIR"GRAPHIC_frame_selection_"$GRAPHIC_VERSION".py" --pattern nomed_ --info_pattern all_info_*nomed* --centering_nsigma 5 --flux_nsigma 5 --psf_width_nsigma 5 -agpm_centre>> $OUTPUT_FILE
# # python $GRAPHIC_DIR"GRAPHIC_frame_selection_"$GRAPHIC_VERSION".py" --pattern nomed_ --info_pattern all_info_*o_NACO* --centering_nsigma 5 --flux_nsigma 1e6 --psf_width_nsigma 1e6 --debug 3 >> $OUTPUT_FILE

# # # Bin the data (if required)
mpirun $GRAPHIC_N_CORES python $GRAPHIC_DIR"GRAPHIC_recenter_cubes_"$GRAPHIC_VERSION".py" --pattern nomed_ --info_pattern all_info_*f* --naxis3 $GRAPHIC_N_TARG --lmax 600 >> $OUTPUT_FILE

# # # Convert the output to the right format for the GRAPHIC PCA and PHISHI
python $GRAPHIC_DIR"GRAPHIC_convert_for_pca.py" --pattern rc --output_dir "../ADI/" >> $OUTPUT_FILE

cd ..
