#!/bin/bash

# The flux frames for the 2013 HD 169142 dataset were taken with separate sky frames, and not dithered.
# So they need to be processed like the AGPM frames

GRAPHIC_VERSION="3.3"

# In case we need to add any command line options to mpirun
GRAPHIC_OPTIONS=""

# Number of processors used to run GRAPHIC
GRAPHIC_N_CORES="" #for running on a cluster
GRAPHIC_N_CORES="-n 6" # for running on a laptop. Make sure n_cores is less than or equal to the number of files!

# Number of sky files
GRAPHIC_N_SKY="$(ls -l Flux_Sky/NACO*.fits | wc -l)"
GRAPHIC_N_CORES_SKY="-n $(($GRAPHIC_N_SKY<${GRAPHIC_N_CORES: -1}?$GRAPHIC_N_SKY:${GRAPHIC_N_CORES: -1}))" # take the minimum of n_cores and n_sky

# Number of flux cubes to bin together to make the final cube (i.e. all of them)
GRAPHIC_N_FLUX="$(ls -l Flux/NACO*.fits | wc -l)" 

# The name of the file used to save the output of each process
OUTPUT_FILE="../flux_output.txt"

##################
# # Make the SKIES
#################
# Create sky frames from the sky directory
cd Flux_Sky

# # Clean up the old txt file
# rm $OUTPUT_FILE

# # # #  # Calculate parallactic angle, but don't look for the star centre
# mpirun $GRAPHIC_N_CORES_SKY $GRAPHIC_OPTIONS python $GRAPHIC_DIR"GRAPHIC_naco_register_"$GRAPHIC_VERSION".py" --pattern NACO -no_psf -naco >> $OUTPUT_FILE

# # # # # Calculate the sky frames 
# mpirun $GRAPHIC_N_CORES_SKY $GRAPHIC_OPTIONS python $GRAPHIC_DIR"GRAPHIC_n_NaNmedian_"$GRAPHIC_VERSION".py" --pattern NACO --num 1 --win 0 --info_pattern all_info_ -nofit >> $OUTPUT_FILE

cd ..

# # clean up the un-needed files
cd Flux_Sky/
# rm o_*.fits
cd ..

#################
# # NOW WORK ON THE DATA
#################
cd Flux

# # # # Register it with graphic but don't calculate the centre yet. Then sky subtract and remove bad pixels
# mpirun $GRAPHIC_N_CORES $GRAPHIC_OPTIONS python $GRAPHIC_DIR"GRAPHIC_naco_register_"$GRAPHIC_VERSION".py" --pattern NACO -naco -no_psf >> $OUTPUT_FILE
# mpirun $GRAPHIC_N_CORES $GRAPHIC_OPTIONS python $GRAPHIC_DIR"GRAPHIC_pca_sky_sub_"$GRAPHIC_VERSION".py" --pattern NACO --sky_pattern med --sky_dir ../Flux_Sky/sky-num --info_dir cube-info --info_pattern all_info_*NACO -nofit --pca_modes 3 >> $OUTPUT_FILE

# # # Now calculate the actual centre using the sky subtracted frames
# mpirun $GRAPHIC_N_CORES_FLUX $GRAPHIC_OPTIONS python $GRAPHIC_DIR"GRAPHIC_naco_register_"$GRAPHIC_VERSION".py" --pattern nomed -naco --psf_width 3 --smooth_width 3 >> $OUTPUT_FILE

# # Run frame selection on the data
# mpirun $GRAPHIC_OPTIONS python $GRAPHIC_DIR"GRAPHIC_frame_selection_"$GRAPHIC_VERSION".py" --pattern nomed --info_pattern all_info_*nomed --centering_nsigma 1e5 --flux_nsigma 5 --psf_width_nsigma 5 --debug 3 -dithered >> $OUTPUT_FILE

# # Bin the data into 1 file
# mpirun $GRAPHIC_N_CORES_FLUX $GRAPHIC_OPTIONS python $GRAPHIC_DIR"GRAPHIC_recenter_cubes_"$GRAPHIC_VERSION".py" --pattern nomed --info_pattern all_info_framesel*nomed --naxis3 $GRAPHIC_N_FLUX -combine_frames --lmax 40 >> $OUTPUT_FILE

# And copy the flux file to the 
# # # Convert the output to the right format for the GRAPHIC PCA and PHISHI
mpirun $GRAPHIC_N_CORES $GRAPHIC_OPTIONS python $GRAPHIC_DIR"GRAPHIC_convert_for_pca.py" --pattern rc --output_dir "../ADI/" -skip_parang --output_file 'flux.fits' >> $OUTPUT_FILE

cd ..
