#!/bin/bash

# module add openmpi/gcc/1.8.4
# module add python/2.7.8

GRAPHIC_VERSION="3.3"

# GRAPHIC_DIR="/obs/sphere/GRAPHIC/production/graphic_git/"$GRAPHIC_VERSION"/"
# GRAPHIC_DIR="/home/spectro/cheetham/code/graphic_git/"$GRAPHIC_VERSION"/"
# GRAPHIC_DIR="/Users/cheetham/code/graphic_git/"$GRAPHIC_VERSION"/"

GRAPHIC_OPTIONS=""

# The name of the file used to save the output of each process
OUTPUT_FILE="../pca_output.txt"
# OUTPUT_FILE="pca_output.txt"


##################
# # Run PCA
#################

cd ADI

# # Clean up the old txt file
rm $OUTPUT_FILE

# # Run PCA
#   PSF is ~5 pix FWHM, and using a 414x414 input file with 2 fwhm per annulus. r_min = 3 pix
#   35 annuli is ok
#   We really only need 100 pix radius, so r_min = 3, r_max=103, n_annuli=10 should be fine
#   arc length of 1e4 means it wont split the annuli up azimuthally.
python $GRAPHIC_DIR"GRAPHIC_pca.py" --output_dir "./GRAPHIC_PCA/" --threads 6 --arc_length 1e4 --n_annuli 6 --n_fwhm 1.25 --n_modes 0.45 --min_reference_frames 0.5 --fwhm 5. --r_min 3. --r_max 63. >> $OUTPUT_FILE
# For a small FoV, use r_min = 3, r_max = 51, n_annuli = 6
# python $GRAPHIC_DIR"GRAPHIC_pca.py" --output_dir "./GRAPHIC_PCA/" --threads 6 --arc_length 1e4 --n_annuli 6 -median_combine --n_fwhm 0.75 --n_modes 12 --min_reference_frames 15 --fwhm 4. --r_min 3. --r_max 51 >> $OUTPUT_FILE

# Now run the detection limits and self-subtraction code
python $GRAPHIC_DIR"GRAPHIC_contrast_pca_3.3.py" --output_dir "./GRAPHIC_PCA/" --threads 6 --plate_scale 0.027 --median_filter_length 25 --smooth_image_length 1.75 --r_max 63. --n_throughput 10 >> $OUTPUT_FILE
# python $GRAPHIC_DIR"GRAPHIC_contrast_pca_3.3.py" --output_dir "./GRAPHIC_PCA/" --threads 6 --plate_scale 0.027 --median_filter_length 20 --smooth_image_length 1.25 --r_max 51. >> $OUTPUT_FILE

cd ..

