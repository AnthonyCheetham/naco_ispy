#!/bin/bash
#SBATCH -p daceg -o my.stdout --mail-user=anthony.cheetham@unige.ch --mail-type=ALL --error=error_output_%j

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

# GRAPHIC_VERSION="3.3"
# GRAPHIC_OPTIONS="-mca btl ^openib"
GRAPHIC_OPTIONS=""

# GRAPHIC_DIR="/obs/sphere/GRAPHIC/production/graphic_git/"$GRAPHIC_VERSION"/"
# GRAPHIC_DIR="/home/spectro/cheetham/code/graphic_git/"$GRAPHIC_VERSION"/"
# GRAPHIC_DIR="/Users/cheetham/code/graphic_git/"$GRAPHIC_VERSION"/"

GRAPHIC_N_CORES="" #for running on a cluster
GRAPHIC_N_CORES="-n 6" # for running on a laptop. Make sure n_cores is less than or equal to the number of files!

GRAPHIC_FLAT_OUTPUT_SIZES="--output_sizes 776 768 600 600 256 264"
GRAPHIC_FLAT_CENTRE_OFFSET="--centre_offset 0 0"

#################
## Prepare the FLAT frames
#################

# cd Flats

Process the flats
mpirun $GRAPHIC_N_CORES $GRAPHIC_OPTIONS python $GRAPHIC_DIR"GRAPHIC_masterflat_naco_"$GRAPHIC_VERSION".py" --pattern NACO --t_max 1e3 -fix_cols --std_max 0.2 -only_lprime $GRAPHIC_FLAT_OUTPUT_SIZES $GRAPHIC_FLAT_CENTRE_OFFSET #>> ../flats_output.txt

cd ..

#################
## Prepare the Astrometric Calibration files
#################

cd AstCal

python $GRAPHIC_DIR"GRAPHIC_astrometric_calib_naco_"$GRAPHIC_VERSION".py" --pattern NACO --t_max 60 -fix_cols -derotate -pca_sky_subtract -find_flat -window -cosmic_ray >> ../astcal_output.txt

cd ..

