#!/bin/bash

GRAPHIC_N_CORES="-n 6"

OUTPUT_FILE="../centring_test.txt"
N_TARG="$(ls -l Targ/NACO*.fits | wc -l)"

# The old centring method:
cd Targ
rm cube-info/all_info_*.rdb # delete the o_ info files and create new ones
rm rc_*.fits
mpirun $GRAPHIC_N_CORES python $GRAPHIC_DIR"GRAPHIC_register_"$GRAPHIC_VERSION".py" --pattern o_ -naco -remove_striping --search_region 60 >> $OUTPUT_FILE
mpirun $GRAPHIC_N_CORES python $GRAPHIC_DIR"GRAPHIC_recenter_cubes_"$GRAPHIC_VERSION".py" --pattern nomed_ --info_pattern all_info_ --naxis3 $N_TARG --lmax 600  >> $OUTPUT_FILE
mpirun python $GRAPHIC_DIR"GRAPHIC_convert_for_pca.py" --pattern rc --output_dir "../ADI_old_agpm/"  >> $OUTPUT_FILE
cd ..

# The old double-gauss method:
cd Targ
rm cube-info/all_info_*.rdb # delete the o_ info files and create new ones
rm rc_*.fits
mpirun $GRAPHIC_N_CORES python $GRAPHIC_DIR"GRAPHIC_register_"$GRAPHIC_VERSION".py" --pattern nomed_ -naco -agpm_centre --search_region 60 >> $OUTPUT_FILE
mpirun $GRAPHIC_N_CORES python $GRAPHIC_DIR"GRAPHIC_recenter_cubes_"$GRAPHIC_VERSION".py" --pattern nomed_ --info_pattern all_info_ --naxis3 $N_TARG --lmax 600  >> $OUTPUT_FILE
mpirun python $GRAPHIC_DIR"GRAPHIC_convert_for_pca.py" --pattern rc --output_dir "../ADI_old_gauss/"  >> $OUTPUT_FILE
cd ..


# The new double-gaussian method:
cd Sky
mpirun $GRAPHIC_N_CORES python $GRAPHIC_DIR"GRAPHIC_naco_agpm_offset_"$GRAPHIC_VERSION".py" --pattern sky-num/med1  >> $OUTPUT_FILE
cd ..
cd Targ
rm cube-info/all_info_*.rdb # delete the o_ info files and create new ones
rm rc_*.fits
mpirun $GRAPHIC_N_CORES python $GRAPHIC_DIR"GRAPHIC_naco_agpm_register_"$GRAPHIC_VERSION".py" --pattern nomed  >> $OUTPUT_FILE
mpirun $GRAPHIC_N_CORES python $GRAPHIC_DIR"GRAPHIC_recenter_cubes_"$GRAPHIC_VERSION".py" --pattern nomed_ --info_pattern all_info_ --naxis3 $N_TARG --lmax 600  >> $OUTPUT_FILE
mpirun python $GRAPHIC_DIR"GRAPHIC_convert_for_pca.py" --pattern rc --output_dir "../ADI_new_gauss/"  >> $OUTPUT_FILE
cd ..

# The new double-gaussian method with 1D fitting:
cd Sky
mpirun $GRAPHIC_N_CORES python $GRAPHIC_DIR"GRAPHIC_naco_agpm_offset_"$GRAPHIC_VERSION".py" --pattern sky-num/med1 -fit_1d  >> $OUTPUT_FILE
cd ..
cd Targ
rm cube-info/all_info_*.rdb # delete the o_ info files and create new ones
rm rc_*.fits
mpirun $GRAPHIC_N_CORES python $GRAPHIC_DIR"GRAPHIC_naco_agpm_register_"$GRAPHIC_VERSION".py" --pattern nomed -fit_1d  >> $OUTPUT_FILE
mpirun $GRAPHIC_N_CORES python $GRAPHIC_DIR"GRAPHIC_recenter_cubes_"$GRAPHIC_VERSION".py" --pattern nomed_ --info_pattern all_info_ --naxis3 $N_TARG --lmax 600  >> $OUTPUT_FILE
mpirun python $GRAPHIC_DIR"GRAPHIC_convert_for_pca.py" --pattern rc --output_dir "../ADI_new_gauss/"  >> $OUTPUT_FILE
cd ..


# The new method, centred on the AGPM
cd Targ
rm rc_*.fits
mpirun $GRAPHIC_N_CORES python $GRAPHIC_DIR"GRAPHIC_recenter_cubes_"$GRAPHIC_VERSION".py" --pattern nomed_ --info_pattern all_info_ --naxis3 $N_TARG --lmax 600 -nofit  >> $OUTPUT_FILE
mpirun python $GRAPHIC_DIR"GRAPHIC_convert_for_pca.py" --pattern rc --output_dir "../ADI_new_agpm/"  >> $OUTPUT_FILE
cd ..

