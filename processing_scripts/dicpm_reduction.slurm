#!/bin/bash

# Reduction script for processing ISPY data with Matthias' DICPM pipeline
# All the actual commands are in the python script, this is just here for consistency
OUTPUT_FILE="../dicpm_output.txt"

cd ADI

# # Clean up the old txt file
rm $OUTPUT_FILE
# Run DICPM
python ~/code/naco_ispy/processing_scripts/dicpm_reduction.py >> $OUTPUT_FILE

cd ..