#!/bin/bash
# Use perldoc <script_name> for documentation on each pipeline component.
# AUTHOR: Ben Bai (u5205339@anu.edu.au)
# DATE: Nov 2013-Feb 2014

# Include the config file, which is in valid bash format
sample_name="$1"
CONFIG_FILE="$2"
source "$CONFIG_FILE"

echo ===== Started cleaning library "[$sample_name]" at $(date) =====

echo clean library "[$sample_name]" at $(date)

CLEAN_DIR="${OUT_DIR}${sample_name}/"
rm $CLEAN_DIR*.fasta
rm $CLEAN_DIR*.blast
rm -r $CLEAN_DIR*/*_call_velvet_assemblies/*_k*
rm -r $CLEAN_DIR*/*_assemble_by_prot
rm -r $CLEAN_DIR*/*_catcontigs
rm -r $CLEAN_DIR/*_mapsnp

### zip
cd $CLEAN_DIR
DIRS=($(ls -d ENSSHAP*))
tar czf ${sample_name}_assemblyfiles.tar.gz ${DIRS[@]}
rm -r ${DIRS[@]}

exit
