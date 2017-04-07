#!/bin/bash
# Use perldoc <script_name> for documentation on each pipeline component.
# AUTHOR: Ben Bai (u5205339@anu.edu.au)
# DATE: Nov 2013-Feb 2014

# Include the config file, which is in valid bash format
sample_name="$1"
CONFIG_FILE="$2"
source "$CONFIG_FILE"

echo ===== Started extract_haplotypes.sh with "[$sample_name]" at $(date) =====

echo vcf2haplofasta "[$sample_name]" at $(date)
perl "$SCRIPT_DIR/pl/vcf2haplofasta.pl" "$sample_name" "$OUT_DIR"

exit
