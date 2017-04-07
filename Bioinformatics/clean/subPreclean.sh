#!/bin/bash
#$ -cwd
#$ -j y
#$ -N preclean 

perl /moritz/jgb/readsFeb2016/LN01/scripts/pre-cleanup.pl
