#!/bin/bash
#$ -cwd
#$ -r y
#$ -N clean_rw 
#$ -l virtual_free=14G,h_vmem=16G
#$ -t 1-56:1
#$ -tc 1 
###$ -q bigmem.q

sleep 15
echo $SGE_TASK_ID

libs=( $(cat /moritz/jgb/readsFeb2016/LN01/liblist.LN01 ) )

filnum=$[SGE_TASK_ID-1]
parlib=${libs[$filnum]}

#parlib=SP07_indexing1
echo $parlib

## http://users.sdsc.edu/~glockwood/comp/installing.php
source ~/python27-jgb/bin/activate


#parlib=SP07_indexing1

od=/moritz/jgb/readsFeb2016/LN01/clean/${parlib}
mkdir $od

f=/moritz/jgb/readsFeb2016/LN01/preclean/
o=$od
a=/moritz/jgb/cleantools/Adapters.fa
b=/moritz/jgb/readsFeb2016/LN01/libinfo.LN01.txt
t=/home/jgb/software/trimmomatic/Trimmomatic-0.22/trimmomatic-0.22.jar
c=/moritz/jgb/cleantools/contaminants.fa
e=200
m=40
i=D9S
d=$parlib
g=0.1
l=100
n=0.6
r=0.5
h=36


nohup perl /moritz/jgb/readsFeb2016/LN01/scripts/scrubReads.pl -f $f -o $o -a $a -b $b -t $t -c $c -e $e -m $m -i $i -d $d -g $g -l $l -n $n -r $r -h $h 1> ${parlib}.1 2> ${parlib}.2
