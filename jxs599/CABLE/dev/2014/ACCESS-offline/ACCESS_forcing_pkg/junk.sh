#!/bin/sh
#PBS -l mem=32Gb
#PBS -l walltime=1200
#PBS -l ncpus=1
#PBS -l jobfs=4GB
#PBS -q normal  
#PBS -j oe
./bins2ncdf_main.py -f part.txt
