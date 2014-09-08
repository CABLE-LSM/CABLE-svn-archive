#!/bin/sh
#PBS -l mem=32Gb
#PBS -l walltime=4800
#PBS -l ncpus=1
#PBS -l jobfs=4GB
#PBS -q normal  
#PBS -j oe

cd /short/p66/jxs599/ACCESS-offline-instance/ACCESS_forcing_pkg_update.py/
./bins2ncdf_main.py -f part.txt
