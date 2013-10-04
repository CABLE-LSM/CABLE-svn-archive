#!/bin/csh
set cwd = `pwd`
set meee = jxs599 
set cable_shine = /short/p66/jxs599/cable_shine 

rsync -avrz $cwd jxs599@vayu.nci.org.au:$cable_shine 

ssh vayu.nci.org.au -l jxs599 -n /opt/pbs/bin/qsub -S /bin/csh $cable_shine/bCABLE-1.9/offline/cable.sh 

