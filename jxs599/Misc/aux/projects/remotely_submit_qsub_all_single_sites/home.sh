#/bin/csh
set cable_shine = /short/p66/jxs599/cable_shine 
rsync -avrz jxs599@vayu.nci.org.au:$cable_shine/bCABLE-1.9/offline/test_offline offline/
