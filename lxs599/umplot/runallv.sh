#!/bin/csh 

######################################################################

set CONV2NC = ~ste69f/umutils/conv2nc.tcl
set CONV2NCTS = ~dix043/src/python/um/um_timeseries.py
set a = a

######################################################################

    set extlist = 'pm pa pe pb pc pi pj pd pf pg ph ps py px'
    foreach ext ( $extlist )
        set filelist = `ls $DIR/$RUNID$a.$ext?????`
        foreach file ( $filelist )
         if (! -e $file.nc) then
           $CONV2NC -i $file -o $file.nc
         else
           echo "File" $file.nc "Already Exists!"
         endif
        end
    end

#   set filelist = `ls $DIR/$RUNID$a.$ext?????`
#   #python ~dix043/src/python/um/um_fields_subset.py -i $file -o $file.sub -v 3236
#   #$CONV2NC -i $file.sub -o $file.nc
#   ##$CONV2NC -i $file.sub -o $file.sub.nc
#   python $CONV2NCTS -i $file -o $file.nc

######################################################################

exit
