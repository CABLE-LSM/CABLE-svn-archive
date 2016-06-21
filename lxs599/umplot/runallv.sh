#!/bin/csh 

######################################################################

if ( $HOSTNAME == ruby ) then
 set CONV2NC = ~ste69f/umutils/conv2nc.tcl
 set CONV2NCTS = ~dix043/src/python/um/um_timeseries.py
else
 set CONV2NC = ~lxs599/umutils/conv2nc.tcl
 set CONV2NCTS = ~mrd599/src/python/um/um_timeseries.py
endif
set a = a

######################################################################

    set extlist = 'pm pa pe pb pc pi pj pd pf pg ph ps py px'
    foreach ext ( $extlist )
        if ($VN == 85) then
        set filelist = `ls $DIR/$RUNID$a.$ext???????`
        else
        set filelist = `ls $DIR/$RUNID$a.$ext?????`
        endif
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
