#!/bin/csh -x

######################################################################

set CONV2NC = ~/umutils/conv2nc.tcl
set CONV2NCTS = ~dix043/src/python/um/um_timeseries.py
set a = a

######################################################################

    set extlist = 'pm pa pd pf pg ph ps py px'
    foreach ext ( $extlist )
        set filelist = `ls $DIR/$RUNID$a.$ext*[0nbrylgptv] $DIR/$RUNID$a.$ext*dec`
        foreach file ( $filelist )
         if (! -e $file.nc) then
           $CONV2NC -i $file -o $file.nc
         endif
        end
    end

   set extlist = 'pb'
    foreach ext ( $extlist )
        set filelist = `ls $DIR/$RUNID$a.$ext*0`
        foreach file ( $filelist )
         if (! -e $file.nc) then
           $CONV2NC -i $file -o $file.nc
           #python ~dix043/src/python/um/um_fields_subset.py -i $file -o $file.sub -v 3236
           #$CONV2NC -i $file.sub -o $file.nc
           ##$CONV2NC -i $file.sub -o $file.sub.nc
         endif
        end
    end

    set extlist = 'pc pe pi pj'
    foreach ext ( $extlist )
        set filelist = `ls $DIR/$RUNID$a.$ext*0`
        foreach file ( $filelist )
         if (! -e $file.nc) then
           python $CONV2NCTS -i $file -o $file.nc
         endif
        end
    end

######################################################################

exit
