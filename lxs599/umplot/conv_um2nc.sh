#!/bin/csh 

######################################################################

# Cherax - ste69f, NCI - lxs599 ----------------------
# Cherax - dix043, NCI - mrd599 ----------------------
if ( $HOSTNAME == cherax ) then
 set CONV2NC   = ~ste69f/umutils/conv2nc.tcl
 set CONV2NCTS = ~dix043/src/python/um/um_timeseries.py
 set FLDSUBSET = ~dix043/src/python/um/um_fields_subset.py
else
 set CONV2NC   = ~lxs599/umutils/conv2nc.tcl
 set CONV2NCTS = ~mrd599/src/python/um/um_timeseries.py
 set FLDSUBSET = ~mrd599/src/python/um/um_fields_subset.py
endif
set a = a

######################################################################

#    cd $DIR
    set extlist = 'pm pa pe pb pc pi pj pd pf pg ph ps py px'
    foreach ext ( $extlist )
         
      if ( $CPL == y ) then
       set fillist = `ls $RUNID.$ext-??????????`
       foreach file ( $fillist )
        if (! -e $file.nc) then
         $CONV2NC -i $file -o $file.nc
        else
         echo "File" $file "Already Exists!"
        endif
       end # fe file
      else # CPL
       set i = 1
       set fillist = `ls $RUNID$a.$ext?????`
       #set fillist = `ls $DIR/$RUNID$a.$ext?????`

        # change months to numbers in .p
        set newlist=`echo $fillist | sed -e 's/jan/001/g'`
        set newlist=`echo $newlist | sed -e 's/feb/002/g'`
        set newlist=`echo $newlist | sed -e 's/mar/003/g'`
        set newlist=`echo $newlist | sed -e 's/apr/004/g'`
        set newlist=`echo $newlist | sed -e 's/may/005/g'`
        set newlist=`echo $newlist | sed -e 's/jun/006/g'`
        set newlist=`echo $newlist | sed -e 's/jul/007/g'`
        set newlist=`echo $newlist | sed -e 's/aug/008/g'`
        set newlist=`echo $newlist | sed -e 's/sep/009/g'`
        set newlist=`echo $newlist | sed -e 's/oct/010/g'`
        set newlist=`echo $newlist | sed -e 's/nov/011/g'`
        set newlist=`echo $newlist | sed -e 's/dec/012/g'`

        foreach file ( $fillist )
         if (! -e $newlist[$i].nc) then
          echo "Converting " $file " to " $newlist[$i]".nc"
          if ($ext == $Ptimes) then
           python $FLDSUBSET -i $file -o $file.sub2 -x 1201 -x 1235 -x 2201 -x 2207
           $CONV2NC -i $file.sub2 -o $newlist[$i]_noswlw.nc
           python $FLDSUBSET -i $file -o $file.sub3 -v 1201 -v 1235 -v 2201 -v 2207
           $CONV2NC -i $file.sub3 -o $newlist[$i]_swlw.nc
           python $FLDSUBSET -i $file -o $file.sub -v 3234 -v 3236 -v 3217 -v 5226 -v 3333
           python $CONV2NCTS -i $file.sub -o $newlist[$i].nc
           rm $file.sub
           rm $file.sub2
           rm $file.sub3
          else
           $CONV2NC -i $file -o $newlist[$i].nc
          endif
         else
          echo "File" $newlist[$i] "Already Exists!"
         endif
         @ i++
        end # fe file

      endif # CPL

    end # foreach
#    cd $DIRW

######################################################################

exit

