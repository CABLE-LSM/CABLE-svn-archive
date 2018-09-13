#!/bin/csh

######################################################################

# Ruby - ste69f, NCI - lxs599 ----------------------
# Ruby - dix043, NCI - mrd599 ----------------------
if ( $HOSTNAME == ruby ) then
 set CONV2NC   = ~ste69f/umutils/conv2nc.tcl
 #set CONV2NC4  = 
 set CONV2NCTS = ~dix043/src/python/um/um_timeseries.py
 set FLDSUBSET = ~dix043/src/python/um/um_fields_subset.py
else
 set CONV2NC   = ~lxs599/umutils/conv2nc.tcl
 set CONV2NC4  = ~access/bin/um2netcdf4.py
 set CONV2NCTS = ~mrd599/src/python/um/um_timeseries.py
 set FLDSUBSET = ~mrd599/src/python/um/um_fields_subset.py
endif
set a = a
@ nom = ($FullYrs * 12) # number of files for no. of years and months
@ nod = ($FullYrs * 12 * 31) # number of files for no. of years, months and days

######################################################################

#    cd $DIR
    #set extlist = 'pm pa pe pb pc pd pf pg ph ps py px'
    set extlist = 'pm pa pe pb pc pi pj pd pf pg ph ps py px'
    foreach ext ( $extlist )
         
      if ( $CPL == y ) then
       set fillist = `ls $RUNID.$ext-?????????? | head -$nom`
       foreach file ( $fillist )
        if (! -e $file.nc) then
         echo "Converting " $file " to " $file".nc"
         $CONV2NC -i $file -o $file.nc
         #python $CONV2NC4 -i $file -o $file.nc
        else
         echo "File" $file "Already Exists!"
        endif
       end # fe file
      else # CPL
       set i = 1
       if ($ext == $Ptimes && $RES == 320) then
       set fillist = `ls $RUNID$a.$ext????? | head -$nod`
       else
       if ($VN >= 85) then
       set fillist = `ls $RUNID$a.$ext??????? | head -$nom`
       else
       set fillist = `ls $RUNID$a.$ext????? | head -$nom`
       endif
       endif
       #set fillist = `ls $DIR/$RUNID$a.$ext?????`
       set cfiles = ${#fillist}

       if ($cfiles > 0) then
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
           #$CONV2NC -i $file -o $newlist[$i]_all.nc
           python $FLDSUBSET -i $file -o $file.sub2 -x 1201 -x 1235 -x 2201 -x 2207 -x 0238 -x 2204 -x 2407
           $CONV2NC -i $file.sub2 -o $newlist[$i]_noswlw.nc
           python $FLDSUBSET -i $file -o $file.sub3 -v 1201 -v 1235 -v 2201 -v 2207 -v 0238 -v 2204 -v 2407
           $CONV2NC -i $file.sub3 -o $newlist[$i]_swlw.nc
           #set cdate=`cdo showdate $newlist[$i]_swlw.nc`
           #set ctime=`cdo showtime $newlist[$i]_swlw.nc`
           #cdo inttime,$cdate[1],$ctime[1],30minutes $newlist[$i]_swlw.nc $newlist[$i]_swlw_intp30min.nc
           python $FLDSUBSET -i $file -o $file.sub -v 3234 -v 3236 -v 3217 -v 5226 -v 3333
           #python $CONV2NCTS -i $file.sub -o $newlist[$i].nc
           $CONV2NC -i $file.sub -o $newlist[$i].nc
           rm $file.sub $file.sub2 $file.sub3
          else # Ptimes
           if ($ext == $Pdaily && $Ptemp1 != $Ptemps ) then
            #$CONV2NC -i $file -o $newlist[$i]_all.nc
            if ($VN > 85) then
            python $FLDSUBSET -i $file -o $file.sub
            else
            python $FLDSUBSET -i $file -o $file.sub -v 3236
            endif
            $CONV2NC -i $file.sub -o $newlist[$i].nc
            #python $CONV2NC4 -i $file.sub -o $newlist[$i].nc
            rm $file.sub
           else # Pdaily
           if ( $TSTEP == 288 && $ext == pa ) then
            python $CONV2NC4 -i $file -o $newlist[$i].nc
           else
            $CONV2NC -i $file -o $newlist[$i].nc
            #python $CONV2NC4 -i $file -o $newlist[$i].nc
           endif
           endif # Pdaily
          endif # Ptimes
         else   # newlist
          echo "File" $newlist[$i] "Already Exists!"
         endif  # newlist
         set ntime=`cdo ntime $newlist[$i].nc`
         if ($VN == 73) then
         if (($ext == pc && $ntime >= 224)) then # 224 = 8ts * 1month
          if (! -e $newlist[$i]_swlw.nc) then
           #$CONV2NC -i $file -o $newlist[$i]_all.nc
           python $FLDSUBSET -i $file -o $file.sub2 -x 1201 -x 1235 -x 2201 -x 2207 -x 0238 -x 2204 -x 2407
           $CONV2NC -i $file.sub2 -o $newlist[$i]_noswlw.nc
           python $FLDSUBSET -i $file -o $file.sub3 -v 1201 -v 1235 -v 2201 -v 2207 -v 0238 -v 2204 -v 2407
           $CONV2NC -i $file.sub3 -o $newlist[$i]_swlw.nc
           #python $FLDSUBSET -i $file -o $file.sub -v 3234 -v 3236 -v 3217 -v 5226 -v 3333
           #$CONV2NC -i $file.sub -o $newlist[$i].nc
           rm $file.sub2 $file.sub3 #$file.sub
          endif
         endif
         if (($ext == pg || $ext == ph)) then
          if (! -e $newlist[$i]_swlw.nc) then
           #$CONV2NC -i $file -o $newlist[$i]_all.nc
           python $FLDSUBSET -i $file -o $file.sub2 -x 1201 -x 1208 -x 1209 -x 1210 -x 1211 -x 2204 -x 2205 -x 2206 -x 2207 -x 2208
           $CONV2NC -i $file.sub2 -o $newlist[$i]_noswlw.nc
           python $FLDSUBSET -i $file -o $file.sub3 -v 1201 -v 1208 -v 1209 -v 1210 -v 1211 -v 2204 -v 2205 -v 2206 -v 2207 -v 2208
           $CONV2NC -i $file.sub3 -o $newlist[$i]_swlw.nc
           rm $file.sub2 $file.sub3
          endif
         endif
         endif
         @ i++
        end # fe file

       endif # cfiles
      endif # CPL

    end # foreach
#    cd $DIRW

######################################################################

exit

