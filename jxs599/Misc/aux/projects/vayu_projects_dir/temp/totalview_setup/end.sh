#!/bin/ksh
 set -x

a=a
c=:
save=/short/p66/jxs599/save
conv2nc=$HOME/umutils/conv2nc.tcl

echo "ARCHIVE START" `date`

cd $UM_DATAM

rid=$RUNID

#for ext in d pa pb pc pd pe pf pg ph pi pj pm py ps; do
#  for ffile in $rid$a.$ext*; do
#    sleep=10
#    for k in 1 2 3 4; do
#      if [[ -s $ffile ]]; then
#        break
#      else
#        echo "Sleeping " $sleep
#        sleep $sleep
#        (( sleep = sleep * 2 ))
#      fi
#    done
#    if [[ -s $ffile ]]; then
#
#      if [[ ! -d $save/$rid ]]; then
#        rm -f $save/$rid
#        mkdir $save/$rid
#      fi
#
#      if [[ ($ext != d ) && ( $ext != ph) && ($ext != pc) && ($ext != pe)  && ($ext != pf)  && ($ext != pg)  && ($ext != pi)  && ($ext != pj) && ($ext != pd) ]]; then
#        $conv2nc -i $ffile -o $ffile.nc
#        mv $ffile.nc $save/$rid
#      else
#        cp $ffile $save/$rid
#      fi
#    fi
#  done
#
#  llist=`ls -1t $rid$a.$ext* | tail -n +3`
## for ffile in $llist; do
##   rm -f $ffile
## done
#done
#
exit
