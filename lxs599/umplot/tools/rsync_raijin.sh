#!/bin/csh

# Raijin
# This script moves data from save directory back to cherax

#####################################################################
#####################################################################

set saveid = $NCI_ID
#set saveid = lxs599
set SAVEHOST = raijin.nci.org.au
#set SAVEHOST = r-dm.nci.org.au
set LOCALSAVE = ~/access
set SAVEBASE = /short/p66/$saveid
#set SAVEBASE = /short/p66/$saveid/save

set c = :
set a = @
set A = a
set t = "~"
set date = `date`

#####################################################################  
# MOVE FORECAST DATA TO DATASTORE 
#####################################################################  

set png = `ssh -l $saveid $SAVEHOST uptime | grep -wo load`
if ( $png != "" ) then
  #set topdirlist=`ssh -l $saveid $SAVEHOST "cd $SAVEBASE; ls -d *"`
  set topdirlist=$RUNID

  foreach topdir ($topdirlist)

    echo "$date - rsync $saveid$a$SAVEHOST$c$SAVEBASE/$topdir"
    if (! -e $LOCALSAVE) mkdir $LOCALSAVE
    #cd $LOCALSAVE
    #if (! -e $topdir) mkdir $topdir
    #cd $topdir

    ##rsync -avu -e ssh $saveid$a$SAVEHOST$c"$SAVEBASE/$topdir/$topdir$A.p[acef]h8*" $topdir/
    #if (-e $DIR) then
    # echo "NOTE: newer files might replace older files in" $DIR
    # rsync -avu -e ssh $saveid$a$SAVEHOST$c"$SAVEBASE/$topdir/$topdir$A.p*" $DIR/
    #else
     if (! -e $LOCALSAVE/$topdir) then
      mkdir $LOCALSAVE/$topdir
     else
      echo " "
      echo "NOTE: newer files might replace older files in" $LOCALSAVE/$topdir
      echo " "
     endif
     rsync -avu -e ssh $saveid$a$SAVEHOST$c"$SAVEBASE/$topdir/$topdir$A.p*" $LOCALSAVE/$topdir/
    #endif

  end # foreach topdir

else  # png
  echo "$date - $SAVEHOST is not responding"
  exit(1)
endif # png

######################################################################
######################################################################  

exit

