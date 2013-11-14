#!/bin/csh

# Raijin
# This script moves data from save directory back to cherax

#####################################################################
#####################################################################

set saveid = $NCI_ID
#set saveid = lxs599
set SAVEHOST = raijin.nci.org.au
#set SAVEHOST = raijin1.nci.org.au
set LOCALSAVE = ~/access
set SAVEBASE = /short/p66/$saveid
#set SAVEBASE = /short/p66/$saveid/save

set c = :
set a = @
set A = a
set t = "~"
set date = `date`
@ nom = ${YR} * 12

set pext = (pa pb pc pd pe pf pg ph pi pj pm pt)
#set pext2 = (ps px py) # 4seas, decadal, yearly

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
      echo "(0)     NOTE: newer files might replace older files in" $LOCALSAVE/$topdir
      echo " "
     endif
     #rsync -avu -e ssh $saveid$a$SAVEHOST$c"$SAVEBASE/$topdir/$topdir$A.p*" $LOCALSAVE/$topdir/
     foreach ext ($pext)
      set filelist=`ssh -l $saveid $SAVEHOST "ls $SAVEBASE/$topdir/$topdir$A.$ext* | head -${nom}"`
      if ( ${#filelist} > 0 ) then
       rsync -avu -e ssh $saveid$a$SAVEHOST$c"$filelist" $LOCALSAVE/$topdir/
      endif
     end
    #endif

  end # foreach topdir

else  # png
  echo "$date - $SAVEHOST is not responding"
  exit(1)
endif # png

######################################################################
######################################################################  

exit

