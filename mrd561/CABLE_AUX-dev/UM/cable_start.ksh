#!/bin/ksh

# As scripts are in ksh, and we need to use ksh env vars, we need to  
# SOURCE ~/.kshrc
. ~/.kshrc


# incase this is not done already, setup necessary links in CABLE-AUX

# IF dir CABLE-AUX doesn't already exist in your $HOME ($CABLE_AUX) directory, ... 
if [[ ! -d $CABLE_AUX/CABLE-AUX ]]; then

   # make it  
   /bin/mkdir $CABLE_AUX/CABLE-AUX

   # cp default CABLE-AUX which copies files ands sets necessary links
   /bin/cp -r /projects/access/CABLE-AUX/UM $CABLE_AUX/CABLE-AUX
   /bin/cp -r /projects/access/CABLE-AUX/core $CABLE_AUX/CABLE-AUX

fi


# A dir /short/p66/`whoami`/$RUNID is generally created as a RUN dir for your UMUI job,
# however in this case we need to create the dir to do some work prior to UM execution
# IF this dir doesn't already exist, ... 
if [[ ! -d /short/w35/`whoami`/$RUNID ]]; then
   
   # make it  
   /bin/mkdir /short/w35/`whoami`/$RUNID

fi


# This is to ensure that each run in the same experiment uses the same cable.nml 

if [[ $TYPE == 'NRUN' ]]; then

   /bin/cp $CABLE_AUX/CABLE-AUX/UM/cable.nml /short/w35/`whoami`/$RUNID/.
   /bin/cp /short/w35/`whoami`/$RUNID/cable.nml /short/w35/`whoami`/$RUNID/cable.nml.`date +%d.%m.%y`

fi

if [[ ! -f /short/w35/`whoami`/$RUNID/libcable-2.0.a ]]; then
   /bin/cp $CABLE_AUX/CABLE-AUX/UM/libcable-2.0.a /short/w35/`whoami`/$RUNID/. 
fi

if [[ -f /short/w35/`whoami`/$RUNID/libcable-2.0.a ]]; then
   /bin/cp $CABLE_AUX/CABLE-AUX/UM/libcable-2.0.a /short/w35/`whoami`/$RUNID
   /bin/cp /short/w35/`whoami`/$RUNID/libcable-2.0.a /short/w35/`whoami`/$RUNID/libcable-2.0.a.`date +%d.%m.%y`
fi


