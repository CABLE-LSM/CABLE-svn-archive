#!/bin/csh

echo ""
echo "=============Job Setup=================="

echo "USERID =" $USERID                       # Directory - where files are 
echo "DIR    =" $DIR                          # Directory - where files are 
echo "DIRW   =" $DIRW                         # Directory - where files are 
echo "RUNID  =" $RUNID                        # Runid or Jobid from UMUI
echo "YR     =" $YR                           # Number of Years - run or otherwise
#echo "SYR    =" $SYR                          # y/n - Not used yet
#echo "SELYR  =" $SELYR                        # [startyr]/[endyr] - Not used yet
echo "REINIT =" $REINIT                       # 1/3 - For processing when output files are 1/3 months in length
#echo "RESUB =" $RESUB                         # Not used

echo "=============Processing================="

echo "NCI_CP =" $NCI_CP                       #
echo "NCI_ID =" $NCI_ID                       #
echo "CNV2NC =" $CNV2NC                       # convert UM ppfiles to NetCDF before running umplot
#echo "PROC   =" $PROC                         # Not used
echo "Pmonth =" $Pmonth                       # Ext. letter for Files containing Monthly Means
echo "Pdaily =" $Pdaily                       # Ext. letter for Files containing Daily Means & Min/Max
echo "Ptemp1 =" $Ptemp1                       # Ext. letter for Files containing Daily Means & Min/Max
echo "Ptemps =" $Ptemps                       # Ext. letter for Files containing Daily Means & Min/Max
echo "Ptimes =" $Ptimes                       # Ext. letter for Files containing Daily Means & Min/Max
echo "PcasaC =" $PcasaC                       # Ext. letter for Files containing Daily Means & Min/Max

echo "=============Temp Varname==============="

echo "tname  =" $tname                        # temp_`number` - variable name for tscrn
echo "txname =" $txname                       # temp_`number` - variable name for tscrn
echo "tiname =" $tiname                       # temp_`number` - variable name for tscrn

echo "=============Model Setup================"

echo "MODEL  =" $MODEL                        # c/m - cable or moses
echo "RES    =" $RES                          # 48/96 - n48 or n96 resolution
echo "MASK   =" $MASK                         # 1/2 - Which land-sea mask ?
echo "CAL    =" $CAL                          # 360/365/366
echo "CASA   =" $CASA                         # CASA-CNP

echo "=============UMPLOT Setup==============="

echo "UMPLT  =" $UMPLT                        # y/n - do you want to run umplot ?
echo "SPLIT  =" $SPLIT                        # y/n - ONLY if YR is a multiple of 5 and  >= 10
echo "BMRK   =" $BMRK                         # y/n
echo "BDIR   =" $BDIR                         # MOSES: ~ste69f/access/xaiyl_m48_20yrs .OR. ~ste69f/access/xagpb_n96_m_20yrs

echo "=============Format Setup==============="

echo "CHFMT  =" $CHFMT                        # Not used yet
echo "FMT    =" $FMT                          # pdf/eps/jpg - Not used yet

#echo "=============Extra Envars==============="

#echo "USERID =" $USERID                       #
#echo "DIRW   =" $DIRW                         #
#echo "TILE   =" $TILE                         #
#echo "SOIL   =" $SOIL                         #
#echo "FullYrs=" $FullYrs                      #
#echo "MOSES  =" $MOSES                        #
#echo "CABLE  =" $CABLE                        #
#echo "BLOCK  =" $BLOCK                        #
#echo "cabsea =" $cabsea                       #
#echo "cabtem =" $cabtem                       #
#echo "cabhyy =" $cabhyy                       #
#echo "cabtum =" $cabtum                       #
#echo "cabhay =" $cabhay                       #
#echo "cabtha =" $cabtha                       #
#echo "cabman =" $cabman                       #
#echo "cablit =" $cablit                       #
#echo "cabwal =" $cabwal                       #
#echo "cabbon =" $cabbon                       #
#echo "cabhar =" $cabhar                       #
#echo "cablob =" $cablob                       #
#echo "cabvie =" $cabvie                       #
#echo "cabnsa =" $cabnsa                       #
#echo "mossea =" $mossea                       #
#echo "mostem =" $mostem                       #
#echo "moshyy =" $moshyy                       #
#echo "mostum =" $mostum                       #
#echo "moshay =" $moshay                       #
#echo "mostha =" $mostha                       #
#echo "mosman =" $mosman                       #
#echo "moslit =" $moslit                       #
#echo "moswal =" $moswal                       #
#echo "mosbon =" $mosbon                       #
#echo "moshar =" $moshar                       #
#echo "moslob =" $moslob                       #
#echo "mosvie =" $mosvie                       #
#echo "mosnsa =" $mosnsa                       #

echo "========================================"
echo ""

exit

