#!/bin/csh

echo ""
echo "=============Job Setup=================="

if (${?USERID} == 1) then
echo "USERID =" $USERID                       # Directory - where files are 
endif
if (${?DIR} == 1) then
echo "DIR    =" $DIR                          # Directory - where files are 
endif
if (${?DIRW} == 1) then
echo "DIRW   =" $DIRW                         # Directory - where files are 
endif
if (${?RUNID} == 1) then
echo "RUNID  =" $RUNID                        # Runid or Jobid from UMUI
endif
if (${?YR} == 1) then
echo "YR     =" $YR                           # Number of Years - run or otherwise
endif
if (${?SYR} == 1) then
#echo "SYR    =" $SYR                          # y/n - Not used yet
endif
if (${?SELYR} == 1) then
#echo "SELYR  =" $SELYR                        # [startyr]/[endyr] - Not used yet
endif
if (${?REINIT} == 1) then
echo "REINIT =" $REINIT                       # 1/3 - For processing when output files are 1/3 months in length
endif
if (${?RESUB} == 1) then
#echo "RESUB =" $RESUB                         # Not used
endif

echo "=============Processing================="

if (${?NCI_CP} == 1) then
echo "NCI_CP =" $NCI_CP                       #
endif
if (${?NCI_ID} == 1) then
echo "NCI_ID =" $NCI_ID                       #
endif
if (${?CNV2NC} == 1) then
echo "CNV2NC =" $CNV2NC                       # convert UM ppfiles to NetCDF before running umplot
endif
if (${?PROC} == 1) then
#echo "PROC   =" $PROC                         # Not used
endif
if (${?Pmonth} == 1) then
echo "Pmonth =" $Pmonth                       # Ext. letter for Files containing Monthly Means
endif
if (${?Pdaily} == 1) then
echo "Pdaily =" $Pdaily                       # Ext. letter for Files containing Daily Means & Min/Max
endif
if (${?Ptemp1} == 1) then
echo "Ptemp1 =" $Ptemp1                       # Ext. letter for Files containing Daily Means & Min/Max
endif
if (${?Ptemps} == 1) then
echo "Ptemps =" $Ptemps                       # Ext. letter for Files containing Daily Means & Min/Max
endif
if (${?Ptimes} == 1) then
echo "Ptimes =" $Ptimes                       # Ext. letter for Files containing Daily Means & Min/Max
endif
if (${?PcasaC} == 1) then
echo "PcasaC =" $PcasaC                       # Ext. letter for Files containing Daily Means & Min/Max
endif

echo "=============Temp Varname==============="

if (${?tname} == 1) then
echo "tname  =" $tname                        # temp_`number` - variable name for tscrn
endif
if (${?txname} == 1) then
echo "txname =" $txname                       # temp_`number` - variable name for tscrn
endif
if (${?tiname} == 1) then
echo "tiname =" $tiname                       # temp_`number` - variable name for tscrn
endif

echo "=============Model Setup================"

if (${?MODEL} == 1) then
echo "MODEL  =" $MODEL                        # c/m - cable or moses
endif
if (${?RES} == 1) then
echo "RES    =" $RES                          # 48/96 - n48 or n96 resolution
endif
if (${?MASK} == 1) then
echo "MASK   =" $MASK                         # 1/2 - Which land-sea mask ?
endif
if (${?CAL} == 1) then
echo "CAL    =" $CAL                          # 360/365/366
endif
if (${?CASA} == 1) then
echo "CASA   =" $CASA                         # CASA-CNP
endif

echo "=============UMPLOT Setup==============="

if (${?UMPLT} == 1) then
echo "UMPLT  =" $UMPLT                        # y/n - do you want to run umplot ?
endif
if (${?SPLIT} == 1) then
echo "SPLIT  =" $SPLIT                        # y/n - ONLY if YR is a multiple of 5 and  >= 10
endif
if (${?BMRK} == 1) then
echo "BMRK   =" $BMRK                         # y/n
endif
if (${?BDIR} == 1) then
echo "BDIR   =" $BDIR                         # MOSES: ~ste69f/access/xaiyl_m48_20yrs .OR. ~ste69f/access/xagpb_n96_m_20yrs
endif

echo "=============Format Setup==============="

if (${?CHFMT} == 1) then
echo "CHFMT  =" $CHFMT                        # Not used yet
endif
if (${?FMT} == 1) then
echo "FMT    =" $FMT                          # pdf/eps/jpg - Not used yet
endif

#echo "=============Extra Envars==============="

#echo "USERID =" $USERID                       #
#echo "DIRW   =" $DIRW                         #
#echo "TILE   =" $TILE                         #
#echo "SOIL   =" $SOIL                         #
#echo "FullYrs=" $FullYrs                      #
#echo "TOP    =" $TOP                          #
#echo "MID    =" $MID                          #
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

