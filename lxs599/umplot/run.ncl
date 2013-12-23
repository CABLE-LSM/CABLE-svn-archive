; Run Plotting Script ======================================================
; Lauren Stevens 21 April 2011
; ==========================================================================

usr  = systemfunc("echo $USERID")
loadscript("~"+usr+"/umplot/variables.ncl")

begin

; INPUTS
  ;(fileDir:string,\   ; Working Directory
  ; Runid:string,\     ; Run Id
  ; numYr:string,\     ; Number of Years
  ; Tiles:string,\     ; Number of Tiles
  ; sLevs:string,\     ; Number of Soil Levels
  ; Res:string,\       ; Resolution
  ; Split:logical,\    ; Split Function i.e. block#
  ; Taylor:logical,\   ; Extra Plots (Special) 
  ; compDir:string,\   ; Comparison Directory
  ; lcable:logical,\   ; Model. True for CABLE, False for MOSES
  ; jpeg:logical)      ; Converts all files to jpeg for PowerPoint

  ; BDIR or ODIR: Comparison Run
  ;Cab9tn48   = "/home/cmar/kow014/access/xaank_run244"         ; N48
  ;Cab17tn96  = "/home/cmar/kow014/access/xaanv_run32"          ; N96
  ;;Mos9tn48   = "/home/cmar/ste69f/access/xafqb_n48_m_20yrs"    ; N48
  ;Mos9tn48   = "/home/cmar/ste69f/access/xaiya_n48_m_20yrs"    ; N48
  ;Mos9tn96   = "/home/cmar/ste69f/access/xagpb_n96_m_20yrs"    ; N96
  ;Mos9tn96v2 = "/home/cmar/ste69f/access/xaiyd_n96_m_20yrs"    ; N96 v2 mask

  ;usr  = systemfunc("echo $USERID")
  dir  = systemfunc("echo $DIR")
  rid  = systemfunc("echo $RUNID")
  year = systemfunc("echo $YR")  
  tile = systemfunc("echo $TILE")
  soil = systemfunc("echo $SOIL")
  resn = systemfunc("echo $RES")

  jpg = systemfunc("echo $CHFMT")
  if (jpg.eq."y") then
  jpeg = True
  else
  jpeg = False
  end if

  splt = systemfunc("echo $SPLIT") 
  if (splt.eq."y") then
   split = True 
  else
   split = False
  end if

  if(.not.split)then
  system("~"+usr+"/umplot/dir.sh")
  end if

  lmod  = systemfunc("echo $MODEL")
  if (lmod.eq."c") then
   lcable = True
  else
   lcable = False
  end if

  extra = systemfunc("echo $BMRK")
  if (extra.eq."y") then
   tay  = True
   odir = systemfunc("echo $BDIR")
  else
   tay  = False
   odir = systemfunc("echo $HOME")
  end if

; Call to variables.ncl
  run_ncl(dir,rid,year,tile,soil,"n"+resn,split,tay,odir,lcable,jpeg)

;Examples:
  ;run_ncl("./","xagpb","20","9" ,"4","n96",False,False,"./",False,False) 
  ;run_ncl("./","xaanv","10","17","6","n96",False,False,"./",False,False)
  ;run_ncl("./","xaank","10","9" ,"4","n48",False,False,"./",False,False)

  ; Currently Split functionality not working and
  ; also Special plots not working for NCL control script

end

