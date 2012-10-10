; Run Plotting Script ======================================================
; Lauren Stevens 21 April 2011
; ==========================================================================

load "/home/cmar/ste69f/umplot/variables.ncl"

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

  ; LDIR or ODIR: Comparison Run
  ;Cab9tn48   = "/home/cmar/kow014/access/xaank_run244"         ; N48
  ;Cab17tn96  = "/home/cmar/kow014/access/xaanv_run32"          ; N96
  ;;Mos9tn48   = "/home/cmar/ste69f/access/xafqb_n48_m_20yrs"    ; N48
  ;Mos9tn48   = "/home/cmar/ste69f/access/xaiya_n48_m_20yrs"    ; N48
  ;Mos9tn96   = "/home/cmar/ste69f/access/xagpb_n96_m_20yrs"    ; N96
  ;Mos9tn96v2 = "/home/cmar/ste69f/access/xaiyd_n96_m_20yrs"    ; N96 v2 mask

  dir  = systemfunc("echo $DIR")
  rid  = systemfunc("echo $RUNID")
  year = systemfunc("echo $YR")  
  tile = systemfunc("echo $TILE")
  soil = systemfunc("echo $SOIL")
  resn = systemfunc("echo $RES")

  jpg = systemfunc("echo $JPEG")
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
  system("/home/cmar/ste69f/umplot/dir.sh")
  end if

  lmod  = systemfunc("echo $MODEL")
  if (lmod.eq."c") then
   lcable = True
  else
   lcable = False
  end if

  extra = systemfunc("echo $TAY")
  if (extra.eq."y") then
   tay  = True
   odir = systemfunc("echo $LDIR")
  else
   tay  = False
   odir = systemfunc("echo $HOME")
  end if

; Call variables.ncl
  run_ncl(dir,rid,year,tile,soil,"n"+resn,split,tay,odir,lcable,jpeg)

;Examples:
  ;dir    = "./" ;"/home/cmar/ste69f/access/xaiip"
  ;dir    = "/home/cmar/kow014/access/xaank"
  ;rid    = "xaiip"
  ;year   = "5"
  ;tile   = "17"
  ;soil   = "6"
  ;resn   = "96"
  ;split  = False ; if (split) then yr = year/5 ; do i =0, yr-1,1 k=i+1 system("mkdir block"+k+"_5yrs") end do end if; system("./cdo_copy2.sh") etc
  ;tay    = False
  ;odir   = "./"
  ;lcable = True
  ;jpeg   = False

  ;run_ncl("./","xagpb","20","9","4","n96",False,False,"./",False,False) 
  ;                                                 ; 20 year n96 Moses Run ( 9 tiles, 4 soil levs )
  ;                                                 ; Run for full 20 years with no special plots
  ;
  ;run_ncl("./","xaiia","5","17","6","n96",False,False,"./",False,False)
  ;run_ncl("./","xaanv","10","17","6","n96",False,False,"./",False,False)
  ;run_ncl("./","xaank","10","9","4","n48",False,False,"./",False,False)

  ; Currently Split functionality not working and
  ; also Special plots not working for NCL control script

end

