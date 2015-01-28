
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

 fyr  = systemfunc("echo $FullYrs")
 yrs  = stringtoint(fyr)
 print("")
 print("Calculating Mod of FullYrs")
 ;print("")
 myrs = mod(yrs,5)
 ;system("echo "+myrs+" > mod_tmp")
 if ( myrs .eq. 0 ) then
  status_exit(38)
 end if

end

