load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

wks = gsn_open_wks("ps","cmap")
;wks = gsn_open_wks("x11","cmap")

;cmapr  = (/"white","black","cornsilk3","skyblue","dodgerblue","darkolivegreen2","chartreuse3","lightgoldenrod1","gold","orange","darkorange2","red","firebrick"/)

cmap = RGBtoCmap("/g/data1/p66/lxs599/umplot/colormap/BlueOrangeRed.rgb")
;cmap = RGBtoCmap("/home/cmar/ste69f/umplot/colormap/BlueOrangeRed.rgb")
;print(cmap)

;gsn_define_colormap(wks,cmapr)
gsn_define_colormap(wks,cmap)

gsn_draw_colormap(wks)

end

