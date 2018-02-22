begin

npts = 50
fin=addfile("landmask.nc","r")
landsea=fin->landsea

landsea_1d = ndtooned(landsea)

gd_pts = ind(landsea_1d .eq. 0)

landsea_1d_all_land = landsea_1d(gd_pts)


landsea_1d_all_land(:)=1

random_setallseed(36484749, 9494848)               ; Set seeds (suggested, NOT required)
low   =  0.0
high  = tofloat(dimsizes(landsea_1d_all_land)-1)
 
rnd_pts   = tointeger(random_uniform(low, high, (/npts/))) 

landsea_1d_all_land(rnd_pts(:)) = 0

landsea_1d(gd_pts)=landsea_1d_all_land

landsea=(/reshape(landsea_1d,(/360,720/))/)

fout=addfile("gswp3_landmask_"+npts+"pts.nc","c")

fout->lat=fin->lat
fout->lon=fin->lon
fout->landsea=landsea
delete(fout)
delete(fin)


end

