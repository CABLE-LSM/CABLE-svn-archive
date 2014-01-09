#!/bin/csh

#foreach i ( 6T00 6B00 )
#   ln -s ../../../../../xaeld/$i.dat d$i.dat
#   ln -s ../../../../../xaeld/$i.bin d$i.bin
#end

foreach i ( 6T 6B )
   ln -s ../../../cat_N_nodes/read/xaeqf/$i.dat f$i.dat
   ln -s ../../../cat_N_nodes/read/xaeqf/$i.bin f$i.bin
end


#foreach i ( 4T00 4B00 )
#   ln -s ../../../../xaeld/$i.dat d$i.dat
#   ln -s ../../../../xaeld/$i.bin d$i.bin
#end
