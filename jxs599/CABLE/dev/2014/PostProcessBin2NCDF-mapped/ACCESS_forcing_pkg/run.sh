#!/bin/csh

# this script calls a legacy PERL script with 3 arguments:
#  1 - the $basename of the variable in the file which should (legacy) also be 
#      the same as the name of the directory
#  2 - this is generically the same as the $basename, however for mapping files
#      the $data_type == 'mapping'. generically defines the output directory 
#  3 - the number of nodes

set nodes = 16
# by default assumes you already have processed the mapping data and don't wish to rerun this
set data_type = 'default'
# uncomment this if you DO wnat to process the mapping data 
set data_type = 'mapping'

#jhan - hack  - duplicated def in perl script
set catted = `pwd`'/data/catted'
set ncdf = `pwd`'/data/ncdf'
set dir_mapping = `pwd`'/data/mapping'

# process the mapping data 
if( $data_type == 'mapping' ) then
     /bin/cp $dir_mapping/longitude001.bin $catted/longitude.bin
     /bin/cp $dir_mapping/longitude001.dat $catted/longitude.dat
     #jhan: commneting as only need 1 file as nprocessor breakdown =1,16
     #"longitude"         \
   foreach basename (   \
      "latitude"          \
      "lat_index"         \
      "lon_index"         \
      "tile_frac"         \
      "tile_index"        \
      )
      echo ""
      echo "concatenating data from $nodes nodes for $basename ..."
      echo "data_type $data_type"
      echo ""
     

      ./fortran2NCDF_map.pl $basename $data_type $nodes 
      
      if( -f $catted/$basename.bin ) then
         echo "Done "
      endif

   end

endif

echo ""
echo ""
echo ""
foreach basename (   \
   tile_frac         \
   )
      echo ""
      echo "concatenating data from $nodes nodes for $basename"
      echo "in directory $basename"
      echo " and generating netcdf file ..."

      ./fortran2NCDF_map.pl $basename $basename $nodes 
      
      if( -f ../ncdf/$basename.nc ) then
         echo "Done "
         echo ""
      endif

end

