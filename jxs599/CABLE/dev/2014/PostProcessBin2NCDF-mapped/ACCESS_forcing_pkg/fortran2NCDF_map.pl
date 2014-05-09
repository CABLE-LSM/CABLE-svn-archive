#!/usr/bin/perl

#./fortran2NCDF_map./pl $basename $data_type $nodes 
#./fortran2NCDF_map./pl latitude mapping 16
#./fortran2NCDF_map./pl metfsd1 metfsd1 16

use Cwd; #needed to get pwd in perl

### the first arg on the command line is the basename of the desired data to process
### i.e. basename.dat, basename.bin
$basename = $ARGV[0];	
$data_type = $ARGV[1];
$nodes =  $ARGV[2];

$dir_pwd = &Cwd::cwd();
# root directories
# =================
# path to binaries cat_Nnodes & ncdf_main
$dir_exec = "$dir_pwd/bin";
#path to root of data directories
$dir_data = "$dir_pwd/data";

# specific data directories 
# =========================

# INPUT
# from ARGV[1] - always cd-s to this dir
# and operates on files therein according to ARGV[0], ARGV[2]
$dir_type= "$dir_data/$data_type/";
#print(" dir_type: ", $dir_type);

# OUTPUT
# where to put the generated files POST concatenation
$dir_catted = "$dir_data/catted/";
# where to put the generated NETCDF files 
$dir_ncdf = "$dir_data/ncdf/";
#rm# $dir_basename = "$dir_data/$basename/";

# changes to this dir to operate on local files therein
# set from calling argument ARGV[0]
chdir("$dir_type");

# first concatenate files/vars as passed (ARGV[0]) in CLI CALL to this script
# =================================================================

# call fortran executable to concatenate per node (ARGV[2]) data files
# and then move output

# in fortran exec - cat_Nnodes(filename, Nnodes)
system( "$dir_exec/cat_Nnodes $basename $nodes" ); 
system("/bin/cp $basename.* $dir_catted" );


# second write netcdf file/var for passed (ARGV[1]) in CLI CALL to this script
if( $data_type eq $basename) {

   #### write a little file for consumption by f90 code, interpreting command line args
   open(FILEPTR,">input.dat");
     print(FILEPTR"$basename\n");
     print(FILEPTR"\"$basename.nc\"\n");
   close(FILEPTR);
  
   # call fortran executable to generate mapped netcdf file
   # ncdf_main (dir_catted)
   system("$dir_exec/ncdf_main $dir_catted");
   system("mv $basename.nc $dir_ncdf" );

}

