#!/usr/bin/perl      


###################################################################################
##### The motivation for creating this "method" is to be able to access variables 
##### from within CABLE (UM only at present), along with their spatial location, 
##### vegetation type and corresponding fraction. Printing to UM output is messy, 
##### expensive, and cumbersome. This "method" gathers all the aforementioned data
##### and creates a netcdf file for viewing via ncl, ferret, etc. It works across 
##### an arbitrary number of processors. It is adaptable in terms of how often it 
##### records data, however keep in mind time is the most expensive dimension here.
####
##### this script is the post-processing counterpart of the restricted cable_diag()
##### fortran subroutine. in this instance the work flow proceeds as:
##### 1. from within( CABLE in the UM)  
##### this script lives in ~jxs599/Public/bin. the user will need to copy this 
##### script for each new SET of jobs. using an example, you may run an experiment
##### USE: ./cable_ncdf.pl basename                                           #####
#####   where basename is the basename of the file/variable you want to convert ###
##### The USER section of this script should also be configured   
###################################################################################


use Cwd; #needed to get pwd in perl
$dir_exec = &Cwd::cwd();

###################################################################################
####   USER SECTION to be configured     ##########################################
###################################################################################

#you can copy the executables to anywhere, but i suggest you use these 
$path_exec = "/home/599/jxs599/Public/bin/";

#in general absolute path to the cable_diag .dat/.bin files you have created
#which is most often in /short/p66/$userID/$modelID

#$dir_mapping00 = "/home/599/jxs599/UM_ROUTDIR/laaaf/1/";
$dir_mapping00 = "$dir_exec/";

#the number of compute nodes the model was run on which produced your .dat/.bin files
$nodes = 1;


#absolute path to the directory you wish to place the mapped netcdf files which put 
#cable vars onto a global map. The files in  "/home/599/jxs599/Public/mapping/" correspond to 
# UM N96 resolution using the version2 landmask. If this is the same configuration which 
#you use then you do not need to compute these again and may leave the below switch
#"$first_mapping" set to "false". otherwise set a new path for yourself an set this 
#switch to "true"to generate these files. Note however that you only need to compute 
#these once for each resolution, NOT for each variable you look at.  
#this may also be in $dir_mapping00, however it is likely to be very messy if you do 
#NB. THIS DIRECTORY INPUT TO ncdf_main.f90 

$dir_mapping = "$dir_exec/mapping/";

#IF this is the first mapping you have done for this resolution then you will need 
#to switch this to true. NOTE: in all likelihood this should be false most of the time,
#and $dir_mapping unchanged 
$first_mapping = "true";


###################################################################################
####   END USER SECTION to be configured     ######################################
###################################################################################




### the first arg on the command line is the basename of the desired data to process
### i.e. basename.dat, basename.bin
$basename = $ARGV[0];	



###IF true - then concatenates files and puts in $dir_mapping
###IF false then ncdfmain() simply uses prepared files in $dir_mapping


if ( $first_mapping eq "true") { 

   chdir($dir_mapping00);
   system("$path_exec/cable_diag_cat_nodes.pl latitude -n $nodes" );
   system("mv latitude.* $dir_mapping" );

   system("$path_exec/cable_diag_cat_nodes.pl lat_index -n $nodes" );
   system("mv lat_index.* $dir_mapping" );
   
   system("$path_exec/cable_diag_cat_nodes.pl lon_index -n $nodes" );
   system("mv lon_index.* $dir_mapping" );
   
   system("$path_exec/cable_diag_cat_nodes.pl tile_index -n $nodes" );
   system("mv tile_index.* $dir_mapping" );
   
   system("$path_exec/cable_diag_cat_nodes.pl tile_frac -n $nodes" );
   system("mv tile_frac.* $dir_mapping" );

system("/bin/cp longitude00.dat longitude.dat" );
system("/bin/cp longitude00.bin longitude.bin" );
system("mv longitude.* $dir_mapping" );
   
}

chdir($dir_mapping00);



system("$path_exec/cable_diag_cat_nodes.pl $basename -n $nodes" );
system("mv $basename.* $dir_exec" );

#might not need this if we put in ~/bin
chdir($dir_exec);


#### write a little file for consumption by f90 code, interpreting command line args
open(FILEPTR,">input.dat");
  print(FILEPTR"$basename\n");
  print(FILEPTR"\"$dir_mapping\"\n");
#  print(FILEPTR"$plot_flag\n");
#  print(FILEPTR"$time_window\n");
#  print(FILEPTR"$landpoint_window\n");
close(FILEPTR);

#might not need this if we put in ~/bin
chdir($dir_exec);
   

### execute f90 program to do all the work desired
system("$path_exec/ncdf_main");


