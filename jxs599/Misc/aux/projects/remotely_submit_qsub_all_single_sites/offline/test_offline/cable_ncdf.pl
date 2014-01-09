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
#$path_exec = "/home/599/jxs599/Public/bin/";
$path_exec = "./";

###################################################################################
####   END USER SECTION to be configured     ######################################
###################################################################################



### the first arg on the command line is the basename of the desired data to process
### i.e. basename.dat, basename.bin
$basename = $ARGV[0];	

#### write a little file for consumption by f90 code, interpreting command line args
open(FILEPTR,">input.dat");
  print(FILEPTR"$basename\n");
close(FILEPTR);

### execute f90 program to do all the work desired
system("$path_exec/ncdf_main");


