#!/usr/bin/perl      

### the first arg on the command line is the basename of the desired data to process
### i.e. basename.dat, basename.bin
#jhan commented this to speed up testing for release
#$basename = $ARGV[0];	
#$basename2 = $ARGV[1];	

$basename = "std";	
$basename2 = "FLUXES00";	

if(-e "std.dat") { 
  print("\nproceed ...\n"); 
}
else{           
  print("\nyou need to link to std files to compare\n"); 
   exit;
}      

#### write a little file for consumption by f90 code, interpreting command line args
open(FILEPTR,">input.dat");
  print(FILEPTR"$basename\n");
  print(FILEPTR"$basename2\n");
close(FILEPTR);

### execute f90 program to do all the work desired
system("diff_main");





































