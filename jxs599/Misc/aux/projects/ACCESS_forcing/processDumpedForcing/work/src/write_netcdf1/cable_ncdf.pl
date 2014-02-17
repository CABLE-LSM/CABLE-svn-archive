#!/usr/bin/perl      

###******this is for ruff'n'ready for writing basename to netcdf

#use: cable_ncdf.pl basename 

### the first arg on the command line is the basename of the desired data to process
### i.e. basename.dat, basename.bin
$basename = $ARGV[0];	
$newfile = 'metfsd1.nc';

### write a little file for consumption by f90 code, interpreting command line args
open(FILEPTR,">input.dat");
  print(FILEPTR"$basename\n");
  print(FILEPTR"$newfile\n");
close(FILEPTR);

### execute f90 program to do all the work desired
system("./ncdf_main");

