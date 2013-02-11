#!/usr/bin/perl      

###############################################################################
###############################################################################

###############################################################################
###############################################################################


### defaults. changed through command line
$n_nodes = 1;

### the first arg on the command line is the basename of the desired data to process
### i.e. basename.dat, basename.bin
$basename = $ARGV[0];	
### check if stuff already exists - overwrite if it does
### get remaining command line args if there are any
### i.e. cable_diag_cat_nodes basename -n N_nodes
foreach $argnum (1 .. $#ARGV) {
   if($ARGV[$argnum] eq '-n' )   { $n_nodes           = $ARGV[$argnum+1]; next; };
}

### write a little file for consumption by f90 code, interpreting command line args
open(FILEPTR,">input.dat");
  print(FILEPTR"$basename\n");
  print(FILEPTR"$n_nodes\n");
close(FILEPTR);

### execute f90 program to do all the work desired
system("/home/599/jxs599/Public/bin/diag_cat_nodes");



































