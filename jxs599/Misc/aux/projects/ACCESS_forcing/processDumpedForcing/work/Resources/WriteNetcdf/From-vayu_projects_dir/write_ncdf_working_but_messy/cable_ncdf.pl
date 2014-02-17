#!/usr/bin/perl      

###******this is for ruff'n'ready for writing basename to netcdf

#use: cable_ncdf.pl basename 

### defaults. changed through command line
#$write_flag = '.false.'; 
#$plot_flag = '.false.'; 
#$time_window = 1;
#$landpoint_window = 1;
#$scratch_flag = '.false.';

### the first arg on the command line is the basename of the desired data to process
### i.e. basename.dat, basename.bin
$basename = $ARGV[0];	
$newfile = 'test.nc';
#### check if stuff already exists - overwrite if it does
#if(-d "scratch" ) { print(" there is already a scratch directory. existing files will be overwritten."); }
#else { system("mkdir scratch"); }
#print("\n");
#if(-d $basename ){ print(" there is already a directory with this basename. existing files will be overwritten.");}
#else { system("mkdir $basename"); }
#print("\n");
#
#### get remaining command line args if there are any
#foreach $argnum (1 .. $#ARGV) {
#   if($ARGV[$argnum] eq '-w' )   { $write_flag        = '.true.'        ; next; }; 
#   if($ARGV[$argnum] eq '-p')    { $plot_flag         = '.true.'        ; next; };		
#   if($ARGV[$argnum] eq '-s')    { $scratch_flag         = '.true.'        ; next; };		
#   if($ARGV[$argnum] eq '-x' )   { $landpoint_window  = $ARGV[$argnum+1]; next; };
#   if($ARGV[$argnum] eq '-t' )   { $time_window       = $ARGV[$argnum+1]; next; };
#}

### write a little file for consumption by f90 code, interpreting command line args
open(FILEPTR,">input.dat");
  print(FILEPTR"$basename\n");
  print(FILEPTR"$newfile\n");
#  print(FILEPTR"$write_flag\n");
#  print(FILEPTR"$plot_flag\n");
#  print(FILEPTR"$time_window\n");
#  print(FILEPTR"$landpoint_window\n");
close(FILEPTR);

### execute f90 program to do all the work desired
system("./ncdf_main");

#### tidy up
#if($write_flag eq '.true.' ) {
#      system("mv $basename.txt $basename/");
#}
#
###get var. names from basename.dat to make an mpeg if we havve to
#open(FILEPTR,"<$basename.dat");
#  $trash = <FILEPTR>;
#  $Nvars = <FILEPTR>;
#  $trash = <FILEPTR>;
#  foreach $i(1 .. $Nvars) {
#      $var[$i] = <FILEPTR>;
#      chomp($var[$i]);
#   }
#close(FILEPTR);
#
#### define the name for each mpeg and call tseries.csh to do the concatenation of .png in scratch
#$ext[0] ="m";
#chdir(scratch);
#if($plot_flag eq '.true.' ) {
#   foreach $i(1 .. $Nvars) {
#      $ext[1] ="$var[$i]";
#      $tseriesname = join("", @ext ); 
#      chomp($tseriesname);
#      $extb[0] = $tseriesname;
#      $extb[1] = ".mpeg";
#      chomp($extb[1]);
#      $mname = join("",@extb );
#      system("tseries.csh $var[$i] $tseriesname");
#      system("mv $mname ../$basename");
#   }
#}
#chdir("../");
#
#### tidy up
#system("rm -f input.dat");
#if($scratch_flag eq '.false.') { system("rm -fr scratch") }
#if($scratch_flag eq '.true.') { 
#   chdir("$basename");
#   system("mv ../scratch .") 
#}

