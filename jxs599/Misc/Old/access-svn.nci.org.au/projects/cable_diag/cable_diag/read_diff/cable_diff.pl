#!/usr/bin/perl      

###############################################################################
###############################################################################
###############################################################################
###                                                                         ###
### perl script called from command line to  interepret binary data written ###
### by GCM using 'call cable_diag(.., .., .., ..,)'. command line usage:    ###
###                                                                         ###
### cable_diag basename [-s -w -p] [-x nx] [-t nt]                          ###
###   optoins:                                                              ###
###            -     with no options specified ciable_diag will make the    ###
###                  directories and exit.                                  ###
###            -w    will write the binary data in basename.bin to a        ### 
###                  textfile basename.txt                                  ###
###            -p    enables plotting/ffmpeg of timeseries data             ###
###                  plotting uses pgplot libraries concatenation of .png   ###  
###                  plot per timestep files uses ffmpeg libraries          ###  
###                  ensure these are available                             ### 
###            -s    keep scratch directory [default=NO]. scratch directory ###             
###                  is used to store plots per timestep, which are then    ###
###                  concatenated to generate mpeg ( varname.mpeg)          ###
###            -x    enables 'smoothing' over x dimension (usually          ###
###                  landpoints)                                            ###
###            nx    must be specified if -x enabled. nx = variance of      ###  
###                  Gaussian filter convoluted with data. nx is in units of###
###                  number of cells                                        ###
###            -t    enables 'smoothing' over t dimension (usually          ###
###                  time)                                                  ###
###            nt    must be specified if -t enabled. nt = number of time-  ###  
###                  steps averaged over before .png plots are done         ###
###                                                                         ###
###############################################################################
###############################################################################
###############################################################################


###******this is for ruff'n'ready diffing of basename/basename2 only

#use: cable_diff.pl basename basename2

### defaults. changed through command line
#$write_flag = '.false.'; 
#$plot_flag = '.false.'; 
#$time_window = 1;
#$landpoint_window = 1;
#$scratch_flag = '.false.';

### the first arg on the command line is the basename of the desired data to process
### i.e. basename.dat, basename.bin
$basename = $ARGV[0];	
$basename2 = $ARGV[1];	

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

#### write a little file for consumption by f90 code, interpreting command line args
open(FILEPTR,">input.dat");
  print(FILEPTR"$basename\n");
  print(FILEPTR"$basename2\n");
#  print(FILEPTR"$write_flag\n");
#  print(FILEPTR"$plot_flag\n");
#  print(FILEPTR"$time_window\n");
#  print(FILEPTR"$landpoint_window\n");
close(FILEPTR);

### execute f90 program to do all the work desired
system("./diff_main");

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




































