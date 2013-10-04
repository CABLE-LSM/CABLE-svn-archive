#!/usr/bin/perl
#defauls
$write_flag = '.false.'; 
$plot_flag = '.false.'; 
$time_window = 1;
$landpoint_window = 1;
$scratch_flag = '.false.';

$basename1 = $ARGV[0];		
$basename2 = $ARGV[1];		
$basename3 = $ARGV[2];		
$basename4 = $ARGV[3];		
if(-d "scratch" ) { print(" there is already a scratch directory. existing files will be overwritten."); }
else { system("mkdir scratch"); }
print("\n");
if(-d $basename1 ){ print(" there is already a directory with this basename. existing files will be overwritten.");}
else { system("mkdir $basename1"); }
print("\n");
if(-d $basename2 ){ print(" there is already a directory with this basename. existing files will be overwritten.");}
else { system("mkdir $basename2"); }
print("\n");

foreach $argnum (4 .. $#ARGV) {
   if($ARGV[$argnum] eq '-w' )   { $write_flag        = '.true.'        ; next; }; 
   if($ARGV[$argnum] eq '-p')    { $plot_flag         = '.true.'        ; next; };		
   if($ARGV[$argnum] eq '-s')    { $scratch_flag         = '.true.'        ; next; };		
   if($ARGV[$argnum] eq '-x' )   { $landpoint_window  = $ARGV[$argnum+1]; next; };
   if($ARGV[$argnum] eq '-t' )   { $time_window       = $ARGV[$argnum+1]; next; };
}

open(FILEPTR,">input.dat");
  print(FILEPTR"$basename1\n");
  print(FILEPTR"$basename2\n");
  print(FILEPTR"$basename3\n");
  print(FILEPTR"$basename4\n");
  print(FILEPTR"$write_flag\n");
  print(FILEPTR"$plot_flag\n");
  print(FILEPTR"$time_window\n");
  print(FILEPTR"$landpoint_window\n");
close(FILEPTR);

system("./diag_main");

if($write_flag eq '.true.' ) {
      system("mv $basename1.txt $basename1/");
      system("mv $basename2.txt $basename2/");
      system("mv $basename3.txt $basename3/");
      system("mv $basename4.txt $basename4/");
}
#open(FILEPTR,"<$basename.dat");
#  $trash = <FILEPTR>;
#  $Nvars = <FILEPTR>;
#  $trash = <FILEPTR>;
#  foreach $i(1 .. $Nvars) {
#      $var[$i] = <FILEPTR>;
#      chomp($var[$i])
#   }
#close(FILEPTR);
#
#$ext[1] ="m";
#chdir(scratch);
#if($plot_flag eq '.true.' ) {
#   foreach $i(1 .. $Nvars) {
#      $ext[2] ="$var[$i]";
#      $tseriesname = join("", @ext ); 
#      system("tseries.csh $var[$i] $tseriesname");
#      system("mv $tseriesname.mpeg ../$basename");
#   }
#}
#chdir("../");
#system("rm -f input.dat");
if($scratch_flag eq '.false.') { system("rm -fr scratch") }
if($scratch_flag eq '.true.') { 
   chdir("$basename");
   system("mv ../scratch .") 
}
#
#


































