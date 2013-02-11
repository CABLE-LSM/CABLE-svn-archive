#!/usr/bin/perl
#6 April 2009
#change the following ABSOLUTE PATH to your CRM directory

#$the_path="/Users/ludi/RSYNCED/Software/f90tohtml/examples/CABLE-test3/";
chomp($the_path = @ARGV[0]);
-e $the_path || die ("$the_path not set right");

$outdir="cable_ls"; #relative path to output directory
mkdir($outdir,0755) if (! -e $outdir);


@the_dirs=(	"core/src/science/",
		"core/src/common/"
			);

$the_files="*.F90";
foreach $dir (@the_dirs){
	@split_dir=split /\//,$dir;
print "\n";
print "\nsplit_dir\n",@split_dir;
print "\n";
	$title=pop @split_dir;
print "\ntitle\n",$title;
print "\n";
	if ($title eq "") { $title="main"};
	$file_name="$outdir/".$title.".ls";
print "\nfilename\n",$file_name;
print "\n";
	$ls_opt=$the_path.$dir.$the_files;
	print "\npreparing $file_name\n\t from $ls_opt\n";
	system "ls $ls_opt >$file_name " || die ("cannot write $file_name \n");
}

$dir=("offline/src/");
$nme=("offline");
$the_files="*.f90";
  @split_dir=split /\//,$dir;
print "\n";
print "\nsplit_dir\n",@split_dir;
print "\n";
  $title="offline";
print "\ntitle\n",$title;
print "\n";
  if ($title eq "") { $title="main"};
  $file_name="$outdir/$nme.ls";
print "\nfilename\n",$file_name;
print "\n";
  $ls_opt=$the_path.$dir.$the_files;
  print "\npreparing $file_name\n\t from $ls_opt\n";
  system "ls $ls_opt >$file_name " || die ("cannot write $file_name \n");



$dir=("UM/src/");
$nme=("UM");
$the_files="*.F90";
  @split_dir=split /\//,$dir;
print "\n";
print "\nsplit_dir\n",@split_dir;
print "\n";
  $title="UM";
print "\ntitle\n",$title;
print "\n";
  if ($title eq "") { $title="main"};
  $file_name="$outdir/$nme.ls";
print "\nfilename\n",$file_name;
print "\n";
  $ls_opt=$the_path.$dir.$the_files;
  print "\npreparing $file_name\n\t from $ls_opt\n";
  system "ls $ls_opt >$file_name " || die ("cannot write $file_name \n");

$dir=("Mk3L/src/");
$nme=("Mk3L");
$the_files="*.f90";
  @split_dir=split /\//,$dir;
print "\n";
print "\nsplit_dir\n",@split_dir;
print "\n";
  $title="Mk3L";
print "\ntitle\n",$title;
print "\n";
  if ($title eq "") { $title="main"};
  $file_name="$outdir/$nme.ls";
print "\nfilename\n",$file_name;
print "\n";
  $ls_opt=$the_path.$dir.$the_files;
  print "\npreparing $file_name\n\t from $ls_opt\n";
  system "ls $ls_opt >$file_name " || die ("cannot write $file_name \n");









#$file_name="$outdir/include.ls";
#unlink($file_name) || warn ("cannot unlink $file_name");
#$the_files="*.h";
#foreach $dir (@the_dirs){
#	$ls_opt=$the_path.$dir.$the_files;
#	print "\nappending to $file_name\n\t from $ls_opt\n";
#	system "ls $ls_opt >>$file_name " || die ("darn, no write $file_name \n");
#}
print "DONE\nNow try typing f90tohtml crm.f2h\n";

