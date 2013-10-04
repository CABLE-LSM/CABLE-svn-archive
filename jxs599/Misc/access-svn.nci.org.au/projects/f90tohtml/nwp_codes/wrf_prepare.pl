#!/usr/bin/perl
#CHANGE THE FOLLOWING PATH TO YOUR WRF DIRECTORY
$the_path="/home/bfiedler/WRF/WRFV3/";

-e $the_path || die ("$the_path not set right");

$outdir="wrf_ls"; #relative path to output directory

mkdir($outdir,0755) if (! -e $outdir);

@the_F_dirs=(	"dyn_em/", "dyn_exp/", "dyn_nmm/", "frame/", "main/", "phys/", "share/","external/io_grib_share/", "external/io_netcdf/","external/io_grib1/","external/io_grib2/");

@the_inc_dirs=("arch/","dyn_nmm/","inc/","share/","tools/","external/io_grib_share/");

@ignores=("module_initialize","module_sf_lsm_nmm","adve_optim\.h","module_cu_kfeta\.F","module_cu_gd.F");
@keeps=("module_initialize_b_wave","addStringHereToOveruleIgnore");

sub matchany{
	my ($string,$matches)=@_;
	foreach $match (@$matches){
		return 1 if $string=~m/$match/;
	}
	return 0;
}
	
	
$the_files="*.F";
foreach $dir (@the_F_dirs){
	@split_dir=split /\//,$dir;
	$title=pop @split_dir;
	$ls_file_name="$outdir/".$title.".ls";
	$ls_opt=$the_path.$dir.$the_files;
	print "\npreparing $ls_file_name\n\t from $ls_opt\n";
	open LSFILE, ">$ls_file_name";
	@fn=glob $ls_opt."*";
	foreach $f90file (@fn){
		print LSFILE $f90file,"\n" unless &matchany($f90file,\@ignores) and !&matchany($f90file,\@keeps);
     }
}

$the_files="*.h";
$ls_file_name="$outdir/"."inc.ls";
open LSFILE, ">$ls_file_name";
foreach $dir (@the_inc_dirs){
	$ls_opt=$the_path.$dir.$the_files;
	@fn=glob $ls_opt."*";
	print "\npreparing $ls_file_name\n\t from $ls_opt\n";
	foreach $hfile (@fn){
		print LSFILE $hfile,"\n" unless &matchany($hfile,\@ignores) and !&matchany($hfile,\@keeps);
	}
}

print "DONE\nNow try typing f90tohtml wrf.f2h\n";
print "(but you may later want to comment out troublesome file names in the .ls files)\n"

