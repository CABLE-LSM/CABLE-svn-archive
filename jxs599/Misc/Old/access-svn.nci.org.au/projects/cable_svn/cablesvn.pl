#!/usr/bin/perl

# Please see the README file at 
#
#  https://access-svn.nci.org.au/trac/cable/browser/README 
#
# for full description of this script
#
# Basically the CABLE repository contains more material than any one user will probably want. 
# cablesvn is a perl script intended to guide you through the CABLE subversion repository and CHECK OUT 
# the WORKING COPY you require. At present the script is restricted to copying the head of the trunk.  
# It is entirely possible for you to access the repository directly if you would prefer at: 
#
# https://access-svn.nci.org.au/trac/cable/browser/README 
#
# or even configure this script directly to your liking. 
# The reality is that you will probably use this script only a few times at most, and all that i
# it really does is bundle together some svn commands that wold otherwise be cumbersome.

# Please, at least configure the top section of this script. 



use Cwd; #module to get current working directory


####################################
### USERS CONFIGURE THIS SECTION ###
####################################

#path to svn binary
$svn = "/opt/subversion/bin/svn"; 

#CABLE user name
$userid = "ludi";

#desired name of working copy.  
#default = "none", which will give working copy same name name as user branch 
# (e.g. CABLE-2.0_r311, where 311 is revison number of copied branch)
$co_name = "none";

#currently the repository contains directories to run CABLE:
#in an offline capacity, and in global models UM and Mk3L. 
#"$usertype" indicates which parts of the repository to checkout. 
#please uncomment the approriate line.
$usertype = "all";
#$usertype ="mk3l";
#$usertype ="offline";
#$usertype = "UM";
#$usertype = "mk3loffline";
#$usertype = "UMoffline";
 
###################################
### END USERS CONFIGURE SECTION ###
###################################


### Static Declarations ###

# repository paths
$svnroot = "https://access-svn.nci.org.au/svn/cable"; 

$svntrunk = "$svnroot/trunk/"; 
$svnbranches= "$svnroot/branches"; 
$svndev= "$svnbranches/dev/"; 
$svnShare = "$svndev/Share/"; 

$svnuser = "$svndev/users/$userid/"; 



### particular revision/version branches 
# added manually for convenience as project speed allows it

$svnSharebeta = "$svnShare/bCABLE-1.9/"; 
$cstamp = 'CABLE-1.9';





### main program

&toplevel($svntarget,$newbranch,$stamp);

if ($co_name eq "none") {
   $label = $stamp; 
} else {
   $label = $co_name; 
}

&svncall($svntarget,$newbranch,$label,$usertype);

### END of main program





sub toplevel($){
   system("clear");

   ### pre define some long stretches of text ###
   $overview = 
"The default behaviour of this script is to copy the trunk at the head revision to the user's branch, tagged with the revision number it was copied from. A local working copy is then checked out for you, depending on your configuration at the top of this script. \n\n Press enter to continue with this default behaviour, or abort and email us at \n\n
      
 for something more complicated that you might need help with.\n\n Please also see the README file at\n\n https://access-svn.nci.org.au/trac/cable/browser/README\n\n";



   print ( "\n\n$overview" );  
   $cont= <STDIN>;

   system( "clear" );
   
   $svntarget = $svntrunk;
   #jhan: for testing only   
   $svntarget = "$svnuser/Aug11/CABLE-test3/";
   $svntarget= $svnSharebeta; 
   
   $head_tr = `$svn info $svntarget -r HEAD | grep Rev | grep Last |cut -d : -f 2`;
   $head_tr = trim($head_tr);   

   $ext[0] = $cstamp; 
   $ext[1] = "_"; 
   $ext[2] = "r"; 
   $ext[3] = $head_tr;
    
   #join ext[] in one string
   $stamp= join("",@ext); 

   $iext[0] = $svnuser; 
   $iext[1] = $stamp; 
    
   $newbranch= join("",@iext); 
   #jhan: for testing only   
   $newbranch= $svnSharebeta; 
}





sub svncall($){ 
   system("clear");

   ### pre define some long stretches of text ###
   $confirmmess = 
   "You are about to copy \n\n$svntarget\n\n to \n\n$newbranch\n\n and check this out to $label. ($label will be created for you)";

   print( "\n\n$confirmmess\n\n" );
   print( "\n\n Press enter to continue. Otherwise abort with Control-C:\n");
   $finalok=<STDIN>;


   #svn log message
   $svnlog = "-m \"copy $svntarget to user branch\"";


   ##copy svn target to user branch
#   system( "svn copy $svntarget $newbranch $svnlog" ); 
      
   if($usertype eq "all") {
      system( "$svn checkout $newbranch $label" );
   } 
   elsif ($usertype eq "mk3l") {
      system( "$svn checkout $newbranch $label --depth empty" ); 
      chdir($label);      
      system( "$svn update --set-depth infinity core" ); 
      system( "$svn update --set-depth infinity Mk3L" ); 
   } 
   elsif ($usertype eq "offline") {
      system( "$svn checkout $newbranch $label --depth empty" ); 
      chdir($label);      
      system( "$svn update --set-depth infinity core" ); 
      system( "$svn update --set-depth infinity offline " ); 
   } 
   elsif ($usertype eq "UM") {
      system( "$svn checkout $newbranch $label --depth empty" ); 
      chdir($label);      
      system( "$svn update --set-depth infinity core" ); 
      system( "$svn update --set-depth infinity UM" ); 
   } 
   elsif ($usertype eq "mk3loffline") {
      system( "$svn checkout $newbranch $label --depth empty" ); 
      chdir($label);      
      system( "$svn update --set-depth infinity core" ); 
      system( "$svn update --set-depth infinity Mk3L " ); 
      system( "$svn update --set-depth infinity offline " ); 
   } 
   elsif ($usertype eq "UMoffline") {
      system( "$svn checkout $newbranch $label --depth empty" ); 
      chdir($label);      
      system( "$svn update --set-depth infinity core" ); 
      system( "$svn update --set-depth infinity UM" ); 
      system( "$svn update --set-depth infinity offline " ); 
   } 
   else { 
      ( "\n\nPlease configure usertype variable at top of this script.\n\n"); 
   }            
}





# Perl trim function to remove whitespace from the start and end of the string
sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}

