# This script is called by nix
# and sets up the evironment
unset PATH
for p in $baseInputs $buildInputs; do
  if [ -d $p/bin ]; then
    export PATH="$p/bin${PATH:+:}$PATH"
  fi
#  if [ -d $p/include ]; then
#    export NIX_CFLAGS_COMPILE="-I $p/include${NIX_CFLAGS_COMPILE:+ }$NIX_CFLAGS_COMPILE"
#  fi
#  if [ -d $p/lib ]; then
#    export NIX_LDFLAGS="-rpath $p/lib -L $p/lib${NIX_LDFLAGS:+ }$NIX_LDFLAGS"
#  fi
done
# we export some environmentvariables needed by our version "build_nix.ksh"# of build.ksh  
# we could also change the makefile to directly incorporate them
export NCDIR="$ncd/lib${NCDIR:+ }$NCDIR"
export NCMOD="$ncd/include${NCMOD:+ }$NCMOD"

# define some functions 
function unpackPhase(){
	echo "nothing to unpack"
}

function buildPhase(){
	echo "#################################"
	echo $ctn
	cd "${ctn}/offline"
	rm -rf ".tmp"
	mksh build_nix.ksh #This will create an executable with name "cable"
	#mksh build_mpi_nix.ksh #This will create an executable with name "cable"
}
function installPhase(){
	mv cable $out
	#mv cable-mpi $out
}
function fixupPhase(){
	#clean_build
	find $out -type f -exec patchelf --shrink-rpath '{}' \;-exec strip '{}' \; 2>/dev/null
}
function genericBuild(){
	unpackPhase
	buildPhase
	installPhase
	#fixupPhase
}
