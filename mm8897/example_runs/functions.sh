function run_gfortran_cable(){
	version=$1
	rm cable.nml
	ln -s "my_cable-${version}.nml" "cable.nml"
	../CABLE-${version}_nix/offline/nix_cable
}
function run_ifort_cable(){
	version=$1
	rm cable.nml
	ln -s "my_cable-${version}.nml" "cable.nml"
	source /home/mm/intel/bin/compilervars.sh -arch intel64
	NCDIR=/home/mm/local/lib
	NCMOD=/home/mm/local/include
	LD_LIBRARY_PATH=${NCDIR}:${LD_LIBRARY_PATH} 
	../CABLE-${version}_nix/offline/ifort_cable
}
