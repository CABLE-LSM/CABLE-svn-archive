#!/bin/sh

#Mark Decker 2017
#This script provides a very simple tool to subset the a land mask

show_help()
{
cat << EOF
usage: $0 options

This script sets all land points outside a bounding box to sea for regional simulations
Example usage:
./create_regional_mask.sh --input=oldmask.nc --output=newmask.nc  --latmin=-25.0 --latmax=2.5 --lonmin=285.0 --lonmax=335.0
OPTIONS:
  -h or --help  Show this message
  --input  The mask file used to create the regional mask file
  --output  Name of the output file
  --latmin  The minimum latitude (-90 to 90) of the desired region
  --latmax  The maximum latitude (-90 to 90) of the desired region
  --lonmin  The minimum longitude (0 to 360) of the desired region
  --lonmax  The maximum longitude (0 to 360) of the desired region
EOF
}


die() {
    printf '%s\n' "$1" >&2
    exit 1
}

# Initialize all the option variables.
# This ensures we are not contaminated by variables from the environment.
file=
verbose=0

while :; do
    case $1 in
        -h|-\?|--help)
            show_help    # Display a usage synopsis.
            exit
            ;;
        -i|--input)       # Takes an option argument; ensure it has been specified.
            if [ "$2" ]; then
                input=$2
                shift
            else
                die 'input is empty, using gswp3 mask file'
            fi
            ;;
        --input=?*)
            input=${1#*=} # Delete everything up to "=" and assign the remainder.
            ;;
        --input=)         # Handle the case of an empty --input=
            die 'input is empty, using gswp3 mask file'
            ;;

        -o|--output)       # Takes an option argument; ensure it has been specified.
            if [ "$2" ]; then
                output=$2
                shift
            else
                echo 'output is empty setting to regional_'${input}
                output='regional_'${input}
            fi
            ;;
        --output=?*)
            output=${1#*=} # Delete everything up to "=" and assign the remainder.
            ;;
        --output=)         # Handle the case of an empty --output=
             echo 'output is empty setting to regional_'${input}
             output='regional_'${input}
            ;;
        --latmin)       # Takes an option argument; ensure it has been specified.
            if [ "$2" ]; then
                latmin=$2
                shift
            else
                die 'latmin must be set'
            fi
            ;;
        --latmin=?*)
            latmin=${1#*=} # Delete everything up to "=" and assign the remainder.
            ;;
        --latmin=)         # Handle the case of an empty --latmin=
            die 'latmin must be set'
            ;;
        --latmax)       # Takes an option argument; ensure it has been specified.
            if [ "$2" ]; then
                latmax=$2
                shift
            else
                die 'latmax must be set'
            fi
            ;;
        --latmax=?*)
            latmax=${1#*=} # Delete everything up to "=" and assign the remainder.
            ;;
        --latmax=)         # Handle the case of an empty --latmax=
            die 'latmax must be set'
            ;;
        --lonmin)       # Takes an option argument; ensure it has been specified.
            if [ "$2" ]; then
                lonmin=$2
                shift
            else
                die 'lonmin must be set'
            fi
            ;;
        --lonmin=?*)
            lonmin=${1#*=} # Delete everything up to "=" and assign the remainder.
            ;;
        --lonmin=)         # Handle the case of an empty --lonmin=
            die 'lonmin must be set'
            ;;
        --lonmax)       # Takes an option argument; ensure it has been specified.
            if [ "$2" ]; then
                lonmax=$2
                shift
            else
                die 'lonmax must be set'
            fi
            ;;
        --lonmax=?*)
            lonmax=${1#*=} # Delete everything up to "=" and assign the remainder.
            ;;
        --lonmax=)         # Handle the case of an empty --lonmax=
            die 'lonmax must be set'
            ;;
        -v|--verbose)
            verbose=$((verbose + 1))  # Each -v adds 1 to verbosity.
            ;;
        --)              # End of all options.
            shift
            break
            ;;
        -?*)
            printf 'WARN: Unknown option (ignored): %s\n' "$1" >&2
            ;;
        *)               # Default case: No more options, so break out of the loop.
            break
    esac

    shift
done

# if --file was provided, open it for writing, else duplicate stdout
echo "Using "${input}" to create "${output}
echo "the bounding box extends from "${latmin}" to "${latmax}" degrees north"
echo "the bounding box extends from "${lonmin}" to "${lonmax}" degrees east"

cat reg_mask.nco > my_regmask.nco
#put values in nco script
sed -i "s/latmin/$latmin/g" my_regmask.nco
sed -i "s/latmax/$latmax/g" my_regmask.nco
sed -i "s/lonmin/$lonmin/g" my_regmask.nco
sed -i "s/lonmax/$lonmax/g" my_regmask.nco

module load nco

ncap2 -O -s 'lat2D[lat,lon]=lat+lon*0.0f' $input tmp.nc 
ncap2 -A -s 'lon2D[lat,lon]=lat*0.0f+lon' tmp.nc tmp.nc
ncap2 -A -S my_regmask.nco tmp.nc tmp.nc
ncks -O -v lat,lon,landsea tmp.nc $output
rm tmp.nc

#ncap2 -O -v -s "where(lat2D < $latmin || lat2D > $latmax || lon2D < $lonmin || lon2D > $lonmax) landsea=1" tmp.nc $output
#
#
#mv tmp_script.nco regional_mask_generator_lats_"$latmin"_"$latmax"_lons_"$lonmin"_"$lonmax".nco
# Rest of the program here.
# If there are input files (for example) that follow the options, they
# will remain in the "$@" positional parameters.

