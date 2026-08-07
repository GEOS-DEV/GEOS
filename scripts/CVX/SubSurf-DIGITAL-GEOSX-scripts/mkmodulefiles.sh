#!/bin/bash 

iroot="/data/saet/mtml/software/x86_64/RHEL7/GEOSX/0.2.0"

for MFF in $(ls -1 $1 ) ; do

    MF=$(basename $MFF .cmake) 
    geosxexes="$(ls -1 ${iroot}/install-${MF}-*/bin/geosx 2>/dev/null)"
    for geosxexe in $geosxexes ; do 
	echo "# $0 : geosxexe=$geosxexe "
	if [ -x $geosxexe ]; then 

	    echo "# Generating module and init files for $geosxexe, $MF.cmake ($MFF)"

	    if [[ "$MF" =~ "CPU"  ]]; then
		if [[ "$MF" =~ "impi"  ]]; then 
		    $ECHO env INIT_ONLY=1 GEOSX=1 WHO=other ./GEOSX-CPU-GCC_10.2.0-impi_2021.06--cache-build-install.sh $MF
		elif [[  "$MF" =~ "hpcx"   ]]; then 
		    $ECHO env INIT_ONLY=1 GEOSX=1 WHO=other ./GEOSX-CPU-GCC_10.2.0-ompi_hpcx--cache-build-install.sh $MF
		fi
	    elif [[ "$MF" =~ "GPU"  ]]; then
		if [[ "$MF" =~ "impi"  ]]; then 
		    $ECHO env INIT_ONLY=1 GEOSX=1 WHO=other ./GEOSX-CPU-GCC_10.2.0-impi_2021.06--cache-build-install.sh $MF
		elif [[  "$MF" =~ "hpcx"   ]]; then 
		    $ECHO env INIT_ONLY=1 GEOSX=1 WHO=other ./GEOSX-GPU-GCC_10.2.0-ompi_hpcx--cache-build-install.sh $MF
		fi
	    fi
	else
	    echo "# geosxexe=$geosxexe does not exist ! "
	fi
    done
done 


## INITONLY=1 WHO=other sroot=${HOME}/src/new_test iroot=${HOME}/src/new_test BUILD_ONLY=1 TPL=1 GEOSX=1 ../GEOSX_mtml/GEOSX-CPU-GCC_10.2.0-impi_2021.06--cache-build-install.sh
