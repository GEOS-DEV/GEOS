#!/bin/bash

# Read input arguments
PYTHON_ROOT=/usr/tce/packages/python/python-3.6.4
VIRTUAL_NAME=geosx
OUTPUT_PATH=$HOME/Python/virtual
MINICONDA_BUILD=""
GEOSX_BUILD=""


while [[ $# > 0 ]]
do
key="$1"

case $key in
    -p|--python_target)
    echo $2
    PYTHON_ROOT="$(dirname $2)/.."
    shift # past argument
    ;;
    -o|--output_path)
    OUTPUT_PATH="$2"
    shift # past argument
    ;;
    -g|--geosx_build)
    GEOSX_BUILD="$2"
    shift # past argument
    ;;
    -m|--miniconda_build)
    MINICONDA_BUILD="$2"
    shift # past argument
    ;;
    -v|--virtual_name)
    VIRTUAL_NAME="$2"
    shift # past argument
    ;;
    -?|--help)
    echo ""
    echo "Virtual Python environment setup options:"
    echo "-m/--miniconda_build \"Fetch and build miniconda for the virtual environment (default = false)\""
    echo "-p/--python_target \"Target parent python (default = system python3.6.4 on LC)\""
    echo "-o/--output_path \"Path to store the new python environment (default = ~/Python/virtual)\""
    echo "-v/--virtual_name \"Virtual environment name (default = geosx)\""
    echo "-g/--geosx_build \"Path of GEOSX build directory for linking console scripts\""
    echo ""
    exit
    ;;
    --default)
    DEFAULT=YES
    ;;
    *)
            # unknown option
    ;;
esac
shift # past argument or value
done


# Setup the virtual environment
mkdir -p $OUTPUT_PATH/$VIRTUAL_NAME


# If requested, build a copy of miniconda
if [ -n "$MINICONDA_BUILD" ]
then
    mkdir -p $MINICONDA_BUILD

    # Install miniconda
    echo "Building a copy of miniconda..."
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -b -p $MINICONDA_BUILD/miniconda
    rm Miniconda3-latest-Linux-x86_64.sh
    PYTHON_ROOT=$MINICONDA_BUILD/miniconda

    # Install virtualenv
    $PYTHON_ROOT/bin/pip install numpy lxml h5py virtualenv
else
    echo "Using an existing python installation..."

    if [ ! -f "$PYTHON_ROOT/bin/virtualenv" ]
    then
      echo "Error: the target python environment must have virtualenv installed!"
      exit 1
    fi
fi


# Create a virtual environment
echo "Creating a virtual python environment..."
$PYTHON_ROOT/bin/virtualenv --system-site-packages $OUTPUT_PATH/$VIRTUAL_NAME


# Write a simple startup script
STARTUP_SCRIPT=$OUTPUT_PATH/start_$VIRTUAL_NAME.sh
echo "Writing a startup script to $STARTUP_SCRIPT"
echo "#!/bin/bash" > $STARTUP_SCRIPT
echo "source $OUTPUT_PATH/$VIRTUAL_NAME/bin/activate" >> $STARTUP_SCRIPT


# Install packages
echo "Installing packages..."
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

source $STARTUP_SCRIPT
pip install $DIR/../src/coreComponents/python/modules/geosx_xml_tools_package
pip install $DIR/../src/coreComponents/python/modules/hdf5_wrapper_package


# Link scripts in the GEOSX build dir
if [ -n "$GEOSX_BUILD" ]
then
    ln -s $OUTPUT_PATH/$VIRTUAL_NAME/bin/preprocess_xml $GEOSX_BUILD/bin/preprocess_xml
    ln -s $OUTPUT_PATH/$VIRTUAL_NAME/bin/format_xml $GEOSX_BUILD/bin/format_xml

    if [ -f "$GEOSX_BUILD/lib/pygeosx.so" ]
    then
        # SITE_PACKAGES=$OUTPUT_PATH/$VIRTUAL_NAME/lib/*/site-packages
        SITE_PACKAGES=$(find $OUTPUT_PATH/$VIRTUAL_NAME/lib/ -mindepth 1 -maxdepth 1 -type d)/site-packages
        ln -s $GEOSX_BUILD/lib/pygeosx.so $SITE_PACKAGES/pygeosx.so
        ln -s $GEOSX_BUILD/lib/libgeosx_core.so $SITE_PACKAGES/libgeosx_core.so
        ln -s $GEOSX_BUILD/lib/pylvarray.so $SITE_PACKAGES/pylvarray.so
        # ln -s $GEOSX_BUILD/lib/libgeosx_core.so $SITE_PACKAGES/libgeosx_core.so
    fi
fi


# Print user-info
echo "Complete!"
echo ""
echo "Notes:"
echo "To load the virtual environment, run:"
echo "  source $OUTPUT_PATH/start_$VIRTUAL_NAME.sh"
echo "To exit the environent, run:"
echo "  deactivate"


