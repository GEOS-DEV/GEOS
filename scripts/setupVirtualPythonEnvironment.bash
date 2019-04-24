#!/bin/bash

# Read input arguments
PYTHON_ROOT=/usr/tce/packages/python/python-3.6.4
VIRTUAL_NAME=geosx
OUTPUT_PATH=$HOME/Python/virtual
MINICONDA_BUILD=""


while [[ $# > 0 ]]
do
key="$1"

case $key in
    -p|--python_root)
    PYTHON_ROOT="$2"
    shift # past argument
    ;;
    -o|--output_path)
    OUTPUT_PATH="$2"
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
    echo "-p/--python_root \"Root path of the parent python environment (default = system python3.6.4 on LC)\""
    echo "-o/--output_path \"Path to store the new python environment (default = ~/Python/virtual)\""
    echo "-v/--virtual_name \"Virtual environment name (default = geosx)\""
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
    $PYTHON_ROOT/bin/pip install virtualenv
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
pip install h5py lxml pyevtk
pip install $DIR/../src/coreComponents/python/modules/pygeos_package


# Print user-info
echo "Complete!"
echo ""
echo "Notes:"
echo "To load the virtual environment, run:"
echo "  source $OUTPUT_PATH/start_$VIRTUAL_NAME.sh"
echo "To exit the environent, run:"
echo "  deactivate"


