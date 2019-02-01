#!/bin/bash

# Read input arguments
PYTHON_ROOT=/usr/tce/packages/python/python-3.6.4
OUTPUT_PATH=$HOME/Python/virtual/geosx

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
    -?|--help)
    echo ""
    echo "Virtual Python environment setup options:"
    echo "-p/--python_root \"Root path of the parent python environment\""
    echo "-o/--output_path \"Path to place the new python environment\""
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
echo "Setting up a virtual python environment..."
mkdir -p $OUTPUT_PATH
$PYTHON_ROOT/bin/virtualenv --system-site-packages $OUTPUT_PATH
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# Install packages
source $OUTPUT_PATH/bin/activate
pip install h5py lxml pyevtk
pip install $DIR/../src/coreComponents/python/modules/pygeos_package

echo "Complete!"
