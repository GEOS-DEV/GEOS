#!/bin/bash


# Configuration
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
PACKAGE_DIR=$SCRIPT_DIR/../src/coreComponents/python/modules
declare -a TARGET_PACKAGES=("geosx_mesh_tools_package"
                            "geosx_xml_tools_package"
                            "hdf5_wrapper_package"
                            "pygeosx_tools_package")
declare -a LINK_SCRIPTS=("preprocess_xml"
                         "format_xml")


# Read input arguments
PYTHON_TARGET="$(which python)"
VIRTUAL_PATH=""
VIRTUAL_NAME="geosx"
MINICONDA_BUILD=""
INSTALL_VIRTUAL=false
BIN_DIR=""

while [[ $# > 0 ]]
do
key="$1"

case $key in
    -p|--python_target)
    PYTHON_TARGET="$2"
    shift # past argument
    ;;
    -m|--miniconda_build)
    MINICONDA_BUILD="$2"
    shift # past argument
    ;;
    -v|--virtual_path)
    VIRTUAL_PATH="$2"
    INSTALL_VIRTUAL=true
    shift # past argument
    ;;
    -b|--bin_dir)
    BIN_DIR="$2"
    shift # past argument
    ;;
    -?|--help)
    echo ""
    echo "Python environment setup options:"
    echo "-p/--python_target \"Target parent python bin\""
    echo "-m/--miniconda_build \"Fetch and build miniconda for the virtual environment (default = false)\""
    echo "-v/--virtual_path \"Path to store the new virtual environment\""
    echo "-b/--bin_dir \"Directory to link new scripts\""
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


# If requested, build a copy of miniconda
if [ -n "$MINICONDA_BUILD" ]
then
    echo "Building a copy of miniconda..."
    mkdir -p $MINICONDA_BUILD
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -b -p $MINICONDA_BUILD/miniconda
    rm Miniconda3-latest-Linux-x86_64.sh
    PYTHON_TARGET=$MINICONDA_BUILD/miniconda/bin/python
fi


# If a virtual environment is not explicitly requested,
# try installing packages directly
if ! $INSTALL_VIRTUAL
then
    echo "Attempting to install packages directly in target python environment..."
    for p in "${TARGET_PACKAGES[@]}"
    do
        echo "  $p"
        RES=$($PYTHON_TARGET -m pip install $PACKAGE_DIR/$p 2>&1)
        if [[ $RES =~ "Error" ]]
        then
            echo "  (cannot install target packes directly)"
            INSTALL_VIRTUAL=true
            break
        fi
    done
fi


# Virtual python installation method
if $INSTALL_VIRTUAL
then
    echo "Attempting to create a virtual python environment..."
    
    # Check to see if virtualenv is installed before continuing
    RES=$($PYTHON_TARGET -m pip list)
    if [[ ! $RES =~ "virtualenv" ]]
    then
        echo "Error: The parent python environment must have virtualenv installed"
        echo "       in order to setup a virtual environment!"
        exit 1
    fi

    # Set the default path if not already defined
    if [ -z "${VIRTUAL_PATH}" ]
    then
        VIRTUAL_PATH="$PWD/virtual_python_environment"
    fi

    # Setup the virtual environment
    mkdir -p $VIRTUAL_PATH/$VIRTUAL_NAME
    $PYTHON_TARGET -m virtualenv --system-site-packages $VIRTUAL_PATH/$VIRTUAL_NAME

    # Install packages
    echo "Installing packages..."
    PYTHON_TARGET=$VIRTUAL_PATH/$VIRTUAL_NAME/bin/python
    for p in "${TARGET_PACKAGES[@]}"
    do
        echo "  $p"
        $PYTHON_TARGET -m pip install $PACKAGE_DIR/$p
    done

    # Print user-info
    echo ""
    echo "Notes:"
    echo "To load the virtual environment, run:"
    echo "  source $VIRTUAL_PATH/$VIRTUAL_NAME/bin/activate"
    echo "To exit the environent, run:"
    echo "  deactivate"
fi


# Link key scripts to the bin directory
if [ -n "${BIN_DIR}" ]
then
    echo "Linking key scripts to bin directory..."
    MOD_PATH="$(dirname $PYTHON_TARGET)"

    for p in "${LINK_SCRIPTS[@]}"
    do
        echo "  $p"

        pp=
        if [ -f "$MOD_PATH/$p" ]
        then
            pp="$MOD_PATH/$p"
        else
            pp="$(which $p)"
        fi

        if [ -z "$pp" ]
        then
            echo "  (could not find where $p is installed)"      
        else
            echo "  (found $p as $pp)"
            ln -s $pp $BIN_DIR/$p 
        fi
    done

    ln -s $SCRIPT_DIR/automatic_xml_preprocess.sh $BIN_DIR/geosx_preprocessed
    ln -s $SCRIPT_DIR/pygeosx_preprocess.py $BIN_DIR/pygeosx_preprocess.py
fi

echo "Done!"

