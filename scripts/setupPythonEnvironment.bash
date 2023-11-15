#!/bin/bash


# Configuration
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
PACKAGE_DIR=$SCRIPT_DIR/../src/coreComponents/python/modules
declare -a TARGET_PACKAGES=("$PACKAGE_DIR/geosx_mesh_tools_package"
                            "$PACKAGE_DIR/geosx_xml_tools_package"
                            "$PACKAGE_DIR/hdf5_wrapper_package"
                            "$PACKAGE_DIR/pygeosx_tools_package"
                            "$SCRIPT_DIR/../integratedTests/scripts/geos_ats_package")
declare -a LINK_SCRIPTS=("preprocess_xml"
                         "format_xml"
                         "convert_abaqus"
                         "run_geos_ats"
                         "setup_ats_environment"
                         "activate"
                         "python")


# Read input arguments
PYTHON_TARGET="$(which python3)"
VIRTUAL_PATH=""
VIRTUAL_NAME="geosx"
MINICONDA_BUILD=""
INSTALL_VIRTUAL=false
BIN_DIR=""
PIP_CMD="pip --disable-pip-version-check"


if [[ -z "${VERBOSE}" ]]
then
    VERBOSE=false
else
    VERBOSE=true
fi


while [[ $# > 0 ]]
do
key="$1"

case $key in
    -p|--python-target)
    PYTHON_TARGET="$2"
    shift # past argument
    ;;
    -m|--miniconda-build)
    MINICONDA_BUILD="$2"
    shift # past argument
    ;;
    -e|--environment-path)
    VIRTUAL_PATH="$2"
    INSTALL_VIRTUAL=true
    shift # past argument
    ;;
    -b|--bin-dir)
    BIN_DIR="$2"
    shift # past argument
    ;;
    -v|--verbose)
    VERBOSE=true
    shift # past argument
    ;;
    -?|--help)
    echo ""
    echo "Python environment setup options:"
    echo "-p/--python-target \"Target parent python bin\""
    echo "-m/--miniconda-build \"Fetch and build miniconda at the target directory (instead of targeting an existing version of python)\""
    echo "-e/--environment-path \"Path to store the new virtual environment\""
    echo "-b/--bin-dir \"Directory to link new scripts\""
    echo "-v/--verbose \"Increase verbosity level\""
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


# Check to make sure that the python target exists
echo "Checking the python target..."
if [ ! -f "$PYTHON_TARGET" ]
then
    echo "The target python executable ($PYTHON_TARGET) cannot be found"

    if [[ "$PYTHON_TARGET" == *"PYGEOSX"* ]]
    then
        echo "If GEOSX is configured to use pygeosx, you may need to run \"make pygeosx\""
        echo "before setting up the geosx_python_tools!"
    else
        echo "Note: if you use the \"-m\" argument, this script will try to build a version of miniconda"
        echo "that can support the target tools"
    fi
    exit 1
fi


# If a virtual environment is not explicitly requested,
# try installing packages directly
if ! $INSTALL_VIRTUAL
then
    echo "Attempting to install packages directly in target python environment..."
    for p in "${TARGET_PACKAGES[@]}"
    do
        if [ -d "$p" ]
        then
            echo "  $p"
            RES=""
            if $VERBOSE
            then
                RES=$($PYTHON_TARGET -m $PIP_CMD install $p)
                echo $RES
            else
                RES=$($PYTHON_TARGET -m $PIP_CMD install $p 2>&1)
            fi
            
            if [[ $RES =~ "HTTP error" ]]
            then
                echo "Could not access package dependencies:"
                echo $RES
                exit 1
            fi

            if [[ $RES =~ "Error" ]]
            then
                echo "Failed to install target package... attempting to do so in a virtual environment"
                INSTALL_VIRTUAL=true
                break
            fi
        else
            echo "Could not find target package: $p"
        fi
    done
fi


# Virtual python installation method
if $INSTALL_VIRTUAL
then
    echo "Attempting to create a virtual python environment..."
    
    # Check to see if virtualenv is installed before continuing
    RES=$($PYTHON_TARGET -m $PIP_CMD list)
    if $VERBOSE
    then
        echo "Available packages in base python environment:"
        echo $RES
    fi

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
    if $VERBOSE
    then
        echo "Virtual environment path:"
        echo $VIRTUAL_PATH
    fi

    mkdir -p $VIRTUAL_PATH/$VIRTUAL_NAME
    $PYTHON_TARGET -m virtualenv --system-site-packages $VIRTUAL_PATH/$VIRTUAL_NAME

    # Install packages
    echo "Installing packages..."
    PYTHON_TARGET=$VIRTUAL_PATH/$VIRTUAL_NAME/bin/python
    for p in "${TARGET_PACKAGES[@]}"
    do
        if [ -d "$p" ]
        then
            echo "  $p"
            RES=""
            if $VERBOSE
            then
                RES=$($PYTHON_TARGET -m $PIP_CMD install $p)
                echo $RES
            else
                RES=$($PYTHON_TARGET -m $PIP_CMD install $p 2>&1)
            fi

            if [[ $RES =~ "HTTP error" ]]
            then
                echo "Could not access package dependencies"
                echo $RES
                exit 1
            fi
            
            if [[ $RES =~ "Error" ]]
            then
                echo "  (failed to install target package)"
                break
            fi
        else
            echo "Could not find target package: $p"
        fi
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
declare -a MOD_SEARCH_PATH=("$(dirname $PYTHON_TARGET)"
                            "$HOME/.local/bin"
                            "$HOME/local/bin")


if [ -n "${BIN_DIR}" ]
then
    echo "Linking key scripts to bin directory..."

    for p in "${LINK_SCRIPTS[@]}"
    do
        echo "  $p"
        if $VERBOSE
        then
            echo "    searching the following paths:"
        fi
        package_found="0"

        for MOD_PATH in "${MOD_SEARCH_PATH[@]}"
        do
            # Check to see if the tool exists
            pp=
            if $VERBOSE
            then
                echo "      $MOD_PATH/$p"
            fi
            
            if [ -f "$MOD_PATH/$p" ]
            then
                pp="$MOD_PATH/$p"
            fi

            # Remove any old links if necessary
            if [ -f "$BIN_DIR/$p" ]
            then
                rm $BIN_DIR/$p
            fi

            # Create links
            if [ ! -z "$pp" ]
            then
                if $VERBOSE
                then
                    echo "    (found $p as $pp)"
                fi
                ln -s $pp $BIN_DIR/$p 
                package_found="1"
                break
            fi
        done

        if [[ "$package_found" == "0" ]]
        then
            echo "    (could not find where $p is installed)" 
        fi
    done

    # Link additional tools from the scripts directory
    echo "Linking additional scripts to the bin directory..."
    if [ ! -f "$BIN_DIR/geosx_preprocessed" ]
    then
        ln -s $SCRIPT_DIR/automatic_xml_preprocess.sh $BIN_DIR/geosx_preprocessed
        ln -s $SCRIPT_DIR/pygeosx_preprocess.py $BIN_DIR/pygeosx_preprocess.py
    fi
fi

echo "Done!"

