#!/bin/bash


# Configuration
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
PYTHON_TARGET=
BIN_DIR=
PACKAGE_DIR=
TMP_CLONE_DIR=
declare -a PIP_INSTALL_ARGS=(--disable-pip-version-check install --upgrade)
PACKAGE_BRANCH=main


declare -a TARGET_PACKAGES=("geos-utils"
                            "geos-mesh"
                            "geos-xml-tools"
                            "hdf5-wrapper"
                            "pygeos-tools"
                            "geos-ats")
declare -a LINK_SCRIPTS=("preprocess_xml"
                         "format_xml"
                         "convert_abaqus"
                         "run_geos_ats"
                         "setup_ats_environment"
                         "geos_ats_log_check"
                         "geos_ats_restart_check"
                         "geos_ats_curve_check"
                         "geos_ats_process_tests_fails"
                         "mesh-doctor"
                         "activate"
                         "python")

function bootstrap_pip () {
    local python_executable="$1"

    if [[ "${GEOS_SKIP_PIP_BOOTSTRAP:-false}" == "true" ]]; then
        echo "Skipping pip bootstrap because GEOS_SKIP_PIP_BOOTSTRAP=true."
        return 0
    fi

    echo "Updating pip"
    if ! "${python_executable}" -m pip install --upgrade pip setuptools wheel; then
        echo "WARNING: pip bootstrap failed for ${python_executable}; continuing with the existing pip so package installation can report the actionable failure."
        "${python_executable}" -m pip --version || true
        return 0
    fi
}


# Read input arguments
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
    -b|--bin-dir)
    BIN_DIR="$2"
    shift # past argument
    ;;
    -d|--pkg-dir)
    PACKAGE_DIR="$2"
    shift # past argument
    ;;
    -r|--python-pkg-branch)
    PACKAGE_BRANCH="$2"
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
    echo "-b/--bin-dir \"Directory to link new scripts\""
    echo "-d/--pkg-dir \"Directory containing target python packages\""
    echo "-r/--python-pkg-branch \"Target branch for geosPythonPackages (default=main) \""
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


# Check to make sure that the python target exists
echo "Checking the python target..."
if [[ -z "${PYTHON_TARGET}" ]]
then
    echo "To setup the python environment, please specify the python executable path"
    exit 1
fi

if [ ! -f "$PYTHON_TARGET" ]
then
    echo "The target python executable ($PYTHON_TARGET) cannot be found"
    exit 1
fi

# Ubuntu 24.04+ marks the system interpreter as EXTERNALLY-MANAGED (PEP 668).
# geosx_python_tools_test already expects packages in ${CMAKE_BINARY_DIR}/python/geosx.
if "${PYTHON_TARGET}" -c 'import os, sysconfig; raise SystemExit(0 if os.path.isfile(os.path.join(sysconfig.get_path("stdlib"), "EXTERNALLY-MANAGED")) else 1)'
then
    if [[ -z "${BIN_DIR}" ]]
    then
        echo "The target python (${PYTHON_TARGET}) is externally managed (PEP 668)."
        echo "Pass -b/--bin-dir so a virtual environment can be created."
        exit 1
    fi

    VENV_DIR="$(cd "$(dirname "${BIN_DIR}")" && pwd)/python/geosx"
    echo "The target python (${PYTHON_TARGET}) is externally managed (PEP 668)."
    echo "Creating a virtual environment at ${VENV_DIR}"
    mkdir -p "$(dirname "${VENV_DIR}")"
    if ! "${PYTHON_TARGET}" -m venv "${VENV_DIR}"
    then
        echo "Failed to create a virtual environment at ${VENV_DIR}."
        echo "On Debian/Ubuntu, install the python3-venv package and retry."
        exit 1
    fi
    PYTHON_TARGET="${VENV_DIR}/bin/python"
fi


# Check for a predefined package directory
echo "Checking for python packages..."
if [[ -z "${PACKAGE_DIR}" ]]
then
    echo "Cloning the GEOS python package repository (branch=$PACKAGE_BRANCH)..."
    TMP_CLONE_DIR=$(mktemp -d)
    PACKAGE_DIR=$TMP_CLONE_DIR/geosPythonPackages
    git clone --depth 1 --branch $PACKAGE_BRANCH --single-branch https://github.com/GEOS-DEV/geosPythonPackages.git $PACKAGE_DIR
elif [ ! -d "${PACKAGE_DIR}/geosx_xml_tools_package" ]
then
    echo "The specified package directory does not contain the expected targets."
    echo "The path specified with -d/--pkg-dir should point to a copy of the geosPythonPackages repository."
    exit 1
fi


# Updating pip
bootstrap_pip "$PYTHON_TARGET"

# Install packages
echo "Installing python packages..."
for p in "${TARGET_PACKAGES[@]}"
do
    if [ -d "$PACKAGE_DIR/$p" ]
    then
        echo "  $p"

        # Try installing the package
        if [[ "${VERBOSE}" == true ]]
        then
            INSTALL_MSG=$("${PYTHON_TARGET}" -m pip "${PIP_INSTALL_ARGS[@]}" "${PACKAGE_DIR}/${p}")
            INSTALL_RC=$?
        else
            INSTALL_MSG=$("${PYTHON_TARGET}" -m pip "${PIP_INSTALL_ARGS[@]}" "${PACKAGE_DIR}/${p}" 2>&1)
            INSTALL_RC=$?
        fi

        if [ $INSTALL_RC -ne 0 ]
        then
            echo $INSTALL_MSG

            if [[ $INSTALL_MSG == *SSLError* || $INSTALL_MSG == *"[CRYPTO]"* ]]
            then
                echo "The setup script failed while pip was using HTTPS."
                echo "On FIPS-enabled hosts running non-FIPS containers, this can indicate an OpenSSL/FIPS policy mismatch."
                echo "Check the earlier Crypto/FIPS summary in the CI log."
            elif [[ $INSTALL_MSG =~ "HTTP error" ]]
            then
                echo "The setup script may have failed to fetch external dependencies"
                echo "Try re-running the \"make ats_environment\" script again on a machine that can access github"
            else
                echo "Failed to install packages"
                echo "Note: if you do not have write access for your target python environment,"
                echo "consider using a virtual python environment.  See these instructions for details:"
                echo "https://docs.python.org/3/library/venv.html"
            fi

            exit $INSTALL_RC
        fi
    else
        echo "Could not find target package: $p"
    fi
done


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


if [[ ! -z "${TMP_CLONE_DIR}" ]]
then
    rm -rf $TMP_CLONE_DIR
fi


echo "Done!"
