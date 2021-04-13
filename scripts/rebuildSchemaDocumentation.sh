
#!/bin/bash

# Read inputs
CONFIG=toss_3_x86_64_ib-clang@6.0.0-NoOPENMP

while [[ $# > 0 ]]
do
key="$1"

case $key in
    -c|--config)
    CONFIG="$2"
    shift # past argument
    ;;
    -?|--help)
    echo ""
    echo "Script options:"
    echo "-c/--config \"Target host config\""
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


# Find and move to the GEOSX root directory
SCRIPTS_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
cd $SCRIPTS_DIR/..


# Create a new GEOSX build
module load cmake/3.12.1
python scripts/config-build.py -bt Release -hc ./host-configs/$CONFIG.cmake
cd build-$CONFIG-release
make -j8
ls

# Generate an updated schema
cd bin
./geosx -s ../../src/coreComponents/schema/schema.xsd

# Build the documentation tables
cd ../../src/coreComponents/schema/
python SchemaToRSTDocumentation.py



