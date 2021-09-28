#!/bin/bash

# Parse and modify input arguments
INPUT=""
NEW_ARGS=""
SCRIPT_DIR=$(dirname "$0")

while [[ $# > 0 ]]
do
    key="$1"
    case $key in
        -i|--input)
        INPUT="$2"
        NEW_ARGS="$NEW_ARGS $key $2.preprocessed"
        shift
        ;;
        *)
        NEW_ARGS="$NEW_ARGS $key $2"
        shift
        ;;
    esac
    shift
done


# Preprocess the input file
FILE=$1     
if [ -f $SCRIPT_DIR/preprocess_xml ]; then
   $SCRIPT_DIR/preprocess_xml $INPUT -o $INPUT.preprocessed -s $SCRIPT_DIR/../../src/coreComponents/schema/schema.xsd
else
   echo "XML preprocessor not found"
   echo "To build it, run \"make geosx_xml_tools\""
fi


# Run the code
$SCRIPT_DIR/geosx $NEW_ARGS


