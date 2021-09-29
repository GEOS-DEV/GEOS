#!/bin/bash

# Parse and modify input arguments
INPUT_ARGS=""
INPUT_COUNTER=0
OUTPUT_NAME=""
PARAMETER_ARGS=""
NEW_ARGS=""
SCRIPT_DIR=$(dirname "$0")

while [[ $# > 0 ]]
do
    key="$1"
    case $key in
        -i|--input)
        INPUT_ARGS="$INPUT_ARGS -i $2"
        INPUT_COUNTER=$(( INPUT_COUNTER + 1 ))
        OUTPUT_NAME=$2.preprocessed
        shift
        ;;
        -p|--parameter)
        echo $2
        PARAMETER_ARGS="$PARAMETER_ARGS -p $2 $3"
        shift
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
if [ -f $SCRIPT_DIR/preprocess_xml ]; then
   if [ "$INPUT_COUNTER" -gt "1" ]
   then
      OUTPUT_NAME="composite.xml.preprocessed"
   fi

   # Preprocess the file
   echo "Preprocessing xml:"
   $SCRIPT_DIR/preprocess_xml $INPUT_ARGS $PARAMETER_ARGS -o $OUTPUT_NAME -s $SCRIPT_DIR/../../src/coreComponents/schema/schema.xsd

   # Continue by running GEOSX
   echo "Running command:"
   echo "$SCRIPT_DIR/geosx $NEW_ARGS -i $OUTPUT_NAME"
   $SCRIPT_DIR/geosx $NEW_ARGS -i $OUTPUT_NAME
else
   echo "Error: XML preprocessor not found"
   echo "To build it, run \"make geosx_xml_tools\""
fi



