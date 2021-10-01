#!/bin/bash

# Parse and modify input arguments
INPUT_ARGS=""
INPUT_COUNTER=0
COMPILED_XML_NAME=""
COMPILED_XML_NAME_OVERRIDE=""
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
        COMPILED_XML_NAME=$2.preprocessed
        shift
        ;;
        -c|--compiled_name)
        COMPILED_XML_NAME_OVERRIDE=$2
        shift
        ;;
        -p|--parameter)
        echo $2
        PARAMETER_ARGS="$PARAMETER_ARGS -p $2 $3"
        shift
        shift
        ;;
        -h|--help)
        echo ""
        echo "Preprocessor options:"
        echo "-i/--input          Input xml file name"
        echo "-p/--parameter      XML parameter overrides (name and value separated by a space)"
        echo "-c/--compiled_name  Compiled xml filename (default=[input_name].preprocessed or composite.xml.preprocessed)"
        echo "Note: multiple -i and -p arguments are allowed"
        echo ""
        echo "GEOSX options:"
        $SCRIPT_DIR/geosx --help
        exit
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
      COMPILED_XML_NAME="composite.xml.preprocessed"
   fi

   if [ ! -z "$COMPILED_XML_NAME_OVERRIDE" ]
   then
      COMPILED_XML_NAME=$COMPILED_XML_NAME_OVERRIDE
   fi

   # Preprocess the file
   echo "Preprocessing xml: $INPUT_ARGS"
   $SCRIPT_DIR/preprocess_xml $INPUT_ARGS $PARAMETER_ARGS -o $COMPILED_XML_NAME -s $SCRIPT_DIR/../../src/coreComponents/schema/schema.xsd

   # Continue by running GEOSX
   echo "Running command: $SCRIPT_DIR/geosx -i $COMPILED_XML_NAME $NEW_ARGS"
   $SCRIPT_DIR/geosx -i $COMPILED_XML_NAME $NEW_ARGS
else
   echo "Error: XML preprocessor not found"
   echo "To build it, run \"make geosx_xml_tools\""
fi



