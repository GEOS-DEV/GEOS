#!/bin/bash

# Parse and modify input arguments
INPUT_ARGS=""
INPUT_COUNTER=0
COMPILED_XML_NAME=""
COMPILED_XML_NAME_OVERRIDE=""
PARAMETER_ARGS=""
NEW_ARGS=""
USE_PYGEOSX=1
SCRIPT_DIR=$(dirname "$0")
PYGEOSX=$SCRIPT_DIR/../lib/PYGEOSX/bin/python

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
        PARAMETER_ARGS="$PARAMETER_ARGS -p $2 $3"
        shift
        shift
        ;;
        -u|--use-pygeosx)
        USE_PYGEOSX=$2
        shift
        ;;
        -h|--help)
        echo ""
        echo "Preprocessor options:"
        echo "-i/--input          Input xml file name"
        echo "-p/--parameter      XML parameter overrides (name and value separated by a space)"
        echo "-u/--use-pygeosx    Use pygeosx for xml preprocessing (default=1)"
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


# Check for pygeosx
if [ "$USE_PYGEOSX" -eq "1" ]
then
   if [ -f $PYGEOSX ]
   then
      echo "Using pygeosx to preprocess the xml file"
   else
      echo "Pygeosx installation not found... reverting to non-pygeosx version"
      USE_PYGEOSX=0
   fi
fi


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


   if [ "$USE_PYGEOSX" -eq "1" ]
   then
      echo "Running command: $PYGEOSX $SCRIPT_DIR/pygeosx_preprocess.py $INPUT_ARGS $PARAMETER_ARGS $NEW_ARGS"
      $PYGEOSX $SCRIPT_DIR/pygeosx_preprocess.py $INPUT_ARGS $PARAMETER_ARGS $NEW_ARGS -s $SCRIPT_DIR/../../src/coreComponents/schema/schema.xsd
   else
      # Preprocess the file
      echo "Preprocessing xml: $INPUT_ARGS"
      $SCRIPT_DIR/preprocess_xml $INPUT_ARGS $PARAMETER_ARGS -o $COMPILED_XML_NAME -s $SCRIPT_DIR/../../src/coreComponents/schema/schema.xsd

      # Continue by running GEOSX
      echo "Running command: $SCRIPT_DIR/geosx -i $COMPILED_XML_NAME $NEW_ARGS"
      $SCRIPT_DIR/geosx -i $COMPILED_XML_NAME $NEW_ARGS
   fi
   
else
   echo "Error: XML preprocessor not found"
   echo "To build it, run \"make geosx_xml_tools\""
fi



