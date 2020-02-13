#!bin/bash

SRC_PATH=$1
VALIDATION_RESULTS=xml_validation_results.log

if hash xmllint &> /dev/null
then
  find $SRC_PATH/../examples -name "*.xml" | xargs xmllint --schema $SRC_PATH/coreComponents/fileIO/schema/schema.xsd --noout &> $VALIDATION_RESULTS
  find $SRC_PATH -name "*.xml" | xargs xmllint --schema $SRC_PATH/coreComponents/fileIO/schema/schema.xsd --noout &>> $VALIDATION_RESULTS
  grep -v validates $VALIDATION_RESULTS >&2

  if grep -q -v validates $VALIDATION_RESULTS
  then
    echo "xml validation failed... See details in $VALIDATION_RESULTS"
    exit 1
  fi
else
  echo "Note: xmllint is required to validate xml files in repository"
fi





