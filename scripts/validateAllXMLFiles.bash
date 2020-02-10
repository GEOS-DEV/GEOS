#!bin/bash

SRC_PATH=$1
VALIDATION_RESULTS=xml_validation_results.log

find $SRC_PATH/../examples -name "*.xml" | xargs xmllint --schema $SRC_PATH/coreComponents/fileIO/schema/schema.xsd --noout &> $VALIDATION_RESULTS
find $SRC_PATH -name "*.xml" | xargs xmllint --schema $SRC_PATH/coreComponents/fileIO/schema/schema.xsd --noout &>> $VALIDATION_RESULTS
grep -v validates $VALIDATION_RESULTS >&2

if grep -q -v validates $VALIDATION_RESULTS
then
exit 1
fi


