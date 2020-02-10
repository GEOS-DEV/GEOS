#!bin/bash

VALIDATION_RESULTS=xml_validation_results.log

# find .. -name "*.xml" | xargs xmllint --schema ../src/coreComponents/fileIO/schema/schema.xsd --noout &> $VALIDATION_RESULTS
find ../examples -name "*.xml" | xargs xmllint --schema ../src/coreComponents/fileIO/schema/schema.xsd --noout &> $VALIDATION_RESULTS
find ../src -name "*.xml" | xargs xmllint --schema ../src/coreComponents/fileIO/schema/schema.xsd --noout &>> $VALIDATION_RESULTS
grep -v validates $VALIDATION_RESULTS >&2

if grep -q -v validates $VALIDATION_RESULTS
then
exit 1
fi


