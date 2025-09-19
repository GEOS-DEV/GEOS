#!/bin/bash

# nothing to do if schema file not given
if [ -z "$1" ]; then
    echo "Usage: $0 <schema> [<path>...]"
    exit
fi

SCHEMA=$1; shift
LOGFILE=$(pwd)/xml_validation_results.log

# "-r" in GNU xargs omits the call if input is empty
# OS X xargs does not support it, but does the same by default
if [ "$(uname)" == "Darwin" ]; then
    XARGS="xargs"
else
    XARGS="xargs -r"
fi

# check if xmllint is present
if ! hash xmllint &> /dev/null; then
    >&2 echo "Error: xmllint is required to validate xml files"
    exit
fi

abs_path ()
{
    if [ "$#" -gt 0 ]; then
        realpath "$@"
    fi
}

list_xml_files_at_path() 
{
    abs_path $(find . -type f -name "*.xml" -not -path "*/\.*")
}

# create/nullify the log file
echo -n > $LOGFILE
# validate each path separately and write results in the log
for path in "$@"; do
    # emit location and check directory
    echo "Validating in directory: $path"
    if [ ! -d "$path" ]; then
        echo "Directory not found: $path" >&2
        exit 1
    fi
    
    cd "$path" || continue
    collected_xml_files=$(list_xml_files_at_path)

    if [ -z "$collected_xml_files" ]; then
        echo "No XML files found in: $path"
        continue
    fi
    
    for xml_file in $collected_xml_files; do
        xmllint --schema $SCHEMA --noout "$xml_file" >> $LOGFILE 2>&1
        if [ $? -ne 0 ]; then
            echo "Validation failed for $xml_file" >> $LOGFILE
        fi
    done    
done


# print any failed validations on the stderr
grep -v validates $LOGFILE >&2

# if there are failed validations, message and return the status
if grep -q -v validates $LOGFILE; then
    >&2 echo "XML validation failed. See details in $LOGFILE"
    exit 1
fi

