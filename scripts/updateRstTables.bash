#!bin/bash

# Generate an updated schema
bin/geosx -s ../src/coreComponents/fileIO/schema/schema.xsd

# Build the documentation tables
cd ../src/coreComponents/fileIO/schema/
python SchemaToRSTDocumentation.py
