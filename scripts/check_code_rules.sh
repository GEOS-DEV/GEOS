#!/bin/bash

echo "Test code rules"
FILE_PATTERNS=(
          "src/coreComponents/codingUtilities/*"
          "src/coreComponents/dataRepository/*"
          "src/coreComponents/denseLinearAlgebra/*"
          "src/coreComponents/discretizationMethods/*"
          "src/coreComponents/events/*"
          "src/coreComponents/fieldSpecification/*"
          "src/coreComponents/fileIO/*"
          "src/coreComponents/finiteElement/*"
          "src/coreComponents/finiteVolume/*"
          "src/coreComponents/functions/*"
          "src/coreComponents/linearAlgebra/*"
          "src/coreComponents/mainInterface/*"
          "src/coreComponents/mesh/*"
          "src/coreComponents/physicsSolvers/*"
        )
  EXCLUDE_PATTERNS=(
        "src/coreComponents/common/Datatype.hpp"
        "src/coreComponents/common/StdContaienrWrappers.hpp"
      )
FIND_CMD="find"
for pattern in "${FILE_PATTERNS[@]}"; do
  FILE_PATH_PATERN+=${pattern}" "
done
FIND_CMD="$FIND_CMD $FILE_PATH_PATERN"' \( -name "*.hpp" -o -name "*.cpp" \)'
echo "FIND_CMD ${FIND_CMD} " 
VIOLATIONS_FOUND=0
FILES=$(eval "$FIND_CMD" 2>/dev/null || echo "");

ARRAY=()
for file in $FILES; do
  SKIP=0
  for exclude in "${EXCLUDE_PATTERNS[@]}"; do
    if [[ "$file" == *"$exclude"* ]]; then
      SKIP=1
      break
    fi
  done
  
  if [ $SKIP -eq 1 ]; then
    continue
  fi
  if [grep -n "std::map\s*<" "$file"] ; then
    STR1="Found forbidden std::map usage in: $file"$'\n'
    STR1+= grep -n "std::map\s*<" "$file" 
    # while read line; do
    # STR1+="   Line: $line"
    # echo "$STR1";
    ARRAY+="$STR1"
    VIOLATIONS_FOUND=1
  fi
done

for element in "${ARRAY[@]}"
do
    echo "$element";
done

if (($VIOLATIONS_FOUND == 1)); then 
  exit 1;
fi
exit 0
