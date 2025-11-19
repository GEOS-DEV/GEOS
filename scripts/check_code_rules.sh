echo "Test if the build works"
eval pwd

if [[ "${TEST_CODE_RULES}" = true ]]; then
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
    FILE_PATH_PATERN+=${GEOS_BUILD_DIR}/${pattern}" "
  done
  FIND_CMD="$FIND_CMD $FILE_PATH_PATERN"' \( -name "*.hpp" -o -name "*.cpp" \)'

  echo "FIND_CMD ${FIND_CMD} " 
  VIOLATIONS_FOUND=0
  FILES=$(eval "$FIND_CMD" 2>/dev/null || echo "");
      
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
    echo "file ${file} " 
    if grep -n "std::map\s*<" "$file" ; then
      echo "Found forbidden std::map usage in: $file"
      grep -n "std::map\s*<" "$file" | while read line; do
        echo "   Line: $line"
      done
      VIOLATIONS_FOUND=1
    fi
  done
  exit 0
fi