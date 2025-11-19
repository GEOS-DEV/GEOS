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
FILES=$(eval "$FIND_CMD" 2>/dev/null || echo "");
VIOLATIONS_FOUND=0

ARRAY_MAP=()
ARRAY_UMAP=()
ARRAY_VECTOR=()
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

  if grep -q "std::map\s*<" "$file" ; then
    STR_MAP="  Found forbidden std::map usage in: $file"$'\n  '
    STR_MAP+= grep -n "std::map\s*<" "$file"
    ARRAY_MAP+="$STR_MAP"
    VIOLATIONS_FOUND=1
  fi
  if grep -q "std::unordered_map\s*<" "$file" ; then
    STR_UMAP="  Found forbidden std::unordered_map usage in: $file"$'\n  '
    STR_UMAP+= grep -n "std::unordered_map\s*<" "$file"
    ARRAY_UMAP+="$STR_UMAP"
    VIOLATIONS_FOUND=2
  fi
  if grep -q "std::vector\s*<" "$file" ; then
    STR_VECTOR="  Found forbidden std::vector usage in: $file"$'\n  '
    STR_VECTOR+= grep -n "std::vector\s*<" "$file"
    ARRAY_VECTOR+="$STR_VECTOR"
    VIOLATIONS_FOUND=3
  fi
done

if (($VIOLATIONS_FOUND != 0)); then 
  echo "----------------------------------------"
  echo "SUMMARY: Code rule violations found"
  echo "----------------------------------------"

  if((VIOLATIONS_FOUND == 1 ));then
    echo $'\nERROR: Forbidden std::map usage detected\n'
  fi

  for element in "${ARRAY_MAP[@]}"
  do
      echo "$element";
  done

  if((VIOLATIONS_FOUND == 2 ));then
    echo $'\nERROR: Forbidden std::unordered_map usage detected\n'
  fi

  for element in "${ARRAY_UMAP[@]}"
  do
      echo "$element";
  done

  if((VIOLATIONS_FOUND == 3 ));then
    echo $'\nERROR: Forbidden std::vector usage detected\n'
  fi

  for element in "${ARRAY_VECTOR[@]}"
  do
      echo "$element";
  done

  echo ""
  exit 1;
fi

echo "No code rule violations found"
exit 0
