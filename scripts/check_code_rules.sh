#!/bin/bash

echo "Test code rules"
FILE_PATTERNS=(
          "src/coreComponents/codingUtilities/*"
          "src/coreComponents/common/*"
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
        "src/coreComponents/common/StdContainerWrappers.hpp"
        "src/coreComponents/dataRepository/BufferOps_inline.hpp"
        "src/coreComponents/dataRepository/BufferOps.hpp"
      )
FIND_CMD="find"
for pattern in "${FILE_PATTERNS[@]}"; do
  FILE_PATH_PATERN+=${pattern}" "
done
FIND_CMD="$FIND_CMD $FILE_PATH_PATERN"' \( -name "*.hpp" -o -name "*.cpp" \)'
FILES=$(eval "$FIND_CMD" 2>/dev/null || echo "");

MAP_VIOLATIONS_FOUND=0
UMAP_VIOLATIONS_FOUND=0
VECTOR_VIOLATIONS_FOUND=0

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
    STR_MAP="  Found forbidden std::map usage in: $file"$'\n'
    while IFS= read -r line; do
      STR_MAP+="    $line"$'\n'
    done < <(grep -n "std::map\s*<" "$file")
    ARRAY_MAP+=("$STR_MAP")
    MAP_VIOLATIONS_FOUND=1
  fi
  if grep -q "std::unordered_map\s*<" "$file" ; then
    STR_UMAP="  Found forbidden std::unordered_map usage in: $file"$'\n'
    while IFS= read -r line; do
      STR_UMAP+="    $line"$'\n'
    done < <(grep -n "std::unordered_map\s*<" "$file")
    ARRAY_UMAP+=("$STR_UMAP")
    UMAP_VIOLATIONS_FOUND=1
  fi
  if grep -q "std::vector\s*<" "$file" ; then
    STR_VECTOR="  Found forbidden std::vector usage in: $file"$'\n'
    while IFS= read -r line; do
      STR_VECTOR+="    $line"$'\n'
    done < <(grep -n "std::vector\s*<" "$file")
    ARRAY_VECTOR+=("$STR_VECTOR")
    VECTOR_VIOLATIONS_FOUND=1
  fi
  done

if (($MAP_VIOLATIONS_FOUND == 1)) || (($UMAP_VIOLATIONS_FOUND == 1 )) || (($VECTOR_VIOLATIONS_FOUND == 1 )); then 
  echo "----------------------------------------"
  echo "SUMMARY: Code rule violations found"
  echo "----------------------------------------"

  if((MAP_VIOLATIONS_FOUND == 1 ));then
    echo $'ERROR: Forbidden std::map usage detected'
    echo "=========================================="
  fi

  for element in "${ARRAY_MAP[@]}"
  do
      echo "$element";
  done

  if((UMAP_VIOLATIONS_FOUND == 1 ));then
    echo $'\nERROR: Forbidden std::unordered_map usage detected'
    echo "======================================================"
  fi

  for element in "${ARRAY_UMAP[@]}"
  do
      echo "$element";
  done

  if((VECTOR_VIOLATIONS_FOUND == 1 ));then
    echo $'\nERROR: Forbidden std::vector usage detected'
    echo "==============================================="
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
