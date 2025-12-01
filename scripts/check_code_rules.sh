#!/bin/bash

# while
# ... done < <(grep -n "std::map\s*<" "$file") 
# Process the result of grep command as a file.
# This allows everything to be handled in the same shell environment.
check_container_usage() {
  local container="$1"
  local -n array_name="$2"
  local -n violation_found="$3"
  local str="  Found forbidden ${container} usage in: $file"$'\n'

  if grep -q "${container}\s*<" "$file"; then
    while IFS= read -r line; do
      str+="    $line"$'\n'
    done < <(grep -n "${container}\s*<" "$file")
    
    array_name+="$str"
    violation_found=1
  fi
}

print_violation()
{
    local violation_found="$1"
    local -n container="$2"
    local targetStd="$3"

    if [ "$violation_found" -eq 1 ];then
      echo "ERROR: Forbidden $targetStd usage detected"
      echo "=========================================="

      for element in "${container[@]}"
      do
          echo "$element";
      done
    fi
}

# Variables initialization 
FILE_PREFIX="src/coreComponents/"
FILE_PATTERNS=(
          "codingUtilities"
          "common"
          "dataRepository"
          "constitutive"
          "denseLinearAlgebra"
          "discretizationMethods"
          "events"
          "fieldSpecification"
          "fileIO"
          "finiteElement"
          "finiteVolume"
          "functions"
          "linearAlgebra"
          "mainInterface"
          "mesh"
          "physicsSolvers"
        )
  EXCLUDE_PATTERNS=(
        "common/Datatype.hpp"     
        "dataRepository/BufferOps_inline.hpp"
        "dataRepository/BufferOps.hpp"
      )


MAP_VIOLATIONS_FOUND=0
UMAP_VIOLATIONS_FOUND=0
VECTOR_VIOLATIONS_FOUND=0

ARRAY_MAP=()
ARRAY_UMAP=()
ARRAY_VECTOR=()

FILE_PATH_ARGS=()
for pattern in "${FILE_PATTERNS[@]}"; do
    if [ -d "${FILE_PREFIX}${pattern}" ]; then
      FILE_PATH_ARGS+=("${FILE_PREFIX}${pattern}")
    fi
done

FILES=$(find "${FILE_PATH_ARGS[@]}" -type f \( -name "*.hpp" -o -name "*.cpp" \) 2>/dev/null)

# Main loop
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

  check_container_usage "std::map" ARRAY_MAP MAP_VIOLATIONS_FOUND
  check_container_usage "std::unordered_map" ARRAY_UMAP UMAP_VIOLATIONS_FOUND
  check_container_usage "std::vector" ARRAY_VECTOR VECTOR_VIOLATIONS_FOUND
done

# Print section
if [ $MAP_VIOLATIONS_FOUND -eq 1 ] || [ $UMAP_VIOLATIONS_FOUND -eq 1 ] || [ $VECTOR_VIOLATIONS_FOUND -eq 1 ]; then 
  echo "----------------------------------------"
  echo "SUMMARY: Code rule violations found"
  echo "----------------------------------------"

  print_violation "$MAP_VIOLATIONS_FOUND" ARRAY_MAP "std::map"
  print_violation "$UMAP_VIOLATIONS_FOUND" ARRAY_UMAP "std::unordered_map"
  print_violation "$VECTOR_VIOLATIONS_FOUND" ARRAY_VECTOR "std::vector"

  echo ""
  exit 1;
fi

echo "No code rule violations found"
exit 0
