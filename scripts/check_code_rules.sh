#!/bin/bash


################################
## I. Function initializations
################################

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

################################
## II. GLOBAL INITIALIZATION  
################################

declare -A CONTAINER_MAP
CONTAINER_MAP=(
    ["std::map"]="MAP"
    ["std::unordered_map"]="UMAP"
    ["std::vector"]="VECTOR"
)

ARRAY_MAP=()
ARRAY_UMAP=()
ARRAY_VECTOR=()

MAP_VIOLATIONS_FOUND=0
UMAP_VIOLATIONS_FOUND=0
VECTOR_VIOLATIONS_FOUND=0

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
          "Datatype.hpp"     
          "StdContainerWrappers.hpp"     
          "BufferOps_inline.hpp"
          "BufferOps.hpp"
          "PVTPackage"
          "hdf5_interface"
    )

FILE_PATH_ARGS=()
for pattern in "${FILE_PATTERNS[@]}"; do
    if [ -d "${FILE_PREFIX}${pattern}" ]; then
      FILE_PATH_ARGS+=("${FILE_PREFIX}${pattern}")
    fi
done

EXCLUDE_EXPRESSION=()
for pattern in "${EXCLUDE_PATTERNS[@]}"; do
    if [[ ! "$pattern" == *".hpp"* ]]; then
      EXCLUDE_EXPRESSION+=( -path "*/$pattern" -prune -o )
    else
      EXCLUDE_EXPRESSION+=( -name "*${pattern}" -prune -o )
    fi
done

echo "======================================================"
echo -e "            TREE LIST OF FILES FOUND            "
echo "======================================================"

if [ ${#FILE_PATH_ARGS[@]} -gt 0 ]; then
    find "${FILE_PATH_ARGS[@]}" "${EXCLUDE_EXPRESSION[@]}" -type d -print  | sort | while IFS= read -r dir; do
        echo -e "->" $(echo "$dir" | sed 's/[]^$*.[]/\\&/g' | sed 's/\// |---/g')
    done
fi

ARRAY_FILES=()
# mapfile used for reading input lines into an array; -d $'\0': Specifies that the delimiter is (\0).
# -print0 : ask find to separate file paths by '\0'
mapfile -d $'\0' ARRAY_FILES < <(find "${FILE_PATH_ARGS[@]}" "${EXCLUDE_EXPRESSION[@]}" \
                                      -type f \( -name "*.hpp" -o -name "*.cpp" -o  -name "*.hpp.template" -o -name "*.cpp.template" \) \
                                      -print0 2>/dev/null)

################################
## III. MAIN LOOP
################################

for file in ${ARRAY_FILES[@]}; do
  for container_name in "${!CONTAINER_MAP[@]}"; do
    prefix="${CONTAINER_MAP[$container_name]}"
    var_violation="${prefix}_VIOLATIONS_FOUND"
    var_name_array="ARRAY_${prefix}"

    check_container_usage "$container_name" "$var_name_array" "$var_violation"
  done
done

# Print section
if [ $MAP_VIOLATIONS_FOUND -eq 1 ] || [ $UMAP_VIOLATIONS_FOUND -eq 1 ] || [ $VECTOR_VIOLATIONS_FOUND -eq 1 ]; then 
  echo "----------------------------------------"
  echo "SUMMARY: Code rule violations found"
  echo "----------------------------------------"
  for container_name in "${!CONTAINER_MAP[@]}"; do
    prefix="${CONTAINER_MAP[$container_name]}"
    var_violation="${prefix}_VIOLATIONS_FOUND"
    var_name_array="ARRAY_${prefix}"
    actual_count="${!var_violation}"

    print_violation "$actual_count" "$var_name_array" "$container_name"
  done
  echo ""
  exit 1;
fi

echo "No code rule violations found"
exit 0
