#!/bin/bash


################################
## I. Function initializations
################################

# while
# ... done < <(grep -n "std::map\s*<" "$file") 
# Process the result of grep command as a file.
# This allows everything to be handled in the same shell environment.
check_container_usage() {
  local container_name="$1"
  local str="  Found forbidden ${container_name} usage in: $file"$'\n'

  if grep -q "${container_name}\s*<" "$file"; then
    while IFS= read -r line; do
      str+="    $line"$'\n'
    done < <(grep -n "${container_name}\s*<" "$file")

    ERRORS_CONTAINER[$container_name]+="$str"
    ((FORBIDDEN_CONTAINER_MAP["$container_name"]++))
  fi
}

print_violation()
{
    local container_name="$1"

    if [ "${FORBIDDEN_CONTAINER_MAP[$container_name]}" -eq 1 ];then
      echo $'\n'
      echo "ERROR: Forbidden $container_name usage detected"
      echo "=========================================="

      printf '%s' "${ERRORS_CONTAINER[$container_name]}"

    fi
}

################################
## II. GLOBAL INITIALIZATION  
################################

FORBIDDEN_EXPRESSIONS=(
          "std::map"     
          "std::unordered_map"
          "std::vector"
)

declare -A FORBIDDEN_CONTAINER_MAP=()
for exp in "${FORBIDDEN_EXPRESSIONS[@]}"; do
    FORBIDDEN_CONTAINER_MAP["$exp"]=0
done

declare -A ERRORS_CONTAINER=()
for exp in "${FORBIDDEN_EXPRESSIONS[@]}"; do
    ERRORS_CONTAINER["$exp"]=""
done

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

EXCLUDED_NAME_PATTERNS=()
for pattern in "${EXCLUDE_PATTERNS[@]}"; do
    if [[ ! "$pattern" == *".hpp"* ]]; then
      EXCLUDED_NAME_PATTERNS+=( -path "*/$pattern" -prune -o )
    else
      EXCLUDED_NAME_PATTERNS+=( -name "*${pattern}" -prune -o )
    fi
done

echo "======================================================"
echo -e "            TREE LIST OF FILES FOUND            "
echo "======================================================"

if [ ${#FILE_PATH_ARGS[@]} -gt 0 ]; then
    find "${FILE_PATH_ARGS[@]}" "${EXCLUDED_NAME_PATTERNS[@]}" -type d -print  | sort | while IFS= read -r dir; do
        echo -e "->" $(echo "$dir" | sed 's/[]^$*.[]/\\&/g' | sed 's/\// |---/g')
    done
fi

ARRAY_FILES=()
# mapfile used for reading inMAPput lines into an array; -d $'\0': Specifies that the delimiter is (\0).
# -print0 : ask find to separate file paths by '\0'
mapfile -d $'\0' ARRAY_FILES < <(find "${FILE_PATH_ARGS[@]}" "${EXCLUDED_NAME_PATTERNS[@]}" \
                                      -type f \( -name "*.hpp" -o -name "*.cpp" -o  -name "*.hpp.template" -o -name "*.cpp.template" \) \
                                      -print0 2>/dev/null)

################################
## III. MAIN LOOP
################################

for file in "${ARRAY_FILES[@]}"; do
  for container_name in "${!FORBIDDEN_CONTAINER_MAP[@]}"; do
    check_container_usage "$container_name"
  done
done

# Print section
for key in "${!FORBIDDEN_CONTAINER_MAP[@]}"; do
    if [[ "${FORBIDDEN_CONTAINER_MAP[$key]}" == "1" ]]; then
      echo $'\n'
      echo "----------------------------------------"
      echo "SUMMARY: Code rule violations found"
      echo "----------------------------------------"

      for container_name in "${!FORBIDDEN_CONTAINER_MAP[@]}"; do
        print_violation "$container_name"
      done
      
      echo ""
      exit 1;
    fi
done

echo $'\n'
echo "No code rule violations found"
exit 0
