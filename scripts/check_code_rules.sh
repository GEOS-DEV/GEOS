#!/bin/bash


################################
## I. Function initializations
################################

# while
# ... done < <(grep -n "std::map\s*<" "$file") 
# Process the result of grep command as a file.
# This allows everything to be handled in the same shell environment.
check_container_usage() {
  local file="$1"

  if grep -qE "$FULL_REGEX" "$file"; then
    while IFS= read -r line; do
      for container in "${FORBIDDEN_EXPRESSIONS[@]}"; do
          local pattern="${container}[[:space:]]*[<\(]"
          if [[ "$line" =~ $pattern ]]; then
            local msg="    $file:$line"$'\n'
            ERRORS_CONTAINER[$container]+="$msg"
            ((FORBIDDEN_CONTAINER_MAP["$container"]++))
          fi
      done
    done < <(grep -nE "$FULL_REGEX" "$file")
  fi

}

print_violation()
{
    local container_name="$1"
    
    if [ "${FORBIDDEN_CONTAINER_MAP[$container_name]}" -ge 1 ];then
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

FILEPATH_PREFIX="src/coreComponents/"
FILEPATH_SUBFOLDERS=(
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

FILEPATH_EXCLUDE_PATTERNS=(
          "Datatype.hpp"     
          "StdContainerWrappers.hpp"     
          "MpiWrapper.cpp"
          "BufferOps_inline.hpp"
          "BufferOps.hpp"
          "PVTPackage"
          "hdf5_interface"
)

FILE_PATH_ARGS=()
for pattern in "${FILEPATH_SUBFOLDERS[@]}"; do
    if [ -d "${FILEPATH_PREFIX}${pattern}" ]; then
      FILE_PATH_ARGS+=("${FILEPATH_PREFIX}${pattern}")
    fi
done

EXCLUDED_NAME_PATTERNS=()
for pattern in "${FILEPATH_EXCLUDE_PATTERNS[@]}"; do
    if [[ ! "$pattern" == *".hpp"*  || "$pattern" == *".cpp"* ]]; then
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

REGEX_PATTERN=$(IFS='|'; echo "${FORBIDDEN_EXPRESSIONS[*]}")
FULL_REGEX="(${REGEX_PATTERN})"

for file in "${ARRAY_FILES[@]}"; do
    check_container_usage "$file"
done

# Print section
for key in "${!FORBIDDEN_CONTAINER_MAP[@]}"; do
    if [[ "${FORBIDDEN_CONTAINER_MAP[$key]}" -ge "1" ]]; then
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
