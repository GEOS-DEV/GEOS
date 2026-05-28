#!/bin/bash
#SBATCH -t 00:30:00
#SBATCH -N 1
#SBATCH -A imcomp
#SBATCH -p pdebug

set -euo pipefail

# =============================================================================
# User settings
# =============================================================================
# List one or more PFW input files to stage and run.  Entries may be either full
# input filenames or case-name suffixes:
#
#   fileNames=( "pfw_input_myCase.py" )
#   fileNames=( "myCase" )                 # interpreted as pfw_input_myCase.py
#
fileNames=(
  "pfw_input_myCase.py"
)

# Directory containing the input files above.  Dependencies tagged in an input
# file with #[pfw_dependency] are resolved relative to this directory unless the
# dependency uses an absolute path or the pfw:<path> prefix.
fileLocation="<PATH_TO_INPUT_FILES>"

# Directory where per-case run directories will be created.  This should be on a
# large parallel file system such as Lustre or a workspace.
runLocation="<PATH_TO_RUN_LOCATION>"

# Directory containing particleFileWriter.py and the shared pfw_*.py modules.
# The default assumes this runClean script lives in the PFW directory.
particleFileWriterPath="${PFW_PATH:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)}"

# Python executable used to run particleFileWriter.py.
pythonCommand="${PFW_PYTHON:-${PFW_PYTHON_COMMAND:-python3}}"

# =============================================================================
# Normal users should not need to edit below this line.
# =============================================================================

normalize_input_file()
{
  local name="$1"

  if [[ "${name}" == *.py ]]; then
    echo "${name}"
  elif [[ "${name}" == pfw_input_* ]]; then
    echo "${name}.py"
  else
    echo "pfw_input_${name}.py"
  fi
}

case_name_from_input()
{
  local inputFile="$1"
  local stem="${inputFile%.py}"
  stem="${stem#pfw_input_}"
  echo "${stem}"
}

sourceDir="$(cd "${fileLocation}" && pwd -P)"
runRoot="$(mkdir -p "${runLocation}" && cd "${runLocation}" && pwd -P)"
userName="$(whoami)"
num_tasks="${1:-1}"

for requestedInput in "${fileNames[@]}"
do
  if [ -z "${requestedInput}" ] || [[ "${requestedInput}" == \<*\> ]]; then
    continue
  fi

  inputFile="$(normalize_input_file "${requestedInput}")"
  caseName="$(case_name_from_input "${inputFile}")"
  inputPath="${sourceDir}/${inputFile}"
  runDir="${runRoot}/${caseName}"

  if [ ! -f "${inputPath}" ]; then
    echo "ERROR: input file not found: ${inputPath}" >&2
    exit 1
  fi

  echo "Preparing ${caseName}"
  echo "  input file: ${inputPath}"
  echo "  run dir:    ${runDir}"
  echo "  tasks:      ${num_tasks}"

  aborted=false
  if [ -d "${runDir}" ] && [ -z "${SLURM_JOBID:-}" ]; then
    echo "Directory ${runDir} exists."
    while true; do
      read -r -p "Do you wish to overwrite? " yn
      case "${yn}" in
        [Yy]* ) echo "Overwriting..."; break ;;
        [Nn]* ) echo "Aborted overwrite..."; aborted=true; break ;;
        * ) echo "Please answer yes (Y/y) or no (N/n)." ;;
      esac
    done
  fi

  if [ "${aborted}" = true ]; then
    continue
  fi

  rm -rf "${runDir}"
  mkdir -p "${runDir}"

  # Copy the input file from the external source directory.  The marker file and
  # environment variable let particleFileWriter.py resolve #[pfw_dependency]
  # entries relative to the original input-file location after staging.
  cp "${inputPath}" "${runDir}/${inputFile}"
  printf '%s\n' "${sourceDir}" > "${runDir}/.pfw_input_source_dir"

  cp "${particleFileWriterPath}/particleFileWriter.py" "${runDir}/"
  for pfwFile in "${particleFileWriterPath}"/pfw_*.py; do
    [ -f "${pfwFile}" ] && cp "${pfwFile}" "${runDir}/"
  done
  for userDefsFile in "${particleFileWriterPath}"/userDefs_*.py; do
    [ -f "${userDefsFile}" ] && cp "${userDefsFile}" "${runDir}/"
  done

  cd "${runDir}"
  export PFW_INPUT_SOURCE_DIR="${sourceDir}"
  if [ "${num_tasks}" = "1" ]; then
    "${pythonCommand}" particleFileWriter.py "${inputFile}"
  else
    srun -n "${num_tasks}" "${pythonCommand}" particleFileWriter.py "${inputFile}"
  fi
  echo

done
