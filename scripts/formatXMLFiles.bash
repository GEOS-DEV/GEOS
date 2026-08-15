#!/usr/bin/env bash
# Format GEOS XML files with the Python format_xml tool (geos-xml-tools).
# Used by `make geosx_format_all_xml_files` and the XML formatting CI check.
# This is not schema validation.

set -uo pipefail

MODE=fix
METHOD=all
MAX_DIFFS=20

usage() {
    cat <<EOF
Usage: $0 [--check|--fix] [-a|--all|-g|--git] <FORMAT_SCRIPT> [<path>...]

  --check   Fail if any XML file is malformed or not in format_xml form
  --fix     Rewrite XML files in place (default)
  -a, --all Search all *.xml files under each path
  -g, --git Search git-tracked *.xml files under each path
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --check)
            MODE=check
            shift
            ;;
        --fix)
            MODE=fix
            shift
            ;;
        -a|--all)
            METHOD=all
            shift
            ;;
        -g|--git)
            METHOD=git
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        --)
            shift
            break
            ;;
        -*)
            echo "Error: unknown option: $1" >&2
            usage >&2
            exit 1
            ;;
        *)
            break
            ;;
    esac
done

if [[ -z "${1:-}" ]]; then
    usage >&2
    exit 1
fi

FORMAT_SCRIPT=$1
shift

if [[ ! -f "${FORMAT_SCRIPT}" && ! -x "${FORMAT_SCRIPT}" ]]; then
    if ! command -v "${FORMAT_SCRIPT}" >/dev/null 2>&1; then
        echo "Error: the format_xml script was not found: ${FORMAT_SCRIPT}" >&2
        exit 1
    fi
fi

if [[ "${METHOD}" = "git" ]] && ! command -v git >/dev/null 2>&1; then
    echo "Error: git is required when -g or --git is specified" >&2
    exit 1
fi

list_xml_files_all() {
    local search_path=$1
    find "${search_path}" -type f -name "*.xml" -not -path "*/.*" -print
}

list_xml_files_git() {
    local search_path=$1
    local git_root prefix
    git_root=$(cd "${search_path}" && git rev-parse --show-toplevel 2>/dev/null) || {
        echo "Error: ${search_path} does not appear to be part of a git repository" >&2
        return 1
    }
    prefix=$(cd "${search_path}" && git rev-parse --show-prefix 2>/dev/null)
    git --git-dir="${git_root}/.git" ls-files "${prefix}" | grep -e '.*[.]xml$' | sed "s|^|${git_root}/|g" || true
}

LOGFILE=xml_formatting_results.log
echo -n > "${LOGFILE}"

status=0
n_checked=0
n_invalid=0
n_unformatted=0
n_fixed=0
n_diffs_shown=0

process_file() {
    local file=$1

    if [[ "${MODE}" == "fix" ]]; then
        if ! "${FORMAT_SCRIPT}" "${file}" &>> "${LOGFILE}"; then
            echo "Invalid XML: ${file}"
            n_invalid=$((n_invalid + 1))
            status=1
        else
            n_fixed=$((n_fixed + 1))
        fi
        return
    fi

    local tmp
    tmp="$(mktemp)"
    cp "${file}" "${tmp}"
    if ! "${FORMAT_SCRIPT}" "${tmp}" &>> "${LOGFILE}"; then
        echo "Invalid XML: ${file}"
        n_invalid=$((n_invalid + 1))
        status=1
    elif ! cmp -s "${file}" "${tmp}"; then
        echo "Incorrect XML formatting: ${file}"
        if [[ ${n_diffs_shown} -lt ${MAX_DIFFS} ]]; then
            diff -u "${file}" "${tmp}" || true
            n_diffs_shown=$((n_diffs_shown + 1))
        fi
        n_unformatted=$((n_unformatted + 1))
        status=1
    fi
    rm -f "${tmp}"
}

if [[ $# -eq 0 ]]; then
    echo "Error: no search paths given" >&2
    exit 1
fi

for path in "$@"; do
    if [[ ! -d "${path}" ]]; then
        echo "Error: directory not found: ${path}" >&2
        exit 1
    fi
    while IFS= read -r file; do
        [[ -z "${file}" ]] && continue
        # examples/ contains git-tracked symlinks into inputFiles/. Skip them
        # so this check does not rewrite or require formatting that tree.
        [[ -L "${file}" ]] && continue
        n_checked=$((n_checked + 1))
        process_file "${file}"
    done < <(list_xml_files_${METHOD} "${path}")
done

if [[ "${MODE}" == "fix" ]]; then
    echo "XML formatting: ${n_checked} file(s) processed, ${n_invalid} invalid."
else
    echo "XML formatting: ${n_checked} file(s) checked, ${n_unformatted} incorrectly formatted, ${n_invalid} invalid."
    if [[ ${status} -ne 0 ]]; then
        echo "Run: make geosx_format_all_xml_files"
    fi
fi

exit "${status}"
