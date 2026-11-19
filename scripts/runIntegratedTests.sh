#!/usr/bin/env bash
#
# Run the GEOS integrated tests in the same Docker image as streak2 CI. The
# repository and tag are read from .devcontainer/Dockerfile and
# .devcontainer/devcontainer.json, respectively.
#
# Hypredrive is enabled for these integrated-test runs.
#
# Usage:
#   scripts/runIntegratedTests.sh
#   scripts/runIntegratedTests.sh --cpus 16 --memory 128g
#   scripts/runIntegratedTests.sh --cpus 32 --memory 256g          # match streak2
#   scripts/runIntegratedTests.sh --baselines /path/to/baselines
#   scripts/runIntegratedTests.sh --filter _iterative
#   scripts/runIntegratedTests.sh --generateBaselines
#   scripts/runIntegratedTests.sh --asan-ubsan --cpus 8 --no-pull
#
# If --baselines is omitted, the archive named in .integrated_tests.yaml is
# downloaded from GCS (public, ~1.5 GB) into .integrated-test-baselines/.
# A failed download does not abort the run: ATS continues without comparing
# against packed baselines.
#
# --generateBaselines builds origin/develop (worktree if this checkout is not
# develop) with 32 cores and 128g, copies every restart/hdf5 artifact the run
# produced into the baseline tree, then packs the archive named in develop's
# .integrated_tests.yaml into this repo's .integrated-test-baselines/. Use that
# when the GCS download is unavailable. ATS -a rebaseline is not used: it
# aborts the rest of an .ats file when one test has no restart output.
#
# After ATS, all GEOS logs that report Hypre solver iterations are harvested to
# .integrated-test-iterations/{hypredrive,legacy}.json. A `--force-legacy` run
# then compares that harvest to the hypredrive JSON for identical log/solve
# coverage and close aggregate iteration profiles. ATS result classifications
# are compared exactly, ignoring only run metadata, and any classification
# difference fails the run.
#
# Resource flags apply to the Docker container and to ATS rank limits.
# Remaining arguments are forwarded to geos_ats.sh.
#
# --asan-ubsan builds GEOS with AddressSanitizer and UndefinedBehaviorSanitizer
# and runs the complete integrated-test suite. TPLs are taken from the image
# when /spack-generated.cmake is present. If the image has no TPL host-config,
# pass --tpl-source-dir pointing at the geos-tpl superbuild; its build and
# install trees are created below the sanitizer root. The default sanitizer
# root and build directory are under /tmp, but both can be overridden.
#
# This script requires a native Linux host. The CI image is linux/amd64 and
# will not run on macOS (including Docker Desktop, Colima, and OrbStack).
#
# LLNL LC (Dane, etc.): `docker` is usually Podman. Image layers cannot live on
# Lustre/NFS/VAST (lsetxattr: operation not supported). Dane compute nodes are
# diskless: /tmp and /var/tmp are tmpfs (RAM), and $TMPDIR is often a Lustre
# symlink. This script then stores images under /tmp/geos-podman-$USER
# (override with GEOS_ITS_CONTAINER_STORE). The store is lost when the job
# ends and counts against node memory.
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GEOS_SRC_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

CPUS="${GEOS_ITS_CPUS:-8}"
MEMORY="${GEOS_ITS_MEMORY:-32g}"
NPROC=""
ATS_MAX_TESTS="${GEOS_ITS_ATS_MAX_TESTS:-}"
OMP_THREADS="${GEOS_ITS_OMP_THREADS:-}"
BASELINES=""
FILTER=""
CUTOFF="45m"
BUILD_DIR_NAME="build-integrated-tests"
PULL=1
GENERATE_BASELINES=0
FORCE_LEGACY=0
SANITIZERS="${GEOS_ITS_SANITIZERS:-0}"
TPL_HOST_CONFIG_PATH="${GEOS_ITS_TPL_HOST_CONFIG:-}"
TPL_SOURCE_PATH="${GEOS_ITS_TPL_SOURCE_DIR:-}"
CPUS_FROM_CLI=0
MEMORY_FROM_CLI=0
EXTRA_ATS_ARGS=()

usage()
{
  sed -n '2,/^$/p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
  cat <<EOF
Options:
  --cpus N              Docker CPU limit and maximum MPI ranks per test. Default ${CPUS}.
                        Under rootless Podman/Slurm the CPU cgroup limit is ignored.
  --memory SIZE         Docker memory limit, e.g. 32g, 128g, 256g (default ${MEMORY}).
                        Under rootless Podman/Slurm this cgroup limit is ignored.
                        A bare number is treated as gigabytes (380 -> 380g).
  --nproc N             Compile parallelism (default: same as --cpus)
  --ats-max-tests N     Maximum concurrent ATS test jobs (default: --nproc).
                        MPI rank capacity is unchanged.
  --baselines DIR       Host directory of ATS baseline tarballs (mounted read-only).
                        If omitted, download the archive from .integrated_tests.yaml.
                        If the download fails, continue without baseline comparison.
  --generateBaselines   Build develop at 32 cores / 128g, copy every restart
                        and hdf5 artifact the run produced, and pack the
                        tarball named in develop's .integrated_tests.yaml into
                        .integrated-test-baselines/. Cannot be used with --baselines.
  --filter EXPR         ATS name filter (tests whose name contains EXPR)
  --force-legacy        Set GEOS_HYPREDRV_FORCE_LEGACY=1 (legacy hypre path)
  --asan-ubsan          Build with ASan+UBSan and run the full test suite.
                        The default sanitizer root/build directory is below
                        /tmp; set GEOS_ITS_SANITIZER_ROOT or --build-dir to
                        use another location.
  --tpl-host-config FILE
                        Use this TPL host-config if the Docker image has none.
  --tpl-source-dir DIR  Build TPLs from this geos-tpl superbuild if the image
                        has no TPL host-config. Build/install paths use the
                        sanitizer root.
  --cutoff TIME         ATS cutoff (default ${CUTOFF})
  --build-dir NAME      GEOS build directory (under the repo normally; under
                        the sanitizer root by default for --asan-ubsan; an
                        absolute path may be used (default ${BUILD_DIR_NAME})
  --no-pull             Do not docker pull the image first
  -h, --help            Show this help
EOF
}

log()
{
  printf '[%s] %s\n' "$(date -u +"%Y-%m-%dT%H:%M:%SZ")" "$*"
}

die()
{
  printf 'ERROR: %s\n' "$*" >&2
  exit 1
}

SKIP_PACKED_BASELINES=0

# Print a boxed notice so a skipped baseline compare is hard to miss.
skip_baselines()
{
  local reason="$1"
  cat >&2 <<EOF

================================================================================
WARNING: Continuing WITHOUT packed-baseline comparison.
  ${reason}
  Tests will still run. geos_ats.sh will not download from GCS. Restart / VTK /
  HDF5 diffs will fail against an empty baseline tree (ATS requires those
  checks to be enabled). Pass --baselines DIR for a real packed archive.
================================================================================

EOF
  BASELINES=""
  SKIP_PACKED_BASELINES=1
}

read_tpl_tag()
{
  local json="${GEOS_SRC_DIR}/.devcontainer/devcontainer.json"
  [[ -f "${json}" ]] || die "missing ${json}"
  if command -v python3 >/dev/null 2>&1; then
    python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["build"]["args"]["GEOS_TPL_TAG"])' "${json}"
  else
    sed -n 's/.*"GEOS_TPL_TAG": *"\([^"]*\)".*/\1/p' "${json}" | head -n1
  fi
}

read_tpl_image_repository()
{
  local dockerfile repository
  if [[ -n "${GEOS_ITS_IMAGE_REPOSITORY:-}" ]]; then
    printf '%s\n' "${GEOS_ITS_IMAGE_REPOSITORY}"
    return 0
  fi

  dockerfile="${GEOS_SRC_DIR}/.devcontainer/Dockerfile"
  repository=""
  if [[ -f "${dockerfile}" ]]; then
    repository="$(sed -n \
      's/^[[:space:]]*FROM[[:space:]]\+\(docker\.io\/\)\?\([^:[:space:]]*\):.*/\2/p' \
      "${dockerfile}" | head -n1)"
  fi
  printf '%s\n' "${repository:-geosx/ubuntu24.04-gcc13}"
}

# Download the public GCS baseline tarball named in .integrated_tests.yaml.
# On any download/parse failure, leave BASELINES empty and keep going.
ensure_baselines()
{
  local yaml bucket baseline tarname dest url cache curl_status
  if [[ -n "${BASELINES}" ]]; then
    [[ -d "${BASELINES}" ]] || die "baselines directory not found: ${BASELINES}"
    return
  fi

  yaml="${GEOS_SRC_DIR}/.integrated_tests.yaml"
  if [[ ! -f "${yaml}" ]]; then
    skip_baselines "missing ${yaml}, cannot determine which archive to fetch."
    return
  fi
  bucket="$(grep -A 5 '^baselines:' "${yaml}" | awk '/bucket:/{print $2; exit}')"
  baseline="$(grep -A 5 '^baselines:' "${yaml}" | awk '/baseline:/{print $2; exit}')"
  if [[ -z "${bucket}" || -z "${baseline}" ]]; then
    skip_baselines "could not parse bucket/baseline from ${yaml}."
    return
  fi

  cache="${GEOS_SRC_DIR}/.integrated-test-baselines"
  tarname="$(basename "${baseline}").tar.gz"
  dest="${cache}/${tarname}"
  url="https://storage.googleapis.com/${bucket}/${baseline}.tar.gz"

  mkdir -p "${cache}"
  if [[ -f "${dest}" && -s "${dest}" ]]; then
    log "Cached baselines: ${dest}"
    BASELINES="${cache}"
    return
  fi

  if ! command -v curl >/dev/null 2>&1; then
    skip_baselines "curl is not on PATH, cannot download ${url}."
    return
  fi

  log "Download ${url}"
  set +e
  curl -fL --retry 3 --retry-delay 2 -o "${dest}.partial" "${url}"
  curl_status=$?
  set -e
  if [[ "${curl_status}" -ne 0 ]]; then
    rm -f "${dest}.partial"
    skip_baselines "failed to download ${url} (curl exit ${curl_status}). The archive may be private, missing from GCS, or the network blocked it."
    return
  fi
  if [[ ! -s "${dest}.partial" ]]; then
    rm -f "${dest}.partial"
    skip_baselines "downloaded ${url} but the file is empty."
    return
  fi
  mv "${dest}.partial" "${dest}"
  log "Saved ${dest}"
  BASELINES="${cache}"
}

# Use origin/develop (or local develop) so the packed archive matches what
# this repo's .integrated_tests.yaml expects. Copies this script into the
# worktree because develop does not contain it.
prepare_develop_for_baselines()
{
  local br wt ref sha
  br="$(git -C "${GEOS_SRC_DIR}" rev-parse --abbrev-ref HEAD 2>/dev/null || true)"
  if [[ "${br}" == "develop" ]]; then
    log "Using current develop checkout"
    return
  fi
  if git -C "${GEOS_SRC_DIR}" fetch origin develop; then
    ref=origin/develop
  elif git -C "${GEOS_SRC_DIR}" rev-parse --verify develop >/dev/null 2>&1; then
    ref=develop
    log "Could not fetch origin/develop; using local develop"
  else
    die "Need origin/develop (or a local develop branch) to generate baselines"
  fi
  sha="$(git -C "${GEOS_SRC_DIR}" rev-parse "${ref}")"
  wt="${ORIGINAL_SRC}-develop-baselines"
  if [[ -e "${wt}/.git" ]]; then
    git -C "${wt}" checkout --detach "${sha}"
  elif [[ -e "${wt}" ]]; then
    die "Refusing to use ${wt}; path exists and is not a git worktree"
  else
    git -C "${ORIGINAL_SRC}" worktree add --detach "${wt}" "${sha}"
  fi
  mkdir -p "${wt}/scripts"
  cp "${SCRIPT_DIR}/runIntegratedTests.sh" "${wt}/scripts/runIntegratedTests.sh"
  GEOS_SRC_DIR="${wt}"
  log "Using develop worktree ${wt} ($(git -C "${wt}" rev-parse --short HEAD))"
  git -C "${GEOS_SRC_DIR}" submodule update --init --recursive
}

resolve_path()
{
  local path="$1"
  if command -v realpath >/dev/null 2>&1; then
    realpath -m "${path}" 2>/dev/null || realpath "${path}" 2>/dev/null || printf '%s\n' "${path}"
  elif readlink -f "${path}" >/dev/null 2>&1; then
    readlink -f "${path}"
  else
    printf '%s\n' "${path}"
  fi
}

is_network_fs()
{
  local path fstype
  path="$(resolve_path "$1")"
  case "${path}" in
    /p/lustre*|/p/gpfs*|/p/vast*|/g/*|/usr/workspace*|/usr/WS*|/nfs/*|/collab/*) return 0 ;;
  esac
  while [[ -n "${path}" && "${path}" != "/" && ! -d "${path}" ]]; do
    path="$(dirname "${path}")"
  done
  [[ -d "${path}" ]] || path=/
  path="$(resolve_path "${path}")"
  fstype="$(df -T "${path}" 2>/dev/null | awk 'NR==2 { print $2 }')"
  case "${fstype}" in
    nfs|nfs4|nfs3|lustre|fuse.lustre|llite|gpfs|cifs|smb|smb3|vast|fuse.vast) return 0 ;;
  esac
  return 1
}

is_ram_fs()
{
  local path fstype
  path="$(resolve_path "$1")"
  fstype="$(df -T "${path}" 2>/dev/null | awk 'NR==2 { print $2 }')"
  case "${fstype}" in
    tmpfs|overlay|ramfs) return 0 ;;
  esac
  return 1
}

fs_supports_user_xattr()
{
  local dir="$1"
  local probe="${dir}/.geos-xattr-probe.$$"
  touch "${probe}" 2>/dev/null || return 1
  if python3 -c 'import os,sys; os.setxattr(sys.argv[1], "user.geos", b"1")' "${probe}" 2>/dev/null; then
    rm -f "${probe}"
    return 0
  fi
  rm -f "${probe}"
  return 1
}

# On LC, $TMPDIR is often /var/tmp/$USER -> symlink onto Lustre. Use the
# local mountpoint itself (e.g. /var/tmp/geos-podman-$USER), not $TMPDIR.
pick_local_podman_store()
{
  local d resolved store user
  user="${USER:-$(id -un)}"
  if [[ -n "${GEOS_ITS_CONTAINER_STORE:-}" ]]; then
    printf '%s\n' "$(resolve_path "${GEOS_ITS_CONTAINER_STORE}")"
    return 0
  fi
  for d in /var/tmp /tmp "${TMPDIR:-}"; do
    [[ -n "${d}" ]] || continue
    mkdir -p "${d}" 2>/dev/null || continue
    resolved="$(resolve_path "${d}")"
    if is_network_fs "${resolved}"; then
      continue
    fi
    store="${resolved}/geos-podman-${user}"
    mkdir -p "${store}" 2>/dev/null || continue
    store="$(resolve_path "${store}")"
    if is_network_fs "${store}"; then
      continue
    fi
    printf '%s\n' "${store}"
    return 0
  done
  return 1
}

# Rootless Podman cannot unpack overlay layers on Lustre/NFS. Point this
# invocation at node-local disk without rewriting ~/.config/containers.
configure_podman_storage()
{
  local graphroot store conf driver
  graphroot="$(${CONTAINER_CMD} info --format '{{.Store.GraphRoot}}' 2>/dev/null || true)"
  if [[ -n "${graphroot}" ]]; then
    graphroot="$(resolve_path "${graphroot}")"
  fi

  if [[ -n "${graphroot}" ]] && ! is_network_fs "${graphroot}"; then
    log "Podman image store: ${graphroot}"
    return
  fi

  store="$(pick_local_podman_store)" \
    || die "Podman needs a node-local image store (default is on a network filesystem${graphroot:+: ${graphroot}}). Set GEOS_ITS_CONTAINER_STORE to a real local directory (not a Lustre symlink), or run enable-podman."

  mkdir -p "${store}/storage" "${store}/run"
  store="$(resolve_path "${store}")"
  if is_network_fs "${store}"; then
    die "Resolved Podman store is still on a network filesystem (${store}). Set GEOS_ITS_CONTAINER_STORE to a node-local path."
  fi

  driver="overlay"
  if ! fs_supports_user_xattr "${store}"; then
    driver="vfs"
    log "Filesystem at ${store} does not support user xattrs; using the vfs storage driver."
  fi

  conf="${store}/storage.conf"
  cat > "${conf}" <<EOF
[storage]
driver = "${driver}"
runroot = "${store}/run"
graphroot = "${store}/storage"
EOF
  if [[ "${driver}" == overlay ]]; then
    cat >> "${conf}" <<EOF
[storage.options.overlay]
ignore_chown_errors = "true"
EOF
  else
    cat >> "${conf}" <<EOF
[storage.options.vfs]
ignore_chown_errors = "true"
EOF
  fi
  export CONTAINERS_STORAGE_CONF="${conf}"
  PODMAN_ROOT="${store}/storage"
  PODMAN_RUNROOT="${store}/run"
  if [[ -n "${graphroot}" ]]; then
    log "Podman default store is on a network filesystem (${graphroot})."
  else
    log "Podman default store could not be determined; assuming a network filesystem."
  fi
  log "Using local image store: ${PODMAN_ROOT}  (driver=${driver})"
  if is_ram_fs "${store}"; then
    log "That store is tmpfs/overlay (node RAM, not a disk). It disappears when the job ends and shares memory with the build and ATS."
  fi
}

normalize_memory()
{
  if [[ "${MEMORY}" =~ ^[0-9]+$ ]]; then
    log "Memory ${MEMORY} has no unit; treating it as ${MEMORY}g"
    MEMORY="${MEMORY}g"
  fi
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --cpus)       CPUS_FROM_CLI=1; CPUS="$2"; shift 2;;
    --memory)     MEMORY_FROM_CLI=1; MEMORY="$2"; shift 2;;
    --nproc)      NPROC="$2"; shift 2;;
    --ats-max-tests) ATS_MAX_TESTS="$2"; shift 2;;
    --baselines)  BASELINES="$2"; shift 2;;
    --generateBaselines) GENERATE_BASELINES=1; shift;;
    --filter)     FILTER="$2"; shift 2;;
    --force-legacy) FORCE_LEGACY=1; shift;;
    --asan-ubsan|--build-asan-ubsan|--build-sanitizers) SANITIZERS=1; shift;;
    --tpl-host-config) TPL_HOST_CONFIG_PATH="$2"; shift 2;;
    --tpl-source-dir) TPL_SOURCE_PATH="$2"; shift 2;;
    --cutoff)     CUTOFF="$2"; shift 2;;
    --build-dir)  BUILD_DIR_NAME="$2"; shift 2;;
    --no-pull)    PULL=0; shift;;
    -h|--help)    usage; exit 0;;
    --)           shift; EXTRA_ATS_ARGS+=("$@"); break;;
    *)            EXTRA_ATS_ARGS+=("$1"); shift;;
  esac
done

if [[ "${GENERATE_BASELINES}" -eq 1 ]]; then
  [[ -z "${BASELINES}" ]] || die "--generateBaselines cannot be used with --baselines"
  [[ "${CPUS_FROM_CLI}" -eq 1 ]] || CPUS=32
  [[ "${MEMORY_FROM_CLI}" -eq 1 ]] || MEMORY=128g
fi
[[ -n "${NPROC}" ]] || NPROC="${CPUS}"
normalize_memory

if [[ -z "${ATS_MAX_TESTS}" ]]; then
  ATS_MAX_TESTS="${NPROC}"
fi

if ! [[ "${ATS_MAX_TESTS}" =~ ^[1-9][0-9]*$ ]]; then
  die "ATS_MAX_TESTS must be a positive integer, got '${ATS_MAX_TESTS}'."
fi

if [[ "${SANITIZERS}" -eq 1 ]]; then
  [[ -z "${FILTER}" ]] || die "--asan-ubsan always runs the full integrated-test suite; omit --filter"
  [[ ${#EXTRA_ATS_ARGS[@]} -eq 0 ]] || die "--asan-ubsan does not accept forwarded ATS filters; it runs the full suite"
  [[ "${GENERATE_BASELINES}" -eq 0 ]] || die "--asan-ubsan cannot be combined with --generateBaselines"

  SANITIZER_ROOT="${GEOS_ITS_SANITIZER_ROOT:-/tmp/geos-integrated-tests-asan-ubsan-${USER:-$(id -u)}}"
  if [[ "${SANITIZER_ROOT}" != /* ]]; then
    SANITIZER_ROOT="${GEOS_SRC_DIR}/${SANITIZER_ROOT}"
  fi
  SANITIZER_ROOT="$(resolve_path "${SANITIZER_ROOT}")"
  if [[ "${BUILD_DIR_NAME}" == "build-integrated-tests" ]]; then
    BUILD_DIR_NAME="${SANITIZER_ROOT}/geos-build"
  elif [[ "${BUILD_DIR_NAME}" != /* ]]; then
    BUILD_DIR_NAME="${GEOS_SRC_DIR}/${BUILD_DIR_NAME}"
  fi
  BUILD_DIR_NAME="$(resolve_path "${BUILD_DIR_NAME}")"
fi

# ---------------------------------------------------------------------------
# Host: launch the streak2 CI image
# ---------------------------------------------------------------------------
if [[ "${GEOS_IN_CONTAINER:-}" != 1 ]]; then
  if [[ "$(uname -s)" == Darwin ]]; then
    die "This script requires a native Linux host. The configured TPL image is linux/amd64 and is not supported on macOS."
  fi

  CONTAINER_CMD=""
  PODMAN_ROOT=""
  PODMAN_RUNROOT=""
  if command -v podman >/dev/null 2>&1; then
    CONTAINER_CMD=podman
  elif command -v docker >/dev/null 2>&1; then
    CONTAINER_CMD=docker
  else
    die "docker or podman is not on PATH"
  fi

  ORIGINAL_SRC="${GEOS_SRC_DIR}"
  BASELINE_CACHE="${ORIGINAL_SRC}/.integrated-test-baselines"
  BASELINE_ARCHIVE=""
  if [[ "${GENERATE_BASELINES}" -eq 1 ]]; then
    SKIP_PACKED_BASELINES=1
    prepare_develop_for_baselines
    yaml="${GEOS_SRC_DIR}/.integrated_tests.yaml"
    [[ -f "${yaml}" ]] || die "missing ${yaml}"
    BASELINE_ARCHIVE="$(grep -A 5 '^baselines:' "${yaml}" | awk '/baseline:/{print $2; exit}')"
    BASELINE_ARCHIVE="$(basename "${BASELINE_ARCHIVE}")"
    [[ -n "${BASELINE_ARCHIVE}" ]] || die "could not parse baseline name from ${yaml}"
    mkdir -p "${BASELINE_CACHE}"
    log "Will pack ${BASELINE_CACHE}/${BASELINE_ARCHIVE}.tar.gz"
  else
    ensure_baselines
    if [[ -n "${BASELINES}" ]]; then
      log "Baselines: ${BASELINES}"
    else
      log "Baselines: none (ATS will not compare against a packed archive)"
    fi
  fi

  TPL_TAG="$(read_tpl_tag)"
  IMAGE_REPOSITORY="$(read_tpl_image_repository)"
  IMAGE="${IMAGE_REPOSITORY}:${TPL_TAG}"

  if [[ "${SANITIZERS}" -eq 1 ]]; then
    [[ -z "${TPL_HOST_CONFIG_PATH}" || -f "${TPL_HOST_CONFIG_PATH}" ]] \
      || die "TPL host-config not found: ${TPL_HOST_CONFIG_PATH}"
    [[ -z "${TPL_SOURCE_PATH}" || -d "${TPL_SOURCE_PATH}" ]] \
      || die "TPL source directory not found: ${TPL_SOURCE_PATH}"
    mkdir -p "${SANITIZER_ROOT}" "${BUILD_DIR_NAME}"
  fi

  log "Image:   ${IMAGE}  (streak2 integrated_tests)"
  log "Runtime: ${CONTAINER_CMD}"
  log "CPUs:    ${CPUS}"
  log "Memory:  ${MEMORY}"
  log "Nproc:   ${NPROC}"
  log "Cutoff:  ${CUTOFF}"
  [[ -n "${FILTER}" ]] && log "Filter:   ${FILTER}"

  if [[ "${CONTAINER_CMD}" == podman ]]; then
    configure_podman_storage
  fi

  runtime_args=()
  if [[ -n "${PODMAN_ROOT}" ]]; then
    runtime_args+=(--root "${PODMAN_ROOT}" --runroot "${PODMAN_RUNROOT}")
  fi
  if [[ "${CONTAINER_CMD}" == podman ]]; then
    # systemd + Slurm: rootless Podman cannot mkdir under slurmstepd.scope.
    runtime_args+=(--cgroup-manager=cgroupfs)
  fi

  if [[ "${PULL}" -eq 1 ]]; then
    log "Pull ${IMAGE}"
    "${CONTAINER_CMD}" "${runtime_args[@]}" pull --platform linux/amd64 "${IMAGE}"
  fi

  docker_args=(
    --rm
    --platform linux/amd64
    --workdir /workspace
    --entrypoint /bin/bash
    --env GEOS_IN_CONTAINER=1
    --env GEOS_ITS_CPUS="${CPUS}"
    --env GEOS_ITS_NPROC="${NPROC}"
    --env GEOS_ITS_ATS_MAX_TESTS="${ATS_MAX_TESTS}"
    --env GEOS_ITS_OMP_THREADS="${OMP_THREADS}"
    --env GEOS_ITS_CUTOFF="${CUTOFF}"
    --env GEOS_ITS_FILTER="${FILTER}"
    --env GEOS_ITS_BUILD_DIR="${BUILD_DIR_NAME}"
    --env GEOS_ITS_SKIP_BASELINES="${SKIP_PACKED_BASELINES}"
    --env GEOS_ITS_GENERATE_BASELINES="${GENERATE_BASELINES}"
    --env GEOS_ITS_BASELINE_ARCHIVE="${BASELINE_ARCHIVE}"
    --env GEOS_ITS_HOST_SRC="${GEOS_SRC_DIR}"
    --env GEOS_ITS_FORCE_LEGACY="${FORCE_LEGACY}"
    --env GEOS_ITS_SANITIZERS="${SANITIZERS}"
    --mount "type=bind,src=${GEOS_SRC_DIR},dst=/workspace"
  )
  if [[ "${SANITIZERS}" -eq 0 && "${BUILD_DIR_NAME}" == /* ]]; then
    mkdir -p "${BUILD_DIR_NAME}"
    docker_args+=(--mount "type=bind,src=${BUILD_DIR_NAME},dst=${BUILD_DIR_NAME}")
  fi
  if [[ "${FORCE_LEGACY}" -eq 1 ]]; then
    docker_args+=(--env GEOS_HYPREDRV_FORCE_LEGACY=1)
    log "Forcing legacy hypre path (GEOS_HYPREDRV_FORCE_LEGACY=1)"
  fi

  if [[ "${SANITIZERS}" -eq 1 ]]; then
    docker_args+=(--env GEOS_ITS_SANITIZER_ROOT="${SANITIZER_ROOT}")
    docker_args+=(--cap-add=SYS_PTRACE --security-opt=seccomp=unconfined)
    docker_args+=(--mount "type=bind,src=${SANITIZER_ROOT},dst=${SANITIZER_ROOT}")
    case "${BUILD_DIR_NAME}" in
      "${SANITIZER_ROOT}"|"${SANITIZER_ROOT}"/*) ;;
      *) docker_args+=(--mount "type=bind,src=${BUILD_DIR_NAME},dst=${BUILD_DIR_NAME}") ;;
    esac
    if [[ -n "${TPL_HOST_CONFIG_PATH}" ]]; then
      TPL_HOST_CONFIG_PATH="$(resolve_path "${TPL_HOST_CONFIG_PATH}")"
      docker_args+=(
        --env GEOS_ITS_TPL_HOST_CONFIG=/tmp/geos-its-tpl-host-config.cmake
        --mount "type=bind,src=${TPL_HOST_CONFIG_PATH},dst=/tmp/geos-its-tpl-host-config.cmake,readonly"
      )
    fi
    if [[ -n "${TPL_SOURCE_PATH}" ]]; then
      TPL_SOURCE_PATH="$(resolve_path "${TPL_SOURCE_PATH}")"
      docker_args+=(
        --env GEOS_ITS_TPL_SOURCE_DIR=/tmp/geos-its-tpl-source
        --mount "type=bind,src=${TPL_SOURCE_PATH},dst=/tmp/geos-its-tpl-source,readonly"
      )
    fi
  fi

  if [[ "${CONTAINER_CMD}" == podman ]]; then
    # Rootless Podman under Slurm cannot create cgroup slices. Do not pass
    # --cpus/--memory. --cgroups=disabled requires crun; runc rejects it
    # ("not compatible with NoCgroups"). ATS still honors NPROC.
    docker_args+=(--security-opt=seccomp=unconfined --systemd=false)
    crun_bin=""
    if command -v crun >/dev/null 2>&1; then
      crun_bin="$(command -v crun)"
    elif [[ -x /usr/bin/crun ]]; then
      crun_bin=/usr/bin/crun
    fi
    if [[ -n "${crun_bin}" ]]; then
      docker_args+=(--runtime "${crun_bin}" --cgroups=disabled)
      log "Podman: crun (${crun_bin}), cgroups disabled; ATS nproc=${NPROC}"
    else
      log "Podman: runc + cgroupfs, no CPU/memory limits (crun not found); ATS nproc=${NPROC}"
    fi
  else
    docker_args+=(--cpus="${CPUS}" --memory="${MEMORY}")
  fi

  if [[ "${GENERATE_BASELINES}" -eq 1 ]]; then
    docker_args+=(--mount "type=bind,src=${BASELINE_CACHE},dst=/tmp/geos/generated-baselines")
    if [[ "${GEOS_SRC_DIR}" != "${ORIGINAL_SRC}" ]]; then
      # Worktree .git is a pointer into the main repo; ATS needs that for history.
      docker_args+=(--mount "type=bind,src=${ORIGINAL_SRC},dst=${ORIGINAL_SRC},readonly")
    fi
  elif [[ -n "${BASELINES}" ]]; then
    docker_args+=(--mount "type=bind,src=$(cd "${BASELINES}" && pwd),dst=/tmp/geos/baselines,readonly")
  fi

  if [[ ${#EXTRA_ATS_ARGS[@]} -gt 0 ]]; then
    exec "${CONTAINER_CMD}" "${runtime_args[@]}" run "${docker_args[@]}" "${IMAGE}" \
      /workspace/scripts/runIntegratedTests.sh "${EXTRA_ATS_ARGS[@]}"
  else
    exec "${CONTAINER_CMD}" "${runtime_args[@]}" run "${docker_args[@]}" "${IMAGE}" \
      /workspace/scripts/runIntegratedTests.sh
  fi
fi

# ---------------------------------------------------------------------------
# Container: configure, build, run ATS (HypreDrive enabled)
# ---------------------------------------------------------------------------
NPROC="${GEOS_ITS_NPROC:-8}"
ATS_MAX_TESTS="${GEOS_ITS_ATS_MAX_TESTS:-${ATS_MAX_TESTS}}"
CUTOFF="${GEOS_ITS_CUTOFF:-45m}"
FILTER="${GEOS_ITS_FILTER:-}"
SANITIZERS="${GEOS_ITS_SANITIZERS:-${SANITIZERS}}"
SANITIZER_ROOT="${GEOS_ITS_SANITIZER_ROOT:-/tmp/geos-integrated-tests-asan-ubsan-${USER:-$(id -u)}}"
TPL_SOURCE_PATH="${GEOS_ITS_TPL_SOURCE_DIR:-${TPL_SOURCE_PATH}}"
OMP_THREADS="${GEOS_ITS_OMP_THREADS:-${OMP_THREADS}}"
if [[ "${GEOS_ITS_BUILD_DIR:-build-integrated-tests}" == /* ]]; then
  BUILD_DIR="${GEOS_ITS_BUILD_DIR}"
else
  BUILD_DIR="/workspace/${GEOS_ITS_BUILD_DIR:-build-integrated-tests}"
fi
VENV_DIR="${BUILD_DIR}/ats-venv"
VENV_PY="${VENV_DIR}/bin/python3"
ATS_WORKING_DIR="${BUILD_DIR}/ats-working"
ATS_BASELINE_DIR="${BUILD_DIR}/ats-baselines"
HOST_CONFIG="${GEOS_ITS_TPL_HOST_CONFIG:-/spack-generated.cmake}"

if [[ "${SANITIZERS}" -eq 1 ]]; then
  if [[ "${SANITIZER_ROOT}" != /* ]]; then
    SANITIZER_ROOT="/workspace/${SANITIZER_ROOT}"
  fi
fi

if ! [[ "${NPROC}" =~ ^[1-9][0-9]*$ ]]; then
  die "NPROC must be a positive integer, got '${NPROC}'."
fi
if ! [[ "${ATS_MAX_TESTS}" =~ ^[1-9][0-9]*$ ]]; then
  die "ATS_MAX_TESTS must be a positive integer, got '${ATS_MAX_TESTS}'."
fi

# Open MPI refuses root launches in CI containers unless both guard variables
# are set. Keep this separate from openmpi_args because repeated geos-ats
# overrides replace values.
export OMPI_ALLOW_RUN_AS_ROOT=1
export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
ATS_ARGUMENTS="--machine openmpi --ats openmpi_mpirun=/usr/bin/mpirun --ats openmpi_procspernode=${NPROC} --ats openmpi_maxprocs=${NPROC} --ats cutoff=${CUTOFF}"
log "ATS concurrent test jobs: ${ATS_MAX_TESTS} (MPI rank capacity: ${NPROC})"

# ATS always calls collect_baselines() on run. If the named archive is already
# marked present (.blob_name), it skips GCS. Seed that marker so a 403 does not
# hang in google.cloud.storage retries.
seed_dummy_ats_baselines()
{
  local yaml blob
  yaml="/workspace/.integrated_tests.yaml"
  blob=""
  if [[ -f "${yaml}" ]]; then
    blob="$(grep -A 5 '^baselines:' "${yaml}" | awk '/baseline:/{print $2; exit}')"
  fi
  blob="$(basename "${blob:-skip-packed-baselines}")"
  mkdir -p "${ATS_BASELINE_DIR}"
  # ATS compares this file with == and does not strip whitespace.
  printf '%s' "${blob}" > "${ATS_BASELINE_DIR}/.blob_name"
  log "ATS will not fetch GCS (seeded ${ATS_BASELINE_DIR}/.blob_name=${blob})"
}

# Mirror restart roots (and their rank hdf5 dirs) plus standalone hdf5
# (curve-check histories) from the ATS working tree into the baseline tree.
copy_run_outputs_to_baselines()
{
  local src="$1" dst="$2" n=0 nh=0
  local root rel dest datadir hdf
  [[ -d "${src}" ]] || die "working directory not found: ${src}"
  mkdir -p "${dst}"
  while IFS= read -r -d '' root; do
    rel="${root#"${src}"/}"
    dest="${dst}/${rel}"
    mkdir -p "$(dirname "${dest}")"
    cp -a "${root}" "${dest}"
    datadir="${root%.root}"
    if [[ -d "${datadir}" ]]; then
      rm -rf "${dest%.root}"
      cp -a "${datadir}" "${dest%.root}"
    fi
    n=$(( n + 1 ))
  done < <(find "${src}" -type f -name '*.root' -print0)
  while IFS= read -r -d '' hdf; do
    rel="${hdf#"${src}"/}"
    dest="${dst}/${rel}"
    mkdir -p "$(dirname "${dest}")"
    cp -a "${hdf}" "${dest}"
    nh=$(( nh + 1 ))
  done < <(find "${src}" -type f -name '*.hdf5' ! -path '*_restart_*/*' -print0)
  log "Copied ${n} restart roots and ${nh} standalone hdf5 files into ${dst}"
}

or_die()
{
  "$@"
  local status=$?
  if [[ ${status} -ne 0 ]]; then
    echo "ERROR ${status} command: $*" >&2
    exit ${status}
  fi
}

SANITIZER_FLAGS="-fsanitize=address,undefined -fno-omit-frame-pointer"
SANITIZER_LINK_FLAGS="-fsanitize=address,undefined"
SANITIZER_CMAKE_ARGS=()

# Build a complete sanitizer-compatible TPL stack when the runtime image does
# not provide its generated host-config. The sanitizer path uses the geos-tpl
# Spack driver because it needs Spack to generate a GEOS host-config carrying
# sanitizer flags. The legacy CMake superbuild also provides Hypredrive for
# non-Spack TPL builds. The source checkout is mounted read-only; every
# generated Spack, build, and install path is below the configurable
# SANITIZER_ROOT (which defaults to /tmp for this script).
build_sanitized_tpls()
{
  local tpl_root="${SANITIZER_ROOT}/tpls"
  local tpl_build="${tpl_root}/build"
  local tpl_install="${tpl_root}/install"
  local spack_env_file="${tpl_root}/spack-sanitizers.yaml"
  local gcc_version target os_name mpi_version
  local candidate
  local -a tpl_spec

  [[ -n "${TPL_SOURCE_PATH}" ]] \
    || die "${HOST_CONFIG} is not available; pass --tpl-source-dir PATH to a geos-tpl checkout"
  [[ -x "${TPL_SOURCE_PATH}/build.sh" ]] \
    || die "TPL source must contain an executable build.sh: ${TPL_SOURCE_PATH}"

  mkdir -p "${tpl_root}" "${tpl_build}" "${tpl_install}"
  command -v gcc >/dev/null 2>&1 || die "gcc is required to build sanitizer TPLs"
  gcc_version="$(gcc -dumpfullversion -dumpversion 2>/dev/null || true)"
  [[ -n "${gcc_version}" ]] || gcc_version="$(gcc -dumpversion)"
  target="$(gcc -dumpmachine | sed 's/-.*//')"
  os_name="linux"
  if [[ -r /etc/os-release ]]; then
    # shellcheck disable=SC1091
    . /etc/os-release
    os_name="${ID:-linux}${VERSION_ID:-}"
  fi
  [[ -x /usr/bin/mpirun ]] || die "Open MPI (/usr/bin/mpirun) is required to build sanitizer TPLs"
  mpi_version="$(/usr/bin/mpirun --version 2>/dev/null | awk '/Open MPI/{print $NF; exit}')"
  [[ -n "${mpi_version}" ]] || die "could not determine the Open MPI version in the image"

  # Supplying an environment file prevents uberenv from copying a generated
  # spack.yaml into the read-only source checkout and records sanitizer flags
  # in Spack's compiler configuration. Those flags are then present in the
  # generated GEOS host-config as well as in the TPL builds.
  cat > "${spack_env_file}" <<EOF
spack:
  view: false
  repos:
    geos_sanitizer:
      destination: /workspace/scripts/sanitizer_spack_repo
  compilers:
  - compiler:
      spec: gcc@${gcc_version}
      paths:
        cc: /usr/bin/gcc
        cxx: /usr/bin/g++
        f77: /usr/bin/gfortran
        fc: /usr/bin/gfortran
      flags:
        cflags: ${SANITIZER_FLAGS}
        cxxflags: ${SANITIZER_FLAGS}
        fflags: ${SANITIZER_FLAGS}
        ldflags: ${SANITIZER_LINK_FLAGS}
      operating_system: ${os_name}
      target: ${target}
      modules: []
      environment: {}
      extra_rpaths: []
  packages:
    all:
      target: [${target}]
    mpi:
      require:
      - openmpi
    openmpi:
      buildable: false
      externals:
      - spec: openmpi@${mpi_version}
        prefix: /usr
EOF

  tpl_spec=(
    "~docs +hypre +hypredrive +vtk +caliper +shared +openmp ~trilinos ~petsc ~pygeosx %gcc ^hypre@3.2.0 ^vtk generator=ninja"
  )
  log "TPL image host-config missing; building sanitizer TPLs under ${tpl_root}"
  export CFLAGS="${SANITIZER_FLAGS}${CFLAGS:+ ${CFLAGS}}"
  export CXXFLAGS="${SANITIZER_FLAGS}${CXXFLAGS:+ ${CXXFLAGS}}"
  export FFLAGS="${SANITIZER_FLAGS}${FFLAGS:+ ${FFLAGS}}"
  export LDFLAGS="${SANITIZER_LINK_FLAGS}${LDFLAGS:+ ${LDFLAGS}}"
  or_die "${TPL_SOURCE_PATH}/build.sh" z5tux \
    --no-modules --spack --build-type RelWithDebInfo \
    --build-dir "${tpl_build}" --install-dir "${tpl_install}" \
    --jobs "${NPROC}" --top-jobs 1 --spec "${tpl_spec[0]}" \
    --log "${tpl_root}/build.log" -- \
    --spack-env-file "${spack_env_file}"

  HOST_CONFIG=""
  while IFS= read -r -d '' candidate; do
    if grep -q 'HYPREDRV_DIR' "${candidate}" && grep -q 'ENABLE_HYPREDRV' "${candidate}"; then
      HOST_CONFIG="${candidate}"
      break
    fi
  done < <(find -L "${tpl_root}" -type f -name '*.cmake' -print0)
  [[ -n "${HOST_CONFIG}" ]] \
    || die "sanitizer TPL build completed without a GEOS host-config containing HYPREDRV_DIR"
  log "Sanitizer TPL host-config: ${HOST_CONFIG}"
}

if [[ "${SANITIZERS}" -eq 1 ]]; then
  if [[ ! -f "${HOST_CONFIG}" ]]; then
    build_sanitized_tpls
  fi
  [[ -f "${HOST_CONFIG}" ]] \
    || die "TPL host-config not found: ${HOST_CONFIG}; pass --tpl-host-config or --tpl-source-dir"
  SANITIZER_CMAKE_ARGS=(
    "-DCMAKE_C_FLAGS:STRING=${SANITIZER_FLAGS}"
    "-DCMAKE_CXX_FLAGS:STRING=${SANITIZER_FLAGS}"
    "-DCMAKE_EXE_LINKER_FLAGS:STRING=${SANITIZER_LINK_FLAGS}"
    "-DCMAKE_SHARED_LINKER_FLAGS:STRING=${SANITIZER_LINK_FLAGS}"
    "-DENABLE_WARNINGS_AS_ERRORS:BOOL=OFF"
  )

  append_sanitizer_options()
  {
    local existing="$1" required="$2"
    if [[ -n "${existing}" ]]; then
      printf '%s:%s' "${existing}" "${required}"
    else
      printf '%s' "${required}"
    fi
  }

  # Append the required settings after caller-provided options so leak
  # detection and the report locations cannot be disabled accidentally.
  export ASAN_OPTIONS="$(append_sanitizer_options "${ASAN_OPTIONS:-}" "abort_on_error=1:detect_leaks=1:leak_check_at_exit=1:fast_unwind_on_malloc=0:print_summary=1:halt_on_error=1:log_path=${SANITIZER_ROOT}/asan")"
  # LeakSanitizer is enabled explicitly as well as through ASAN_OPTIONS. Keep
  # its reports separate so leak findings remain visible in the final report.
  export LSAN_OPTIONS="$(append_sanitizer_options "${LSAN_OPTIONS:-}" "detect_leaks=1:leak_check_at_exit=1:fast_unwind_on_malloc=0:print_suppressions=0:suppressions=/workspace/scripts/lsan.supp:log_path=${SANITIZER_ROOT}/lsan")"
  # UBSan recovers so all tests in the full suite can run. The report scanner
  # below still turns any UBSan or LeakSanitizer diagnostic into a failing
  # sanitizer run.
  export UBSAN_OPTIONS="$(append_sanitizer_options "${UBSAN_OPTIONS:-}" "halt_on_error=0:print_stacktrace=1:log_path=${SANITIZER_ROOT}/ubsan")"
  OMP_THREADS="${OMP_THREADS:-1}"
  if ! [[ "${OMP_THREADS}" =~ ^[1-9][0-9]*$ ]]; then
    die "GEOS_ITS_OMP_THREADS must be a positive integer, got '${OMP_THREADS}'."
  fi
  export OMP_NUM_THREADS="${OMP_THREADS}"
  export OMP_DYNAMIC=FALSE
  log "Sanitizers: ASAN_OPTIONS=${ASAN_OPTIONS}"
  log "Sanitizers: LSAN_OPTIONS=${LSAN_OPTIONS}"
  log "Sanitizers: UBSAN_OPTIONS=${UBSAN_OPTIONS}"
  log "Sanitizers: OMP_NUM_THREADS=${OMP_NUM_THREADS}"
fi

log "Create ATS virtualenv"
if [[ ! -x "${VENV_PY}" ]]; then
  or_die python3 -m venv --system-site-packages "${VENV_DIR}"
fi
git config --global --add safe.directory /workspace 2>/dev/null || true
"${VENV_PY}" -m pip install --upgrade pip setuptools wheel || true
"${VENV_PY}" -m pip install --disable-pip-version-check --upgrade matplotlib || true

log "Configure (ENABLE_HYPREDRV=ON)"
BUILD_TYPE=Release
if [[ "${SANITIZERS}" -eq 1 ]]; then
  BUILD_TYPE=RelWithDebInfo
fi
or_die cmake -S "${GEOS_SRC_DIR}/src" -B "${BUILD_DIR}" -G Ninja \
  -C "${HOST_CONFIG}" \
  "-DCMAKE_BUILD_TYPE=${BUILD_TYPE}" \
  -DGEOS_INSTALL_SCHEMA=ON \
  -DENABLE_HYPRE=ON \
  -DENABLE_HYPREDRV=ON \
  -DENABLE_HYPRE_DEVICE=CPU \
  -DENABLE_TRILINOS=OFF \
  -DGEOS_LA_INTERFACE:PATH=Hypre \
  -DGEOS_ENABLE_BOUNDS_CHECK=ON \
  "-DBLT_MPI_COMMAND_APPEND=--allow-run-as-root;--oversubscribe" \
  "-DATS_ARGUMENTS:STRING=${ATS_ARGUMENTS}" \
  -DPython3_EXECUTABLE="${VENV_PY}" \
  -DATS_BASELINE_DIR="${ATS_BASELINE_DIR}" \
  -DATS_WORKING_DIR="${ATS_WORKING_DIR}" \
  "${SANITIZER_CMAKE_ARGS[@]}"

ENABLE_HYPREDRV_VAL="$(sed -n 's/^ENABLE_HYPREDRV:[^=]*=//p' "${BUILD_DIR}/CMakeCache.txt" | head -n1)"
log "CMake ENABLE_HYPREDRV=${ENABLE_HYPREDRV_VAL}"
if [[ "${ENABLE_HYPREDRV_VAL}" != "ON" ]]; then
  die "ENABLE_HYPREDRV is ${ENABLE_HYPREDRV_VAL}; expected ON"
fi

export ATS_FILTER="np<=${NPROC}"
log "ATS_FILTER=${ATS_FILTER}"

log "Build geosx"
or_die cmake --build "${BUILD_DIR}" --target geosx --parallel "${NPROC}"
log "Build ats_environment"
or_die cmake --build "${BUILD_DIR}" --target ats_environment --parallel "${NPROC}"

if [[ "${SANITIZERS}" -eq 1 ]]; then
  # The sanitizer report must describe this invocation only. These paths are
  # generated below the configurable sanitizer root and are not source or
  # baseline data.
  log "Clear previous sanitizer ATS outputs"
  rm -rf -- "${ATS_WORKING_DIR}" "${BUILD_DIR}/integratedTests/TestResults"
  rm -f -- "${SANITIZER_ROOT}"/asan.* "${SANITIZER_ROOT}"/lsan.* "${SANITIZER_ROOT}"/ubsan.* \
    "${SANITIZER_ROOT}/sanitizer-reports.txt" "${SANITIZER_ROOT}/sanitizer-reports.txt.tmp"
fi

cd "${BUILD_DIR}"
ATS_CMD=(integratedTests/geos_ats.sh)
if [[ "${GEOS_ITS_GENERATE_BASELINES:-0}" == 1 ]]; then
  seed_dummy_ats_baselines
  log "Generate mode: run ATS without packed baselines, then copy outputs and pack"
elif [[ "${GEOS_ITS_SKIP_BASELINES:-0}" == 1 || ! -d /tmp/geos/baselines ]]; then
  seed_dummy_ats_baselines
  cat >&2 <<EOF

================================================================================
WARNING: No packed baselines are mounted.
  geos_ats.sh will not download from GCS. Tests still run with restart / VTK /
  HDF5 checks enabled (ATS rejects -c none). Those checks will fail until you
  pass --baselines DIR on the host.
================================================================================

EOF
else
  ATS_CMD+=(--baselineCacheDirectory /tmp/geos/baselines --delete-old-baselines)
  log "Using mounted baselines at /tmp/geos/baselines"
fi
if [[ -n "${FILTER}" ]]; then
  # geos_ats uses -f for "allow failed tests".  Pass the requested test-name
  # expression through to ATS' own --filter option instead.
  ATS_CMD+=(--ats "filter=\"${FILTER}\" in SELF.name")
fi
if [[ ${#EXTRA_ATS_ARGS[@]} -gt 0 ]]; then
  ATS_CMD+=("${EXTRA_ATS_ARGS[@]}")
fi

# This route configures geos-ats' Open MPI machine explicitly. Capture the
# selected launcher's identity once, and fail before tests if an image change
# would make the Open MPI placement environment below ineffective.
log "Diagnostic: mpirun --version"
MPIRUN_VERSION_STATUS=0
MPIRUN_VERSION_OUTPUT="$(/usr/bin/mpirun --version 2>&1)" || MPIRUN_VERSION_STATUS=$?
printf '%s\n' "${MPIRUN_VERSION_OUTPUT}"
log "Diagnostic 'mpirun --version' exit status: ${MPIRUN_VERSION_STATUS}"
if [[ "${MPIRUN_VERSION_STATUS}" -ne 0 ]]; then
  die "Open MPI integrated tests require /usr/bin/mpirun --version to succeed (exit ${MPIRUN_VERSION_STATUS})."
fi
if [[ "${MPIRUN_VERSION_OUTPUT}" != *"Open MPI"* ]]; then
  die "Open MPI integrated tests require /usr/bin/mpirun to identify itself as Open MPI."
fi

# ATS accounts for the rank capacity of concurrent tests but does not assign
# disjoint CPU affinities to their independent mpirun processes. Disable each
# launcher's local rank binding so Linux can schedule ranks across the
# container's allowed CPU set; this policy does not guarantee exclusive cores.
export OMPI_MCA_hwloc_base_binding_policy=none
log "Open MPI rank binding disabled with OMPI_MCA_hwloc_base_binding_policy=${OMPI_MCA_hwloc_base_binding_policy}."

# Report the resulting placement policy in each launch's retained stderr.
export OMPI_MCA_hwloc_base_report_bindings=1
log "Open MPI binding reports enabled with OMPI_MCA_hwloc_base_report_bindings=${OMPI_MCA_hwloc_base_report_bindings}."

if [[ "${ATS_MAX_TESTS}" -lt "${NPROC}" ]]; then
  export GEOS_ITS_ATS_MAX_TESTS="${ATS_MAX_TESTS}"
  export PYTHONPATH="${GEOS_SRC_DIR}/scripts${PYTHONPATH:+:${PYTHONPATH}}"
  log "ATS concurrency patch enabled via ${GEOS_SRC_DIR}/scripts/sitecustomize.py"
fi

log "Run ${ATS_CMD[*]}"
if [[ "${GEOS_ITS_FORCE_LEGACY:-0}" == 1 || -n "${GEOS_HYPREDRV_FORCE_LEGACY:-}" ]]; then
  export GEOS_HYPREDRV_FORCE_LEGACY="${GEOS_HYPREDRV_FORCE_LEGACY:-1}"
  log "GEOS_HYPREDRV_FORCE_LEGACY=${GEOS_HYPREDRV_FORCE_LEGACY} (legacy hypre path)"
fi
set +e
"${ATS_CMD[@]}"
ATS_STATUS=$?
set -e

if [[ "${GEOS_ITS_GENERATE_BASELINES:-0}" == 1 ]]; then
  [[ -n "${GEOS_ITS_BASELINE_ARCHIVE:-}" ]] || die "GEOS_ITS_BASELINE_ARCHIVE is empty"
  log "Copy run outputs into the baseline tree"
  # ATS -a rebaseline stops the rest of an .ats file when one test has no
  # restart (e.g. np=27 never wrote a file). That dropped later tests in the
  # same file. Copy every restart/hdf5 artifact the run produced instead.
  copy_run_outputs_to_baselines "${ATS_WORKING_DIR}" "${ATS_BASELINE_DIR}"
  mkdir -p /tmp/geos/generated-baselines
  pack_dest="/tmp/geos/generated-baselines/${GEOS_ITS_BASELINE_ARCHIVE}.tar.gz"
  log "Pack ${pack_dest}"
  or_die integratedTests/geos_ats.sh -a pack_baselines --baselineArchiveName "${pack_dest}"
  log "Baselines archive: ${pack_dest}"
  exit 0
fi

if [[ -x bin/geos_ats_log_check && -f integratedTests/TestResults/test_results.ini ]]; then
  log "geos_ats_log_check"
  bin/geos_ats_log_check integratedTests/TestResults/test_results.ini \
    -y /workspace/.integrated_tests.yaml || true
fi

SANITIZER_REPORT_STATUS=0
if [[ "${SANITIZERS}" -eq 1 ]]; then
  sanitizer_report="${SANITIZER_ROOT}/sanitizer-reports.txt"
  sanitizer_report_tmp="${sanitizer_report}.tmp"
  : > "${sanitizer_report_tmp}"
  for sanitizer_input in integratedTests/TestResults "${ATS_WORKING_DIR}"; do
    if [[ -e "${sanitizer_input}" ]]; then
      grep -R -n -I -E 'AddressSanitizer|UndefinedBehaviorSanitizer|LeakSanitizer|runtime error:|detected memory leaks|SUMMARY:.*(leak|leaked)|Direct leak:|Indirect leak:' \
        "${sanitizer_input}" >> "${sanitizer_report_tmp}" || true
    fi
  done
  for sanitizer_input in "${SANITIZER_ROOT}"/asan.* "${SANITIZER_ROOT}"/lsan.* "${SANITIZER_ROOT}"/ubsan.*; do
    if [[ -f "${sanitizer_input}" ]]; then
      grep -n -I -E 'AddressSanitizer|UndefinedBehaviorSanitizer|LeakSanitizer|runtime error:|detected memory leaks|SUMMARY:.*(leak|leaked)|Direct leak:|Indirect leak:' \
        "${sanitizer_input}" >> "${sanitizer_report_tmp}" || true
    fi
  done
  sort -u "${sanitizer_report_tmp}" -o "${sanitizer_report_tmp}"
  mv "${sanitizer_report_tmp}" "${sanitizer_report}"
  if [[ -s "${sanitizer_report}" ]]; then
    SANITIZER_REPORT_STATUS=1
    log "Sanitizer diagnostics found; full report: ${sanitizer_report}"
    sed -n '1,160p' "${sanitizer_report}"
  else
    log "No ASan/UBSan/LSan diagnostics found; report: ${sanitizer_report}"
  fi
fi

harvest_hypre_logs() {
  local label="$1"
  local dest="${GEOS_SRC_DIR}/.integrated-test-iterations/${label}.json"
  local data_dir="integratedTests/TestResults/test_data"
  if [[ ! -d "${data_dir}" ]]; then
    log "No ${data_dir} to harvest (${label})"
    return 0
  fi
  mkdir -p "$(dirname "${dest}")"
  python3 "${GEOS_SRC_DIR}/scripts/compareLinearSolverIterations.py" harvest \
    "${data_dir}" --geosx-only --strip-ats-prefix -o "${dest}"
}

HYPRE_COMPARE_STATUS=0
if [[ "${GEOS_ITS_FORCE_LEGACY:-0}" == 1 ]]; then
  LEGACY_HARVEST_STATUS=0
  harvest_hypre_logs legacy || LEGACY_HARVEST_STATUS=$?
  HYPREDRIVE_ITERS="${GEOS_SRC_DIR}/.integrated-test-iterations/hypredrive.json"
  LEGACY_ITERS="${GEOS_SRC_DIR}/.integrated-test-iterations/legacy.json"
  if [[ "${LEGACY_HARVEST_STATUS}" -ne 0 ]]; then
    log "Hypre legacy iteration harvest failed"
    HYPRE_COMPARE_STATUS=1
  elif [[ -f "${HYPREDRIVE_ITERS}" && -f "${LEGACY_ITERS}" ]]; then
    log "Comparing all Hypre linear-solver iteration profiles"
    python3 "${GEOS_SRC_DIR}/scripts/compareLinearSolverIterations.py" compare \
      "${HYPREDRIVE_ITERS}" "${LEGACY_ITERS}" --require-identical-keys \
      --total-rel-tol 0.01 --total-abs-tol 2 --max-solve-abs-tol 1 \
      || HYPRE_COMPARE_STATUS=$?
    if [[ "${HYPRE_COMPARE_STATUS}" -ne 0 ]]; then
      log "hypredrive vs legacy Hypre iteration profiles do not match"
    fi
  else
    log "Missing Hypre iteration harvest; both ${HYPREDRIVE_ITERS} and ${LEGACY_ITERS} are required"
    HYPRE_COMPARE_STATUS=1
  fi
  HYPREDRIVE_RESULTS="${GEOS_SRC_DIR}/.integrated-test-iterations/hypredrive-results.ini"
  if [[ -f "${HYPREDRIVE_RESULTS}" && -f integratedTests/TestResults/test_results.ini ]]; then
    log "Comparing ATS result classifications for the hypredrive and legacy passes"
    python3 "${GEOS_SRC_DIR}/scripts/compareIntegratedTestResults.py" \
      "${HYPREDRIVE_RESULTS}" integratedTests/TestResults/test_results.ini \
      || HYPRE_COMPARE_STATUS=$?
    cp integratedTests/TestResults/test_results.ini \
      "${GEOS_SRC_DIR}/.integrated-test-iterations/legacy-results.ini"
  else
    log "Missing ATS result harvest; both ${HYPREDRIVE_RESULTS} and legacy test_results.ini are required"
    HYPRE_COMPARE_STATUS=1
  fi
else
  harvest_hypre_logs hypredrive || true
  if [[ -f integratedTests/TestResults/test_results.ini ]]; then
    mkdir -p "${GEOS_SRC_DIR}/.integrated-test-iterations"
    cp integratedTests/TestResults/test_results.ini \
      "${GEOS_SRC_DIR}/.integrated-test-iterations/hypredrive-results.ini"
  fi
fi

log "ATS exit status: ${ATS_STATUS}"
log "Results (in container): ${BUILD_DIR}/integratedTests/TestResults"
if [[ "${GEOS_ITS_BUILD_DIR:-build-integrated-tests}" == /* ]]; then
  RESULTS_ON_HOST="${GEOS_ITS_BUILD_DIR}/integratedTests/TestResults"
else
  RESULTS_ON_HOST="${GEOS_ITS_HOST_SRC:-.}/${GEOS_ITS_BUILD_DIR:-build-integrated-tests}/integratedTests/TestResults"
fi
log "Results (on host):      ${RESULTS_ON_HOST}"
if [[ "${SANITIZER_REPORT_STATUS}" -ne 0 && "${ATS_STATUS}" -eq 0 ]]; then
  exit "${SANITIZER_REPORT_STATUS}"
fi
if [[ "${HYPRE_COMPARE_STATUS}" -ne 0 && "${ATS_STATUS}" -eq 0 ]]; then
  exit "${HYPRE_COMPARE_STATUS}"
fi
exit "${ATS_STATUS}"
