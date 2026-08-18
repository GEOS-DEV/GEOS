#!/usr/bin/env bash
#
# Run the GEOS integrated tests in the same Docker image as streak2 CI:
#   geosx/ubuntu24.04-gcc12:<GEOS_TPL_TAG>
# where GEOS_TPL_TAG comes from .devcontainer/devcontainer.json.
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
# After ATS, `_iterative` geosx logs are harvested to
# .integrated-test-iterations/{hypredrive,legacy}.json. A `--force-legacy` run
# then compares that harvest to the hypredrive JSON with exact per-solve
# sequences (zero slack) and fails if they differ.
#
# Resource flags apply to the Docker container and to ATS rank limits.
# Remaining arguments are forwarded to geos_ats.sh.
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
BASELINES=""
FILTER=""
CUTOFF="45m"
BUILD_DIR_NAME="build-integrated-tests"
PULL=1
GENERATE_BASELINES=0
FORCE_LEGACY=0
CPUS_FROM_CLI=0
MEMORY_FROM_CLI=0
EXTRA_ATS_ARGS=()

usage()
{
  sed -n '2,/^$/p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
  cat <<EOF
Options:
  --cpus N              ATS max ranks (and Docker CPU limit). Default ${CPUS}.
                        Under rootless Podman/Slurm the CPU cgroup limit is ignored.
  --memory SIZE         Docker memory limit, e.g. 32g, 128g, 256g (default ${MEMORY}).
                        Under rootless Podman/Slurm this cgroup limit is ignored.
                        A bare number is treated as gigabytes (380 -> 380g).
  --nproc N             Compile parallelism (default: same as --cpus)
  --baselines DIR       Host directory of ATS baseline tarballs (mounted read-only).
                        If omitted, download the archive from .integrated_tests.yaml.
                        If the download fails, continue without baseline comparison.
  --generateBaselines   Build develop at 32 cores / 128g, copy every restart
                        and hdf5 artifact the run produced, and pack the
                        tarball named in develop's .integrated_tests.yaml into
                        .integrated-test-baselines/. Cannot be used with --baselines.
  --filter EXPR         ATS name filter (tests whose name contains EXPR)
  --force-legacy        Set GEOS_HYPREDRV_FORCE_LEGACY=1 (legacy hypre path)
  --cutoff TIME         ATS cutoff (default ${CUTOFF})
  --build-dir NAME      Build directory name under the repo (default ${BUILD_DIR_NAME})
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
    --baselines)  BASELINES="$2"; shift 2;;
    --generateBaselines) GENERATE_BASELINES=1; shift;;
    --filter)     FILTER="$2"; shift 2;;
    --force-legacy) FORCE_LEGACY=1; shift;;
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

# ---------------------------------------------------------------------------
# Host: launch the streak2 CI image
# ---------------------------------------------------------------------------
if [[ "${GEOS_IN_CONTAINER:-}" != 1 ]]; then
  if [[ "$(uname -s)" == Darwin ]]; then
    die "This script requires a native Linux host. The geosx/ubuntu24.04-gcc12 image is linux/amd64 and is not supported on macOS."
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
  IMAGE="geosx/ubuntu24.04-gcc12:${TPL_TAG}"

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
    --env GEOS_ITS_CUTOFF="${CUTOFF}"
    --env GEOS_ITS_FILTER="${FILTER}"
    --env GEOS_ITS_BUILD_DIR="${BUILD_DIR_NAME}"
    --env GEOS_ITS_SKIP_BASELINES="${SKIP_PACKED_BASELINES}"
    --env GEOS_ITS_GENERATE_BASELINES="${GENERATE_BASELINES}"
    --env GEOS_ITS_BASELINE_ARCHIVE="${BASELINE_ARCHIVE}"
    --env GEOS_ITS_HOST_SRC="${GEOS_SRC_DIR}"
    --env GEOS_ITS_FORCE_LEGACY="${FORCE_LEGACY}"
    --mount "type=bind,src=${GEOS_SRC_DIR},dst=/workspace"
  )
  if [[ "${FORCE_LEGACY}" -eq 1 ]]; then
    docker_args+=(--env GEOS_HYPREDRV_FORCE_LEGACY=1)
    log "Forcing legacy hypre path (GEOS_HYPREDRV_FORCE_LEGACY=1)"
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
CUTOFF="${GEOS_ITS_CUTOFF:-45m}"
FILTER="${GEOS_ITS_FILTER:-}"
BUILD_DIR="/workspace/${GEOS_ITS_BUILD_DIR:-build-integrated-tests}"
VENV_DIR="${BUILD_DIR}/ats-venv"
VENV_PY="${VENV_DIR}/bin/python3"
ATS_WORKING_DIR="${BUILD_DIR}/ats-working"
ATS_BASELINE_DIR="${BUILD_DIR}/ats-baselines"
HOST_CONFIG=/spack-generated.cmake

if ! [[ "${NPROC}" =~ ^[1-9][0-9]*$ ]]; then
  die "NPROC must be a positive integer, got '${NPROC}'."
fi

# Open MPI refuses root launches in CI containers unless both guard variables
# are set. Keep this separate from openmpi_args because repeated geos-ats
# overrides replace values.
export OMPI_ALLOW_RUN_AS_ROOT=1
export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
ATS_ARGUMENTS="--machine openmpi --ats openmpi_mpirun=/usr/bin/mpirun --ats openmpi_procspernode=${NPROC} --ats openmpi_maxprocs=${NPROC} --ats cutoff=${CUTOFF}"

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

log "Create ATS virtualenv"
if [[ ! -x "${VENV_PY}" ]]; then
  or_die python3 -m venv --system-site-packages "${VENV_DIR}"
fi
git config --global --add safe.directory /workspace 2>/dev/null || true
"${VENV_PY}" -m pip install --upgrade pip setuptools wheel || true
"${VENV_PY}" -m pip install --disable-pip-version-check --upgrade matplotlib || true

log "Configure (ENABLE_HYPREDRV=ON)"
or_die cmake -S /workspace/src -B "${BUILD_DIR}" -G Ninja \
  -C "${HOST_CONFIG}" \
  -DCMAKE_BUILD_TYPE=Release \
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
  -DATS_WORKING_DIR="${ATS_WORKING_DIR}"

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
  ATS_CMD+=(-f "${FILTER}")
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

harvest_iterative_logs() {
  local label="$1"
  local dest="${GEOS_SRC_DIR}/.integrated-test-iterations/${label}.json"
  local data_dir="integratedTests/TestResults/test_data"
  if [[ ! -d "${data_dir}" ]]; then
    log "No ${data_dir} to harvest (${label})"
    return 0
  fi
  mkdir -p "$(dirname "${dest}")"
  python3 "${GEOS_SRC_DIR}/scripts/compareLinearSolverIterations.py" harvest \
    "${data_dir}" --iterative -o "${dest}"
}

ITER_COMPARE_STATUS=0
if [[ "${GEOS_ITS_FORCE_LEGACY:-0}" == 1 ]]; then
  harvest_iterative_logs legacy || true
  HYPREDRIVE_ITERS="${GEOS_SRC_DIR}/.integrated-test-iterations/hypredrive.json"
  LEGACY_ITERS="${GEOS_SRC_DIR}/.integrated-test-iterations/legacy.json"
  if [[ -f "${HYPREDRIVE_ITERS}" && -f "${LEGACY_ITERS}" ]]; then
    log "Comparing _iterative linear-solver iteration sequences (zero slack)"
    python3 "${GEOS_SRC_DIR}/scripts/compareLinearSolverIterations.py" compare \
      "${HYPREDRIVE_ITERS}" "${LEGACY_ITERS}" --exact-sequence || ITER_COMPARE_STATUS=$?
    if [[ "${ITER_COMPARE_STATUS}" -ne 0 ]]; then
      log "hypredrive vs legacy _iterative iteration sequences do not match"
    fi
  else
    log "Skipping iteration compare (need both ${HYPREDRIVE_ITERS} and ${LEGACY_ITERS})"
  fi
else
  harvest_iterative_logs hypredrive || true
fi

log "ATS exit status: ${ATS_STATUS}"
log "Results (in container): ${BUILD_DIR}/integratedTests/TestResults"
log "Results (on host):      ${GEOS_ITS_HOST_SRC:-.}/${GEOS_ITS_BUILD_DIR:-build-integrated-tests}/integratedTests/TestResults"
if [[ "${ITER_COMPARE_STATUS}" -ne 0 && "${ATS_STATUS}" -eq 0 ]]; then
  exit "${ITER_COMPARE_STATUS}"
fi
exit "${ATS_STATUS}"
