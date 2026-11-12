#!/usr/bin/env bash
set -euo pipefail

# Create a GEOS MPM-scoped review bundle from a git branch/ref.
#
# Defaults:
#   review branch: mpm
#   scope:         wide
#   output:        /tmp/geos-mpm-review-<branch>-<scope>-<timestamp>.zip
#
# Examples:
#   scripts/make_mpm_review_bundle.sh
#   scripts/make_mpm_review_bundle.sh --dev
#   scripts/make_mpm_review_bundle.sh HEAD --minimal --list-files
#   scripts/make_mpm_review_bundle.sh HEAD --minimal -o /tmp/geos-mpm-minimal.zip
#   scripts/make_mpm_review_bundle.sh mpmdev /tmp/geos-mpmdev-review.zip --wide
#   scripts/make_mpm_review_bundle.sh my-feature --reference mpm --lean
#   scripts/make_mpm_review_bundle.sh my-feature --reference mpmdev --wide -o /tmp/my-feature-review.zip

DEFAULT_RELEASE_BRANCH="mpm"
DEFAULT_DEV_BRANCH="mpmdev"

REVIEW_INPUT="$DEFAULT_RELEASE_BRANCH"
REFERENCE_INPUT=""
SCOPE="wide"
OUT=""
INCLUDE_WORKTREE="auto"
INCLUDE_UNTRACKED="0"
LIST_FILES="0"
POSITIONAL=()

usage() {
  cat <<'USAGE'
Usage:
  scripts/make_mpm_review_bundle.sh [BRANCH] [OUTPUT] [options]

Create a zip bundle containing MPM-relevant GEOS files from BRANCH.
The default BRANCH is mpm. The default scope is --wide.

Branch selection:
  BRANCH                         Branch/ref to review. Default: mpm
  --branch BRANCH, -b BRANCH      Same as positional BRANCH
  --release                       Review mpm
  --dev                           Review mpmdev

Comparison metadata:
  --reference REF, --base REF, -r REF
                                 Include an MPM-scoped diff against REF.
                                 REF may be mpm, mpmdev, a local branch,
                                 a remote branch, or a commit SHA.

Scopes:
  --minimal                       Implementation context for MPM solver and solid
                                  constitutive development, core PFW modules,
                                  three small PFW inputs, and one XML smoke case.
                                  No full test/docs/example trees or repo inventory.
  --lean                          Solver, constitutive models, particle mesh,
                                  MPM events, setup/host-configs, MPM inputs,
                                  and the full particle-file-writer tree.
  --wide                          Everything normally needed for an MPM code
                                  review: lean scope plus mesh/spatial partition,
                                  Silo/VTK/restart output, XML/schema/parsing,
                                  field specifications, table functions, and
                                  related tests/docs. Default.

Output:
  --output FILE, -o FILE           Output archive path. Default is /tmp/*.zip
  --list-files                    Print the selected file paths and exit without
                                  creating an archive. Honors worktree options.

Worktree handling:
  --include-worktree               Overlay staged/modified MPM-relevant files
                                  from the current worktree when reviewing HEAD.
  --no-worktree                    Do not include worktree overlays.
  --include-untracked              With --include-worktree, also overlay
                                  untracked MPM-relevant files.

Notes:
  - --lean and --wide include the full tracked particle-file-writer tree.
    --minimal includes only selected implementation files and sample inputs.
  - --minimal is reading/editing context, not a standalone build checkout.
    See scripts/mpm_minimal_bundle.md for the file groups and exclusions.
  - The script uses git archive, so it can bundle branches that are not checked out.
  - If BRANCH or REF is not found locally, origin/BRANCH or origin/REF is tried.
USAGE
}

die() {
  echo "ERROR: $*" >&2
  exit 1
}

need_arg() {
  local opt="$1"
  local val="${2:-}"
  [ -n "$val" ] || die "$opt requires an argument"
}

while [ "$#" -gt 0 ]; do
  case "$1" in
    --branch|-b)
      need_arg "$1" "${2:-}"
      REVIEW_INPUT="$2"
      shift 2
      ;;
    --release)
      REVIEW_INPUT="$DEFAULT_RELEASE_BRANCH"
      shift
      ;;
    --dev)
      REVIEW_INPUT="$DEFAULT_DEV_BRANCH"
      shift
      ;;
    --reference|--base|--compare|-r)
      need_arg "$1" "${2:-}"
      REFERENCE_INPUT="$2"
      shift 2
      ;;
    --minimal)
      SCOPE="minimal"
      shift
      ;;
    --lean)
      SCOPE="lean"
      shift
      ;;
    --wide)
      SCOPE="wide"
      shift
      ;;
    --output|-o)
      need_arg "$1" "${2:-}"
      OUT="$2"
      shift 2
      ;;
    --list-files)
      LIST_FILES="1"
      shift
      ;;
    --include-worktree)
      INCLUDE_WORKTREE="1"
      shift
      ;;
    --no-worktree)
      INCLUDE_WORKTREE="0"
      shift
      ;;
    --include-untracked)
      INCLUDE_UNTRACKED="1"
      INCLUDE_WORKTREE="1"
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    --)
      shift
      while [ "$#" -gt 0 ]; do
        POSITIONAL+=("$1")
        shift
      done
      ;;
    -*)
      die "unknown option: $1"
      ;;
    *)
      POSITIONAL+=("$1")
      shift
      ;;
  esac
done

if [ "${#POSITIONAL[@]}" -gt 0 ]; then
  REVIEW_INPUT="${POSITIONAL[0]}"
fi
if [ "${#POSITIONAL[@]}" -gt 1 ]; then
  OUT="${POSITIONAL[1]}"
fi
if [ "${#POSITIONAL[@]}" -gt 2 ]; then
  die "too many positional arguments"
fi

if ! git rev-parse --show-toplevel >/dev/null 2>&1; then
  die "run this script from inside a GEOS git checkout"
fi

ROOT="$(git rev-parse --show-toplevel)"
cd "$ROOT"

alias_ref_name() {
  case "$1" in
    release) echo "$DEFAULT_RELEASE_BRANCH" ;;
    dev|development) echo "$DEFAULT_DEV_BRANCH" ;;
    *) echo "$1" ;;
  esac
}

resolve_ref() {
  local raw aliased candidate
  raw="$1"
  aliased="$(alias_ref_name "$raw")"

  # HEAD means this checkout, even when origin/HEAD names a different branch.
  if [ "$aliased" = "HEAD" ]; then
    git rev-parse --verify 'HEAD^{commit}'
    return
  fi

  # Prefer fully-qualified ref names so branch names such as "mpm" do not
  # collide with paths such as the repo-root ./mpm convenience symlink.
  local candidates=(
    "refs/heads/$aliased"
    "refs/remotes/$aliased"
    "refs/remotes/origin/$aliased"
    "refs/tags/$aliased"
  )
  for candidate in "${candidates[@]}"; do
    if git rev-parse --verify --quiet "${candidate}^{commit}" >/dev/null; then
      echo "$candidate"
      return 0
    fi
  done

  # Fall back to arbitrary commit-ish input, including a raw SHA. Print the
  # resolved SHA rather than the raw token to keep later git calls unambiguous.
  if candidate="$(git rev-parse --verify --quiet "${aliased}^{commit}" 2>/dev/null)"; then
    echo "$candidate"
    return 0
  fi

  die "could not resolve git ref '$raw' or 'origin/$aliased'"
}

REVIEW_REF="$(resolve_ref "$REVIEW_INPUT")"
REVIEW_SHA="$(git rev-parse --verify "${REVIEW_REF}^{commit}")"
CURRENT_HEAD="$(git rev-parse HEAD 2>/dev/null || true)"
CURRENT_BRANCH="$(git branch --show-current 2>/dev/null || true)"

REFERENCE_REF=""
REFERENCE_SHA=""
DIFF_RANGE=""
DIFF_RANGE_KIND=""
if [ -n "$REFERENCE_INPUT" ]; then
  REFERENCE_REF="$(resolve_ref "$REFERENCE_INPUT")"
  REFERENCE_SHA="$(git rev-parse --verify "${REFERENCE_REF}^{commit}")"
  if git merge-base "$REFERENCE_SHA" "$REVIEW_SHA" >/dev/null 2>&1; then
    DIFF_RANGE="${REFERENCE_SHA}...${REVIEW_SHA}"
    DIFF_RANGE_KIND="merge-base three-dot diff"
  else
    DIFF_RANGE="${REFERENCE_SHA}..${REVIEW_SHA}"
    DIFF_RANGE_KIND="two-dot diff; no merge base found"
  fi
fi

if [ "$INCLUDE_WORKTREE" = "auto" ]; then
  if [ -n "$CURRENT_HEAD" ] && [ "$CURRENT_HEAD" = "$REVIEW_SHA" ]; then
    INCLUDE_WORKTREE="1"
  else
    INCLUDE_WORKTREE="0"
  fi
fi

safe_ref="$(printf '%s' "$REVIEW_INPUT" | tr '/: @' '____' | tr -cd '[:alnum:]_.-')"
[ -n "$safe_ref" ] || safe_ref="branch"
if [ -z "$OUT" ]; then
  OUT="/tmp/geos-mpm-review-${safe_ref}-${SCOPE}-$(date +%Y%m%d-%H%M%S).zip"
fi
# zip runs from the temporary directory; resolve relative output paths first.
case "$OUT" in
  /*) ;;
  *) OUT="$ROOT/$OUT" ;;
esac

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
BUNDLE="$TMP/geos-mpm-review"
mkdir -p "$BUNDLE/meta" "$BUNDLE/files"

# Filter a list of git paths on stdin down to the selected MPM review scope.
filter_mpm_paths() {
  awk -v scope="$SCOPE" '
    function pfw(p) {
      return p ~ /^scripts\/preProcessing\/materialPointMethodParticleFileWriter(\/|$)/ ||
             p ~ /^scripts\/preprocessing\/materialPointMethodParticleFileWriter(\/|$)/
    }
    function minimal(p, rel) {
      # Keep the writer and its support code, not every input, plot or suite.
      if (pfw(p)) {
        rel = p
        sub(/^scripts\/pre[Pp]rocessing\/materialPointMethodParticleFileWriter\//, "", rel)
        if (rel ~ /^(particleFileWriter|pfw_(geometryObjects|materials|axisymmetric|reference_geometry|simplex_clipping|surface_refinement|particle_layout|tracerPoints|check|inputFixer|material_point)|userDefs_username)\.py$/) return 1
        if (rel == "runClean_username.sh") return 1
        if (rel ~ /^material_point_driver\/(__init__|pfw_material_point)\.py$/) return 1
        return rel == "examples/unifiedGeometry/pfw_input_unified_geometry.py" ||
               rel == "examples/true2D/pfw_input_true_2d.py" ||
               rel == "examples/axisymmetricCylinder/pfw_input_axisymmetricCylinder.py"
      }
      # One small, complete XML/particle/ATS case illustrates the input contract.
      if (p ~ /^inputFiles\/materialPointMethod\/singleParticle\/(mpm_singleParticle\.xml|mpmParticleFile_singleParticle|solidMechanicsMPM_singleParticle\.ats)$/) return 1
      if (p ~ /^src\/cmake\/(GeosxOptions|GeosxConfig)\.cmake$/) return 1
      if (p ~ /^src\/coreComponents\/(physicsSolvers|physicsSolvers\/solidMechanics|constitutive|events|mesh)\/CMakeLists\.txt$/) return 1
      if (p ~ /^src\/coreComponents\/physicsSolvers\/PhysicsSolverBase\.(cpp|hpp)$/) return 1
      if (p ~ /^src\/coreComponents\/physicsSolvers\/solidMechanics\/(SolidMechanicsMPM|MPM[^/]*|SolidMechanicsFields|LogLevelsInfo)\.(cpp|hpp)$/) return 1
      if (p ~ /^src\/coreComponents\/physicsSolvers\/solidMechanics\/kernels\/(ExplicitMPM|ExplicitFiniteStrain|ExplicitSmallStrain)\.hpp$/) return 1
      if (p ~ /^src\/coreComponents\/finiteElement\/(Kinematics\.h|FiniteElementBase\.hpp|FiniteElementDiscretization\.hpp|FiniteElementDiscretizationManager\.hpp|finiteElementBaseKernels\.hpp)$/) return 1
      if (p == "src/coreComponents/finiteElement/kernelInterface/KernelBase.hpp") return 1
      if (p == "src/coreComponents/finiteElement/elementFormulations/FiniteElementBase.hpp") return 1
      # Top-level solid implementations and registration; no porosity/test/docs trees.
      if (p ~ /^src\/coreComponents\/constitutive\/(ConstitutiveBase|ConstitutiveManager|ConstitutivePassThru|ConstitutivePassThruHandler|ContinuumBase|NullModel|ExponentialRelation)\.(cpp|hpp)$/) return 1
      if (p == "src/coreComponents/constitutive/ExternalConstitutiveModels.hpp.in") return 1
      if (p ~ /^src\/coreComponents\/constitutive\/(solid|cohesiveZone|contact|gas|liquid)\/[^/]+\.(cpp|hpp)$/) return 1
      if (p ~ /^src\/coreComponents\/constitutive\/materialPointDriver\/(CMakeLists\.txt|[^/]+\.(cpp|hpp))$/) return 1
      if (p ~ /^src\/coreComponents\/events\/mpmEvents\/[^/]+\.(cpp|hpp)$/) return 1
      if (p ~ /^src\/coreComponents\/mesh\/(Particle[^/]*|ToParticleRelation|CohesiveZone[^/]*|AxisymmetricGeometry)\.(cpp|hpp)$/) return 1
      if (p ~ /^src\/coreComponents\/mesh\/particleGenerators\/[^/]+\.(cpp|hpp)$/) return 1
      # Interfaces keep field registration, mesh access and MPI calls readable.
      if (p ~ /^src\/coreComponents\/mesh\/(DomainPartition|MeshBody|MeshLevel|ObjectManagerBase|NodeManager|MeshFields)\.hpp$/) return 1
      if (p ~ /^src\/coreComponents\/mesh\/mpiCommunications\/(CommunicationTools|MPI_iCommData|NeighborCommunicator|PartitionBase|SpatialPartition)\.hpp$/) return 1
      if (p ~ /^src\/coreComponents\/mesh\/utilities\/ComputationalGeometry\.hpp$/) return 1
      if (p ~ /^src\/coreComponents\/dataRepository\/(Group|Wrapper|WrapperBase|ExecutableGroup|RestartFlags|InputFlags|ReferenceWrapper|xmlWrapper)\.hpp$/) return 1
      if (p ~ /^src\/coreComponents\/common\/(DataTypes|GeosxMacros|GEOS_RAJA_Interface|TypeDispatch|TimingMacros|MpiWrapper)\.hpp$/) return 1
      if (p == "src/coreComponents/common/format/EnumStrings.hpp") return 1
      if (p ~ /^src\/coreComponents\/fieldSpecification\/(FieldSpecificationBase|FieldSpecificationManager|TractionBoundaryCondition)\.hpp$/) return 1
      if (p ~ /^src\/coreComponents\/functions\/(FunctionBase|FunctionManager|TableFunction)\.hpp$/) return 1
      return 0
    }
    function lean(p, lp) {
      lp = tolower(p)
      if (pfw(p)) return 1
      if (p == "setupMPM" || p == ".gitignore" || p == ".gitmodules") return 1
      if (p == "scripts/make_mpm_review_bundle.sh" || p == "scripts/mpm_codex_review_prompt.txt") return 1
      if (p ~ /^host-configs\/profiles\/mpm/) return 1
      if (p ~ /^host-configs\// && lp ~ /mpm/) return 1
      if (p ~ /^src\/cmake\/(GeosxOptions|GeosxConfig)\.cmake$/) return 1
      if (p ~ /^src\/cmake\/thirdparty\/SetupGeosxThirdParty\.cmake$/) return 1
      if (p ~ /^src\/coreComponents\/physicsSolvers\/solidMechanics\// && p ~ /(MPM|MaterialPoint|ExplicitMPM|CMakeLists\.txt)/) return 1
      if (p ~ /^src\/coreComponents\/physicsSolvers\/solidMechanics\/kernels\// && p ~ /(MPM|MaterialPoint)/) return 1
      if (p ~ /^src\/coreComponents\/events\/(CMakeLists\.txt|mpmEvents\/)/) return 1
      if (p ~ /^src\/coreComponents\/mesh\// && p ~ /(CMakeLists\.txt|Particle|particleGenerators|ToParticle)/) return 1
      if (p ~ /^src\/coreComponents\/constitutive\/(CMakeLists\.txt|[^\/]+\.cpp|[^\/]+\.hpp|ExternalConstitutiveModels\.hpp\.in)/) return 1
      if (p ~ /^src\/coreComponents\/constitutive\/(solid|cohesiveZone|contact)\//) return 1
      if (p ~ /^src\/coreComponents\/constitutive\/gas\/Gas\./) return 1
      if (p ~ /^src\/coreComponents\/constitutive\/docs\/(ExternalMaterialModels\.rst|solid\/)/) return 1
      if (p ~ /^inputFiles\/materialPointMethod\//) return 1
      if (p ~ /^src\/coreComponents\/schema\/docs\/.*(MPM|Particle|particleRegionsGroup|Constitutive|Solid|Elastic|Plastic|Damage|Cohesive|Gas).*\.rst$/) return 1
      return 0
    }
    function wide(p, lp) {
      lp = tolower(p)
      if (lean(p)) return 1
      if (p ~ /^src\/coreComponents\/physicsSolvers\/solidMechanics\//) return 1
      if (p ~ /^src\/coreComponents\/events\//) return 1
      if (p ~ /^src\/coreComponents\/mesh\//) return 1
      if (p ~ /^src\/coreComponents\/fileIO\//) return 1
      if (p ~ /^src\/coreComponents\/constitutive\//) return 1
      if (p ~ /^src\/coreComponents\/dataRepository\// && p ~ /(CMakeLists\.txt|ConduitRestart|RestartFlags|xmlWrapper|DataContext|Group|Wrapper|Utilities|BufferOps|LogLevels|ExecutableGroup)/) return 1
      if (p ~ /^src\/coreComponents\/schema\//) return 1
      if (p ~ /^src\/coreComponents\/codingUtilities\/(Parsing|RTTypes|CMakeLists)/) return 1
      if (p ~ /^src\/coreComponents\/fieldSpecification\//) return 1
      if (p ~ /^src\/coreComponents\/functions\//) return 1
      if (p ~ /^src\/coreComponents\/common\// && p ~ /(CMakeLists|DataTypes|FieldSpecificationOps|MpiWrapper|Path|Units|Logger|ErrorHandling|format\/|initializeEnvironment)/) return 1
      if (p ~ /^src\/coreComponents\/integrationTests\// && p ~ /(MPM|Particle|meshTests|xmlTests|dataRepositoryTests|fieldSpecificationTests|tableFunctionsFileTests)/) return 1
      if (p ~ /^src\/docs\/sphinx\// && lp ~ /(mpm|materialpoint|material point|particle|xml|output|spatialpartition|constitutive)/) return 1
      if (p ~ /^examples\// && lp ~ /(mpm|materialpoint|material point|particle)/) return 1
      return 0
    }
    {
      if (scope == "minimal") {
        if (minimal($0)) print $0
      } else if (scope == "lean") {
        if (lean($0)) print $0
      } else {
        if (wide($0)) print $0
      }
    }
  '
}

git ls-tree -r --name-only "$REVIEW_SHA" > "$TMP/all-branch-files.txt"
filter_mpm_paths < "$TMP/all-branch-files.txt" | sort -u > "$TMP/review-files.txt"

# Plan the overlay once so --list-files and the archive have identical contents.
: > "$BUNDLE/meta/worktree-files-mpm.txt"
if [ "$INCLUDE_WORKTREE" = "1" ] && [ "$CURRENT_HEAD" = "$REVIEW_SHA" ]; then
  git diff --name-only HEAD > "$TMP/worktree-files.txt"
  if [ "$INCLUDE_UNTRACKED" = "1" ]; then
    git ls-files --others --exclude-standard >> "$TMP/worktree-files.txt"
  fi
  filter_mpm_paths < "$TMP/worktree-files.txt" | sort -u > "$BUNDLE/meta/worktree-files-mpm.txt"
fi
if [ "$LIST_FILES" = "1" ]; then
  {
    comm -23 "$TMP/review-files.txt" "$BUNDLE/meta/worktree-files-mpm.txt"
    while IFS= read -r path; do
      if [ -f "$path" ]; then
        printf '%s\n' "$path"
      fi
    done < "$BUNDLE/meta/worktree-files-mpm.txt"
  } | sort -u
  exit 0
fi

# Comparison metadata. The all-file diff names are included for transparency,
# but the main patch is MPM-scoped to avoid wasting review time.
: > "$TMP/changed-files-all.txt"
: > "$TMP/changed-files-mpm.txt"
if [ -n "$DIFF_RANGE" ]; then
  if [ "$SCOPE" != "minimal" ]; then
    git diff --stat "$DIFF_RANGE" > "$BUNDLE/meta/diff-stat-all.txt" || true
    git diff --name-status "$DIFF_RANGE" > "$BUNDLE/meta/diff-name-status-all.txt" || true
  fi
  git diff --name-only "$DIFF_RANGE" | sort -u > "$TMP/changed-files-all.txt" || true
  filter_mpm_paths < "$TMP/changed-files-all.txt" | sort -u > "$TMP/changed-files-mpm.txt"

  : > "$BUNDLE/meta/diff-name-status-mpm.txt"
  : > "$BUNDLE/meta/diff-stat-mpm.txt"
  : > "$BUNDLE/meta/diff-mpm.patch"
  while IFS= read -r path; do
    [ -n "$path" ] || continue
    git diff --name-status "$DIFF_RANGE" -- "$path" >> "$BUNDLE/meta/diff-name-status-mpm.txt" || true
    git diff --stat "$DIFF_RANGE" -- "$path" >> "$BUNDLE/meta/diff-stat-mpm.txt" || true
    git diff --binary "$DIFF_RANGE" -- "$path" >> "$BUNDLE/meta/diff-mpm.patch" || true
  done < "$TMP/changed-files-mpm.txt"
  sort -u "$BUNDLE/meta/diff-name-status-mpm.txt" -o "$BUNDLE/meta/diff-name-status-mpm.txt"

  if [ "$SCOPE" != "minimal" ]; then
    comm -23 "$TMP/changed-files-all.txt" "$TMP/changed-files-mpm.txt" > "$BUNDLE/meta/diff-files-outside-mpm-scope.txt" || true
  fi
else
  if [ "$SCOPE" != "minimal" ]; then
    : > "$BUNDLE/meta/diff-stat-all.txt"
    : > "$BUNDLE/meta/diff-name-status-all.txt"
    : > "$BUNDLE/meta/diff-files-outside-mpm-scope.txt"
  fi
  : > "$BUNDLE/meta/diff-name-status-mpm.txt"
  : > "$BUNDLE/meta/diff-stat-mpm.txt"
  : > "$BUNDLE/meta/diff-mpm.patch"
fi

cp "$TMP/review-files.txt" "$BUNDLE/meta/review-files.txt"
cp "$TMP/changed-files-mpm.txt" "$BUNDLE/meta/diff-files-mpm.txt"

if [ "$SCOPE" != "minimal" ]; then
  cp "$TMP/all-branch-files.txt" "$BUNDLE/meta/branch-files-all.txt"
  cp "$TMP/changed-files-all.txt" "$BUNDLE/meta/diff-files-all.txt"
  git status --short > "$BUNDLE/meta/current-worktree-status-short.txt" || true
  git submodule status --recursive > "$BUNDLE/meta/current-submodule-status.txt" 2>/dev/null || true
  git log --oneline --decorate -n 20 "$REVIEW_SHA" -- > "$BUNDLE/meta/review-branch-log.txt" || true
fi

cat > "$BUNDLE/meta/README.txt" <<INFO
GEOS MPM review bundle
Created: $(date)
Repository: $ROOT
Review input: $REVIEW_INPUT
Resolved review ref: $REVIEW_REF
Review commit: $REVIEW_SHA
Current checkout branch: ${CURRENT_BRANCH:-detached-or-unknown}
Current checkout HEAD: ${CURRENT_HEAD:-unknown}
Scope: $SCOPE
Reference input: ${REFERENCE_INPUT:-none}
Resolved reference ref: ${REFERENCE_REF:-none}
Reference commit: ${REFERENCE_SHA:-none}
Diff range: ${DIFF_RANGE:-none}
Diff range kind: ${DIFF_RANGE_KIND:-none}
Worktree overlay: $INCLUDE_WORKTREE
Untracked overlay: $INCLUDE_UNTRACKED

Start review with:
  meta/codex_mpm_review_prompt.txt
  meta/review-files.txt
  meta/diff-mpm.patch
  files/

The minimal scope includes selected PFW implementation files and three sample
inputs. Lean and wide include the full tracked PFW tree.
INFO

# Extract selected tracked files from the requested branch/ref.
if [ -s "$TMP/review-files.txt" ]; then
  mapfile -t review_paths < "$TMP/review-files.txt"
  git archive --format=tar "$REVIEW_SHA" -- "${review_paths[@]}" | tar -xf - -C "$BUNDLE/files"
fi

# Optionally overlay staged/modified files from the current worktree when the
# reviewed ref is the current HEAD. This lets developers review local edits
# without committing them first, while keeping branch snapshots clean by default.
if [ "$INCLUDE_WORKTREE" = "1" ]; then
  if [ -n "$CURRENT_HEAD" ] && [ "$CURRENT_HEAD" = "$REVIEW_SHA" ]; then
    # The minimal bundle already contains the overlaid source. Avoid duplicate
    # patches and unrelated worktree changes in this scope.
    if [ "$SCOPE" != "minimal" ]; then
      git diff --binary HEAD > "$BUNDLE/meta/worktree-unstaged.patch" || true
      git diff --binary --cached HEAD > "$BUNDLE/meta/worktree-staged.patch" || true
    fi

    while IFS= read -r path; do
      [ -n "$path" ] || continue
      if [ ! -f "$path" ]; then
        rm -f "$BUNDLE/files/$path"
        continue
      fi
      mkdir -p "$BUNDLE/files/$(dirname "$path")"
      cp -p "$path" "$BUNDLE/files/$path"
    done < "$BUNDLE/meta/worktree-files-mpm.txt"
  else
    echo "Worktree overlay requested but skipped: review ref is not current HEAD." > "$BUNDLE/meta/worktree-overlay-skipped.txt"
    : > "$BUNDLE/meta/worktree-files-mpm.txt"
  fi
else
  : > "$BUNDLE/meta/worktree-files-mpm.txt"
fi

# Include the Codex prompt from the repo if present, otherwise use an embedded copy.
PROMPT_IN_REPO="$ROOT/scripts/mpm_codex_review_prompt.txt"
if [ "$SCOPE" = "minimal" ]; then
  cat > "$BUNDLE/meta/codex_mpm_review_prompt.txt" <<'PROMPT'
Use this bundle to implement the requested GEOS MPM solver or solid constitutive
feature. Start with meta/README.txt, meta/files-copied.txt and meta/file-sizes.txt.
Read only files relevant to the feature; do not load the entire bundle at once.
Source under files/ preserves repository paths and includes worktree overlays
when enabled. meta/diff-mpm.patch is populated only when --reference is supplied.

Entry points (relative to files/):
- src/coreComponents/physicsSolvers/solidMechanics/SolidMechanicsMPM.{hpp,cpp}:
  solver configuration and time stepping; MPMSolverFields.hpp and
  MPMSolverEnums.hpp define fields and options; MPM*.hpp and kernels/ExplicitMPM.hpp
  contain mapping, geometry, surface and constitutive update helpers.
- src/coreComponents/constitutive/: ContinuumBase, solid/SolidBase and
  solid/SolidUtilities define update contracts. ConstitutivePassThru.hpp,
  ExternalConstitutiveModels.hpp.in and CMakeLists.txt control registration/build.
  Read the relevant solid or cohesiveZone model, plus materialPointDriver/ when
  developing or checking a constitutive update.
- src/coreComponents/mesh/: Particle*, particleGenerators/, CohesiveZone* and
  AxisymmetricGeometry.hpp define particle storage, geometry and input formats.
- src/coreComponents/events/mpmEvents/: event registration and runtime updates.
- scripts/preProcessing/materialPointMethodParticleFileWriter/: particleFileWriter.py
  forwards solver/material options to XML; pfw_geometryObjects.py, pfw_materials.py
  and the geometry/refinement modules define particle generation. Three small
  examples cover Cartesian, plane-strain and axisymmetric inputs. The
  pfw_material_point.py entry point and material_point_driver/ support material tests.
- inputFiles/materialPointMethod/singleParticle/: one XML/particle/ATS smoke case.

This is implementation context, not a standalone build checkout. General framework
code is represented mostly by interfaces. Full tests, documentation, other inputs,
output backends, third-party code and externally configured material repositories
are omitted. Ask for specific omitted files if the feature depends on them.
Do not infer missing APIs or claim compilation/testing from this bundle alone.
Supply a focused patch and a feature-specific validation plan.
PROMPT
elif [ -f "$PROMPT_IN_REPO" ]; then
  cp "$PROMPT_IN_REPO" "$BUNDLE/meta/codex_mpm_review_prompt.txt"
else
  cat > "$BUNDLE/meta/codex_mpm_review_prompt.txt" <<'PROMPT'
You are reviewing a GEOS MPM-scoped bundle. Start with meta/README.txt,
meta/review-files.txt, and meta/diff-mpm.patch. Limit your review to MPM-relevant
code and context included under files/.

Focus on setupMPM, MPM host-configs, MPM solver code, MPM events, particle mesh
infrastructure, constitutive models used by MPM, MPM XML/input/schema handling,
field specifications, table functions, spatial partition/communication touched by
MPM, and Silo/VTK/restart output paths used by MPM. The entire
scripts/preProcessing/materialPointMethodParticleFileWriter tree is in scope.

Do not spend time reviewing unrelated GEOS subsystems unless they are changed in
meta/diff-mpm.patch, directly referenced by an MPM code path, or required to
explain a build/runtime error. Report findings with path, line number when
available, severity, why it matters, and a concrete fix.
PROMPT
fi

find "$BUNDLE/files" -type f | sed "s#^$BUNDLE/files/##" | sort > "$BUNDLE/meta/files-copied.txt"
if [ "$SCOPE" = "minimal" ]; then
  while IFS= read -r path; do
    printf '%12d %s\n' "$(wc -c < "$BUNDLE/files/$path")" "$path"
  done < "$BUNDLE/meta/files-copied.txt" | sort -nr > "$BUNDLE/meta/file-sizes.txt"
fi
{
  echo "review files: $(wc -l < "$BUNDLE/meta/review-files.txt")"
  echo "copied files:  $(wc -l < "$BUNDLE/meta/files-copied.txt")"
  echo "MPM changed files in comparison: $(wc -l < "$BUNDLE/meta/diff-files-mpm.txt")"
  echo "all changed files in comparison: $(wc -l < "$TMP/changed-files-all.txt")"
  echo "worktree MPM files overlaid: $(wc -l < "$BUNDLE/meta/worktree-files-mpm.txt")"
  if [ "$SCOPE" = "minimal" ]; then
    awk '{ bytes += $1 } END {
      printf "source bytes: %.0f\n", bytes
      printf "rough source tokens (bytes / 4; not tokenizer-measured): %.0f\n", bytes / 4
    }' "$BUNDLE/meta/file-sizes.txt"
  fi
} > "$BUNDLE/meta/file-counts.txt"

mkdir -p "$(dirname "$OUT")"
case "$OUT" in
  *.zip)
    if command -v zip >/dev/null 2>&1; then
      # A fresh zip prevents files from an older, wider bundle surviving a rerun.
      ( cd "$TMP" && zip -qr "$TMP/bundle.zip" geos-mpm-review )
      mv "$TMP/bundle.zip" "$OUT"
    elif command -v python3 >/dev/null 2>&1; then
      python3 - "$OUT" "$BUNDLE" <<'PY'
import os
import sys
import zipfile
out = sys.argv[1]
root = sys.argv[2]
base = os.path.dirname(root)
with zipfile.ZipFile(out, "w", compression=zipfile.ZIP_DEFLATED) as zf:
    for dirpath, dirnames, filenames in os.walk(root):
        for name in filenames:
            path = os.path.join(dirpath, name)
            arcname = os.path.relpath(path, base)
            zf.write(path, arcname)
PY
    else
      die "neither zip nor python3 is available; choose an output ending in .tar.gz"
    fi
    ;;
  *.tar.gz|*.tgz)
    tar -C "$TMP" -czf "$OUT" geos-mpm-review
    ;;
  *)
    die "unsupported output extension '$OUT'; use .zip or .tar.gz"
    ;;
esac

echo "$OUT"
