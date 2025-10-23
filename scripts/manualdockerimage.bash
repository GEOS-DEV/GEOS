#!/usr/bin/env bash
set -euo pipefail

# --------- CONFIG (edit or pass as env vars) -----------------
: "${APP_IMAGE_REPO:=geosx/geosx_ubuntu22.04-clang15}"      # Docker Hub repo to push to
: "${GIT_URL:=git@github.com:GEOS-DEV/GEOS.git}"         # SSH repo URL
: "${GIT_REF:=stdMap}"                              # branch/tag/sha (optional; will try to checkout)
: "${BASE_IMAGE:=geosx/ubuntu22.04-clang15:320-766}"         # TPL image to build in
: "${CMAKE_BUILD_TYPE:=Release}"
: "${HOSTCFG:=/spack-generated.cmake}"
: "${ENABLE_HYPRE:=ON}"
: "${ENABLE_TRILINOS:=OFF}"
: "${GEOS_ENABLE_BOUNDS_CHECK:=OFF}"
: "${PLATFORM:=}"  # e.g. set to linux/amd64 on Apple Silicon
# -------------------------------------------------------------

# sanity: need an ssh-agent socket to forward
if [[ -z "${SSH_AUTH_SOCK:-}" || ! -S "$SSH_AUTH_SOCK" ]]; then
  echo "ERROR: SSH_AUTH_SOCK not set or not a socket. Run:  eval \"\$(ssh-agent -s)\" && ssh-add ~/.ssh/id_ed25519" >&2
  exit 1
fi

NAME="geosx_build_$(date +%s)"

echo "Base image:      ${BASE_IMAGE}"
echo "Repo URL:        ${GIT_URL}"
echo "Git ref:         ${GIT_REF}"
echo "Output registry: ${APP_IMAGE_REPO}"

docker rm -f "$NAME" >/dev/null 2>&1 || true
docker run -it ${PLATFORM:+--platform="${PLATFORM}"} --name "$NAME" --cap-add=SYS_PTRACE \
  -v "$SSH_AUTH_SOCK:/ssh-agent" -e SSH_AUTH_SOCK=/ssh-agent \
  -e ENABLE_HYPRE="${ENABLE_HYPRE}" \
  -e ENABLE_TRILINOS="${ENABLE_TRILINOS}" \
  -e GEOS_BUILD_SHARED_LIBS=ON \
  --volume /Users/settgast1/Codes/geos/GEOS_Develop:/src/geosx \
  "$BASE_IMAGE"
  


  #   # Build & install using your existing CI script
    /src/geosx/scripts/ci_build_and_test_in_container.sh --repository /src/geosx --cmake-build-type Release --host-config /spack-generated.cmake --enable-hypre ON --enable-trilinos OFF --geos-enable-bounds-check OFF --install-dir-basename geosx-install --no-run-unit-tests --nproc 4 --ninja

  #   # Hygiene: drop git metadata so creds aren't baked into the committed image
  #   rm -rf /src/geosx/.git /src/geosx/.gitmodules || true
  # "

# # Grab the short sha that the container wrote
# docker cp "$NAME":/tmp/shortsha ./shortsha || true
# if [ -s ./shortsha ]; then
#   SHORT_SHA="$(cat ./shortsha)"
# else
#   # Fallback if file missing: use the ref or container ID tail
#   SHORT_SHA="${GIT_REF:-$(docker inspect --format='{{.Id}}' "$NAME" | cut -c1-7)}"
# fi
# rm -f ./shortsha

# TAG="${APP_IMAGE_REPO}:${SHORT_SHA}"

# echo "Logging into Docker Hub to push ${TAG} ..."
# docker login
# docker commit --change='ENTRYPOINT []' --change='CMD ["/bin/bash"]' "$NAME" "$TAG"
# docker push "$TAG"
# docker rm "$NAME" || true

# echo "Pushed ${TAG}"
