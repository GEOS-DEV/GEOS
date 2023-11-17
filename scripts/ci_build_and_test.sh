#!/bin/bash
env

# The linux build relies on two environment variables DOCKER_REPOSITORY and GEOSX_TPL_TAG to define the TPL version.
# And another CMAKE_BUILD_TYPE to define the build type we want for GEOSX.
# Optional BUILD_AND_TEST_ARGS to pass arguments to build_test_helper.sh script.
#
# We extract the location of the GEOSX_TPL from the container...
GEOSX_TPL_DIR=$(docker run --rm ${DOCKER_REPOSITORY}:${GEOSX_TPL_TAG} /bin/bash -c 'echo ${GEOSX_TPL_DIR}' | tail -1)
# ... so we can install GEOSX alongside. This is assumed for bundling the binaries, so consider modifying with care.
GEOSX_DIR=${GEOSX_TPL_DIR}/../GEOSX-${GITHUB_SHA:0:7}
# We need to get the build directory
BUILD_DIR=${GITHUB_WORKSPACE}
# We need to know where the code folder is mounted inside the container so we can run the script at the proper location!
# Since this information is repeated twice, we use a variable.
BUILD_DIR_MOUNT_POINT=/tmp/GEOSX

# We need to keep track of the building container (hence the `CONTAINER_NAME`)
# so we can extract the data from it later (if needed). Another solution would have been to use a mount point,
# but that would not have solved the problem for the TPLs (we would require extra action to copy them to the mount point).
CONTAINER_NAME=geosx_build

dargs=()  # docker arguments

if [ ${USE_SCCACHE} = true ]; then
  dargs+=(--volume=/opt/gcs:/opt/gcs)
  SCCACHE_CLI="--use-sccache"

  echo "Creating the configuration file for sccache..."
  mkdir -p /opt/gcs
  cp -rp ${GOOGLE_GHA_CREDS_PATH} /opt/gcs/credentials.json
fi

dargs+=(--name=${CONTAINER_NAME})
dargs+=(--volume=${BUILD_DIR}:${BUILD_DIR_MOUNT_POINT})
dargs+=(--cap-add=ALL)

dargs+=(-e HOST_CONFIG=${HOST_CONFIG:-host-configs/environment.cmake})
dargs+=(-e GEOSX_DIR=${GEOSX_DIR})
dargs+=(-e ENABLE_HYPRE=${ENABLE_HYPRE:-OFF})
dargs+=(-e ENABLE_HYPRE_DEVICE=${ENABLE_HYPRE_DEVICE:-CPU})
dargs+=(-e ENABLE_TRILINOS=${ENABLE_TRILINOS:-ON})
dargs+=(-e CMAKE_BUILD_TYPE)

# Now we can build GEOSX.
docker run "${dargs[@]}" \
${DOCKER_REPOSITORY}:${GEOSX_TPL_TAG} \
${BUILD_DIR_MOUNT_POINT}/scripts/ci_build_and_test_in_container.sh \
${BUILD_AND_TEST_ARGS} ${SCCACHE_CLI};
