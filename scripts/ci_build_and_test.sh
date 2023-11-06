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

SCCACHE_CLI="--no-use-sccache"
SCCACHE_VOLUME_MOUNT=""
if [ ${USE_SCCACHE} = true ]; then
  echo "Creating the configuration file for sccache..."
  mkdir -p /opt/gcs
  cp -rp ${GOOGLE_GHA_CREDS_PATH} /opt/gcs/credentials.json
  SCCACHE_CLI=""
  SCCACHE_VOLUME_MOUNT="--volume=/opt/gcs:/opt/gcs"
fi

# We need to keep track of the building container (hence the `CONTAINER_NAME`)
# so we can extract the data from it later (if needed). Another solution would have been to use a mount point,
# but that would not have solved the problem for the TPLs (we would require extra action to copy them to the mount point).
CONTAINER_NAME=geosx_build
# Now we can build GEOSX.
docker run \
  --name=${CONTAINER_NAME} \
  --volume=${BUILD_DIR}:${BUILD_DIR_MOUNT_POINT} ${SCCACHE_VOLUME_MOUNT} \
  --cap-add=ALL \
  -e ENABLE_HYPRE=${ENABLE_HYPRE:-OFF} \
  -e ENABLE_HYPRE_DEVICE=${ENABLE_HYPRE_DEVICE:-CPU} \
  -e ENABLE_TRILINOS=${ENABLE_TRILINOS:-ON} \
  ${DOCKER_REPOSITORY}:${GEOSX_TPL_TAG} \
  ${BUILD_DIR_MOUNT_POINT}/scripts/ci_build_and_test_in_container_args.sh \
    --cmake-build-type ${CMAKE_BUILD_TYPE} \
    --install-dir ${GEOSX_DIR} \
    --host-config ${HOST_CONFIG:-host-configs/environment.cmake} ${BUILD_AND_TEST_ARGS} ${SCCACHE_CLI};
