#!/bin/bash
env

git submodule update --init --recursive src/cmake/blt
git submodule update --init --recursive src/coreComponents/LvArray
git submodule update --init --recursive src/coreComponents/constitutive/PVTPackage
git submodule update --init src/coreComponents/mesh/PAMELA
git submodule update --init --recursive src/coreComponents/fileIO/coupling/hdf5_interface
# The linux build relies on two environment variables DOCKER_REPOSITORY and GEOSX_TPL_TAG to define the TPL version.
# And another CMAKE_BUILD_TYPE to define the build type we want for GEOSX.
# Optional BUILD_AND_TEST_ARGS to pass arguments to build_test_helper.sh script.
#
# We extract the location of the GEOSX_TPL from the container...
GEOSX_TPL_DIR=$(docker run --rm ${DOCKER_REPOSITORY}:${GEOSX_TPL_TAG} /bin/bash -c 'echo ${GEOSX_TPL_DIR}')
# ... so we can install GEOSX alongside. This is assumed for bundling the binaries, so consider modifying with care.
GEOSX_DIR=${GEOSX_TPL_DIR}/../GEOSX-INSTALL
# We need to get the build directory, which is different between Travis and Azure Pipelines.
BUILD_DIR=${TRAVIS_BUILD_DIR:-$BUILD_SOURCESDIRECTORY}
# We need to know where the code folder is mounted inside the container so we can run the script at the proper location!
# Since this information is repeated twice, we use a variable.
BUILD_DIR_MOUNT_POINT=/tmp/GEOSX
# We need to keep track of the building container (hence the `CONTAINER_NAME`)
# so we can extract the data from it later (if needed). Another solution would have been to use a mount point,
# but that would not have solved the problem for the TPLs (we would require extra action to copy them to the mount point).
CONTAINER_NAME=geosx_build
# Now we can build GEOSX.
while sleep 5m; do echo "... still building ..."; done & 
docker run \
--name=${CONTAINER_NAME} \
--volume=${BUILD_SOURCESDIRECTORY}:${BUILD_DIR_MOUNT_POINT} \
--cap-add=ALL \
-e HOST_CONFIG=${HOST_CONFIG:-host-configs/environment.cmake} \
-e CMAKE_BUILD_TYPE \
-e GEOSX_DIR=${GEOSX_DIR} \
-e ENABLE_HYPRE=${ENABLE_HYPRE:-OFF} \
-e ENABLE_HYPRE_CUDA=${ENABLE_HYPRE_CUDA:-OFF} \
-e ENABLE_TRILINOS=${ENABLE_TRILINOS:-ON} \
${DOCKER_REPOSITORY}:${GEOSX_TPL_TAG} \
${BUILD_DIR_MOUNT_POINT}/scripts/build_test_helper.sh ${BUILD_AND_TEST_ARGS};

