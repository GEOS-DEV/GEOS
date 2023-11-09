#!/bin/bash
env

echo "Running CLI $0 $@"

function usage () {
>&2 cat << EOF
Usage: $0
  [ --docker-repository ]
  [ --docker-tag ]
  [ --use-sccache (true|false) ]
  [ -h | --help ]
EOF
exit 1
}

args=$(getopt -a -o h --long docker-repository:,docker-tag:,use-sccache:,help -- "$@")
if [[ $? -gt 0 ]]; then
  echo "Error after getopt"
  usage
fi

USE_SCCACHE=true
# CMAKE_BUILD_TYPE=""
eval set -- ${args}
while :
do
  case $1 in
    # --cmake-build-type)  CMAKE_BUILD_TYPE=$2;  shift 2;;
    --docker-repository) DOCKER_REPOSITORY=$2; shift 2;;
    --docker-tag)        DOCKER_TAG=$2;        shift 2;;
    --use-sccache)       USE_SCCACHE=$2;     shift 2;;
    -h | --help)         usage;                shift;;
    # -- means the end of the arguments; drop this, and break out of the while loop
    --) shift; break;;
    *) >&2 echo Unsupported option: $1
       usage;;
  esac
done

ADDITIONAL_ARGS=$@

# if [[ -z "${CMAKE_BUILD_TYPE}" ]]; then
#   echo "Variable \"CMAKE_BUILD_TYPE\" is undefined or empty. Define it using '--cmake-build-type'."
#   exit 1
# fi

# The linux build relies on the two variables DOCKER_REPOSITORY and DOCKER_TAG to define the proper image version.
DOCKER_IMAGE=${DOCKER_REPOSITORY}:${DOCKER_TAG}#
# We extract the location of the GEOSX_TPL from the container...
GEOSX_TPL_DIR=$(docker run --rm ${DOCKER_IMAGE} /bin/bash -c 'echo ${GEOSX_TPL_DIR}' | tail -1)
# ... so we can install GEOSX alongside. This is assumed for bundling the binaries, so consider modifying with care.
GEOSX_DIR=${GEOSX_TPL_DIR}/../GEOSX-${GITHUB_SHA:0:7}
# We need to know where the code folder is mounted inside the container so we can run the script at the proper location!
# Since this information is repeated twice, we use a variable.
GITHUB_WORKSPACE_MOUNT_POINT=/tmp/geos

SCCACHE_VOLUME_MOUNT=""
if [ ${USE_SCCACHE} = true ]; then
  echo "Creating the configuration file for sccache..."
  mkdir -p /opt/gcs
  cp -rp ${GOOGLE_GHA_CREDS_PATH} /opt/gcs/credentials.json
  SCCACHE_VOLUME_MOUNT="--volume=/opt/gcs:/opt/gcs"
fi

# We need to keep track of the building container (hence the `CONTAINER_NAME`)
# so we can extract the data from it later (if needed). Another solution would have been to use a mount point,
# but that would not have solved the problem for the TPLs (we would require extra action to copy them to the mount point).
CONTAINER_NAME=geos_build
# Now we can build GEOS.
echo "docker run \
  --name=${CONTAINER_NAME} \
  --volume=${GITHUB_WORKSPACE}:${GITHUB_WORKSPACE_MOUNT_POINT} ${SCCACHE_VOLUME_MOUNT} \
  --cap-add=ALL \
  -e ENABLE_HYPRE=${ENABLE_HYPRE:-OFF} \
  -e ENABLE_HYPRE_DEVICE=${ENABLE_HYPRE_DEVICE:-CPU} \
  -e ENABLE_TRILINOS=${ENABLE_TRILINOS:-ON} \
  ${DOCKER_IMAGE} \
  ${GITHUB_WORKSPACE_MOUNT_POINT}/scripts/ci_build_and_test_in_container_args.sh \
    --install-dir ${GEOSX_DIR} \
    --host-config ${HOST_CONFIG:-host-configs/environment.cmake} \
    --use-sccache ${USE_SCCACHE} ${ADDITIONAL_ARGS};"

docker run \
  --name=${CONTAINER_NAME} \
  --volume=${GITHUB_WORKSPACE}:${GITHUB_WORKSPACE_MOUNT_POINT} ${SCCACHE_VOLUME_MOUNT} \
  --cap-add=ALL \
  -e ENABLE_HYPRE=${ENABLE_HYPRE:-OFF} \
  -e ENABLE_HYPRE_DEVICE=${ENABLE_HYPRE_DEVICE:-CPU} \
  -e ENABLE_TRILINOS=${ENABLE_TRILINOS:-ON} \
  ${DOCKER_IMAGE} \
  ${GITHUB_WORKSPACE_MOUNT_POINT}/scripts/ci_build_and_test_in_container_args.sh \
    --install-dir ${GEOSX_DIR} \
    --host-config ${HOST_CONFIG:-host-configs/environment.cmake} \
    --use-sccache ${USE_SCCACHE} ${ADDITIONAL_ARGS};



    # --cmake-build-type ${CMAKE_BUILD_TYPE} \
