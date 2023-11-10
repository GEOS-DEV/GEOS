#!/bin/bash

printenv

echo "Running CLI $0 $@"

function usage () {
>&2 cat << EOF
Usage: $0
  [ --docker-repository ... ]
  [ --docker-tag ... ]
  [ -h | --help ]
EOF
exit 1
}

args=$(getopt -a -o h --long docker-repository:,docker-tag:,help -- "$@")
if [[ $? -gt 0 ]]; then
  echo "Error after getopt"
  usage
fi

# USE_SCCACHE=true
eval set -- ${args}
while :
do
  case $1 in
    --docker-repository) DOCKER_REPOSITORY=$2; shift 2;;
    --docker-tag)        DOCKER_TAG=$2;        shift 2;;
    -h | --help)         usage;                shift;;
    # -- means the end of the arguments; drop this, and break out of the while loop
    --) shift; break;;
    *) >&2 echo Unsupported option: $1
       usage;;
  esac
done

ADDITIONAL_ARGS=$@
echo "Additional arguments '${ADDITIONAL_ARGS}' will be transfered to the final build."

# The linux build relies on the two variables DOCKER_REPOSITORY and DOCKER_TAG to define the proper image version.
DOCKER_IMAGE=${DOCKER_REPOSITORY}:${DOCKER_TAG}
# We extract the location of the GEOSX_TPL from the container...
GEOSX_TPL_DIR=$(docker run --rm ${DOCKER_IMAGE} /bin/bash -c 'echo ${GEOSX_TPL_DIR}' | tail -1)
# ... so we can install GEOSX alongside. This is assumed for bundling the binaries, so consider modifying with care.
GEOSX_DIR=${GEOSX_TPL_DIR}/../GEOSX-${GITHUB_SHA:0:7}
# We need to know where the code folder is mounted inside the container so we can run the script at the proper location!
# Since this information is repeated twice, we use a variable.
GITHUB_WORKSPACE_MOUNT_POINT=/tmp/geos

# SCCACHE_VOLUME_MOUNT=""
# if [ ${USE_SCCACHE} = true ]; then
#   echo "Creating the configuration file for sccache..."
#   echo "File path is ${GOOGLE_GHA_CREDS_PATH}"
#   echo "GITHUB_WORKSPACE is ${GITHUB_WORKSPACE}"
#   ls /home/runner/work/GEOS/GEOS  # -> it's the root of the repository
#   mkdir -p /opt/gcs
#   cp -rp ${GOOGLE_GHA_CREDS_PATH} /opt/gcs/credentials.json
#   SCCACHE_VOLUME_MOUNT="--volume=/opt/gcs:/opt/gcs"
# fi

# We need to keep track of the building container (hence the `CONTAINER_NAME`)
# so we can extract the data from it later (if needed). Another solution would have been to use a mount point,
# but that would not have solved the problem for the TPLs (we would require extra action to copy them to the mount point).
CONTAINER_NAME=geos_build
# Now we can build GEOS.
docker run \
  --name=${CONTAINER_NAME} \
  --volume=${GITHUB_WORKSPACE}:${GITHUB_WORKSPACE_MOUNT_POINT} \
  --cap-add=SYS_PTRACE \
  -e ENABLE_HYPRE=${ENABLE_HYPRE:-OFF} \
  -e ENABLE_HYPRE_DEVICE=${ENABLE_HYPRE_DEVICE:-CPU} \
  -e ENABLE_TRILINOS=${ENABLE_TRILINOS:-ON} \
  ${DOCKER_IMAGE} \
  ${GITHUB_WORKSPACE_MOUNT_POINT}/scripts/ci_build_and_test_in_container_args.sh \
    --install-dir ${GEOSX_DIR} \
    ${ADDITIONAL_ARGS}
