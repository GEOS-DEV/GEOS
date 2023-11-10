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

args=$(getopt -a -o h --long docker-repository:,docker-tag:,exchange:,help -- "$@")
if [[ $? -gt 0 ]]; then
  echo "Error after getopt"
  usage
fi

# DATA_EXCHANGE=/tmp/exchange  # TODO get from outside

eval set -- ${args}
while :
do
  case $1 in
    --docker-repository) DOCKER_REPOSITORY=$2; shift 2;;
    --docker-tag)        DOCKER_TAG=$2;        shift 2;;
    --exchange)          DATA_EXCHANGE=$2;     shift 2;;
    -h | --help)         usage;                shift;;
    # -- means the end of the arguments; drop this, and break out of the while loop
    --) shift; break;;
    *) >&2 echo Unsupported option: $1
       usage;;
  esac
done

ADDITIONAL_ARGS=$@
echo "Additional arguments '${ADDITIONAL_ARGS}' will be transfered to the final build."

if [[ -z "${DATA_EXCHANGE}" ]]; then
  echo "Variable DATA_EXCHANGE is either empty or not defined. Please define it using '--exchange'."
  exit 1
fi
DATA_EXCHANGE_MOUNT_POINT=/tmp/exchange-in-docker

# We need to know where the code folder is mounted inside the container so we can run the script at the proper location!
# Since this information is repeated twice, we use a variable.
GITHUB_WORKSPACE_MOUNT_POINT=/tmp/geos

# We need to keep track of the building container (hence the `CONTAINER_NAME`)
# so we can extract the data from it later (if needed). Another solution would have been to use a mount point,
# but that would not have solved the problem for the TPLs (we would require extra action to copy them to the mount point).
# CONTAINER_NAME=geos_build

# Now we can build GEOS.
docker run \
  --volume=${GITHUB_WORKSPACE}:${GITHUB_WORKSPACE_MOUNT_POINT} \
  --volume=${DATA_EXCHANGE}:${DATA_EXCHANGE_MOUNT_POINT} \
  --cap-add=SYS_PTRACE \
  -e ENABLE_HYPRE=${ENABLE_HYPRE:-OFF} \
  -e ENABLE_HYPRE_DEVICE=${ENABLE_HYPRE_DEVICE:-CPU} \
  -e ENABLE_TRILINOS=${ENABLE_TRILINOS:-ON} \
  ${DOCKER_REPOSITORY}:${DOCKER_TAG} \
  ${GITHUB_WORKSPACE_MOUNT_POINT}/scripts/ci_build_and_test_in_container_args.sh \
    --repository ${GITHUB_WORKSPACE_MOUNT_POINT} \
    --exchange ${DATA_EXCHANGE_MOUNT_POINT} \
    ${ADDITIONAL_ARGS}

  # --name=${CONTAINER_NAME} \
    # --install-dir ${GEOSX_DIR} \
