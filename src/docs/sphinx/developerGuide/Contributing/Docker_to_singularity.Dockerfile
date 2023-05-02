# The docker base image
# Work on TPL docker image pull-able from https://hub.docker.com/u/geosx
#
# Main Layout (eventually):
#               /opt/GEOSX/
#                   |- GEOSX_TPL-xxx/
#                   |- GEOSX-version/
#                   |- inputFiles/
#
# install GEOSX in /opt/GEOSX with HYPRE ON (hardcoded)
# build-arg:
#   DOCKER_ROOT_IMAGE has to be chosen among different geosx tpl images on dockerhub
#   GEOSX_TPL_DIR a.k.a can then be get via
#           > export GEOSX_TPL_DIR=$(docker run --rm ${DOCKER_ROOT_IMAGE} /bin/bash -c 'echo ${GEOSX_TPL_DIR}')
#   HOST_CONFIG has to be set to e.g. environment.cmake

ARG DOCKER_ROOT_IMAGE
ARG HOST_CONFIG=environment.cmake
ARG GEOSX_TPL_DIR

FROM ${DOCKER_ROOT_IMAGE} as tpl_toolchain_intersect_geosx_toolchain

RUN apt-get -y update
RUN yes "y"|  apt-get install --fix-missing --reinstall ca-certificates
RUN apt-get install -y \
    make \
    python3 \
    zlib1g-dev\
    curl

ARG CMAKE_VERSION=3.20.3
RUN curl -sL https://cmake.org/files/v${CMAKE_VERSION%.[0-9]*}/cmake-${CMAKE_VERSION}-linux-x86_64.tar.gz | tar --directory=/usr/local --strip-components=1 -xzf -

FROM tpl_toolchain_intersect_geosx_toolchain AS geosx_toolchain

ARG HOST_CONFIG

RUN apt-get -y install \
    texlive \
    graphviz \
    libxml2 \
    git \
	git-lfs

ARG GEOSX_SRC_DIR=/tmp/src
ARG GEOSX_BUILD_DIR=${GEOSX_SRC_DIR}/build

RUN git clone https://github.com/GEOS-DEV/GEOS.git ${GEOSX_SRC_DIR}
WORKDIR ${GEOSX_SRC_DIR}
RUN git pull
RUN git submodule update --init src/cmake/blt
RUN	git submodule update --init src/coreComponents/LvArray
RUN	git submodule update --init src/coreComponents/fileIO/coupling/hdf5_interface
RUN	git submodule update --init src/coreComponents/constitutive/PVTPackage

RUN git submodule update

ENV ENABLE_HYPRE=ON
ENV ENABLE_CUDA=OFF
RUN python3 scripts/config-build.py --hostconfig=./host-configs/${HOST_CONFIG} --buildtype=Release --buildpath=${GEOSX_BUILD_DIR} --installpath=/opt/GEOSX/GEOSX-version
WORKDIR ${GEOSX_BUILD_DIR}
RUN make -j 8 && make install

RUN mv ${GEOSX_SRC_DIR}/inputFiles/ /opt/GEOSX/
RUN rm -rf /tmp/src

CMD ["/opt/GEOSX/GEOSX-version/bin/geosx"]