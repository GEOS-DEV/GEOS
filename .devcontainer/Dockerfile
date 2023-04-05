# Replace `DEFINE_ME` by the value of `GEOSX_TPL_TAG` in the `.travis.yml` file.
FROM docker.io/geosx/ubuntu20.04-gcc10:DEFINE_ME

RUN apt-get update
RUN apt-get install -y \
    ca-certificates \
    curl \
    git-lfs \
    vim \
    ninja-build \
    gdb \
# for paraview headless
    libcurl4 \
    libglapi-mesa
