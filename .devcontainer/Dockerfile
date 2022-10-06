FROM docker.io/geosx/ubuntu20.04-gcc10:203-819

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
