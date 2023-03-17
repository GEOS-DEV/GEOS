FROM docker.io/geosx/ubuntu20.04-gcc10:203-915

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

WORKDIR /workspaces/GEOSX
RUN ln -s .vscode-codespaces .vscode
