FROM docker.io/geosx/ubuntu20.04-gcc10:false-794

RUN apt-get install -y \
    ca-certificates \
    curl \
    git-lfs \
    vim \
    ninja-build \
# for paraview headless
    libcurl4 \
    libglapi-mesa

ARG PARAVIEW_VERSION=5.10.1
RUN curl -fsSL "https://www.paraview.org/paraview-downloads/download.php?submit=Download&version=v${PARAVIEW_VERSION%.[0-9]*}&type=binary&os=Linux&downloadFile=ParaView-${PARAVIEW_VERSION}-osmesa-MPI-Linux-Python3.9-x86_64.tar.gz" | tar --directory=/usr/local --strip-components=1 -xzf -
