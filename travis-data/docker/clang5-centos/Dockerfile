FROM geosx/compiler:centos7.5

LABEL maintainer="Randolph Settgast <settgast1@llnl.gov>"

ENV CC=/opt/rh/llvm-toolset-7/root/usr/bin/clang
ENV CXX=/opt/rh/llvm-toolset-7/root/usr/bin/clang++
ENV OMPI_CC=/opt/rh/llvm-toolset-7/root/usr/bin/clang
ENV OMPI_CXX=/opt/rh/llvm-toolset-7/root/usr/bin/clang++
ENV MPICC=/usr/lib64/openmpi/bin/mpicc
ENV MPICXX=/usr/lib64/openmpi/bin/mpicxx

ADD setup.sh /home/geosx/setup.sh
ADD buildtpl.sh /home/geosx/buildtpl.sh

RUN /home/geosx/setup.sh
RUN /home/geosx/buildtpl.sh

USER geosx
WORKDIR /home/geosx