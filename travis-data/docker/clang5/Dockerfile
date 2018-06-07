FROM geosx/compiler:ubuntu-clang-base

LABEL maintainer="Randolph Settgast <settgast1@llnl.gov>"

ENV llvmtar=clang+llvm-5.0.2-x86_64-linux-gnu-ubuntu-16.04
ENV tarext=.tar.xz
ENV CC=clang
ENV CXX=clang++
ENV OMPI_CC=clang
ENV OMPI_CXX=clang++
ENV MPICC=mpicc
ENV MPICXX=mpicxx

ADD setup.sh /home/geosx/setup.sh

RUN /home/geosx/setup.sh

USER geosx
WORKDIR /home/geosx