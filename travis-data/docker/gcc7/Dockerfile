FROM geosx/compiler:ubuntu-clang-base

LABEL maintainer="Randolph Settgast <settgast1@llnl.gov>"

ENV gccver=7
ENV CC=gcc-7
ENV CXX=g++-7
ENV OMPI_CC=gcc-7
ENV OMPI_CXX=g++-7
ENV MPICC=mpicc
ENV MPICXX=mpicxx

ADD setup.sh /home/geosx/setup.sh

RUN /home/geosx/setup.sh

USER geosx
WORKDIR /home/geosx