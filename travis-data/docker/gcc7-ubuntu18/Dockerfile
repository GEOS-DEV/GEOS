FROM geosx/compiler:ubuntu18

LABEL maintainer="Randolph Settgast <settgast1@llnl.gov>"

ENV gccver=7
ENV CC=gcc
ENV CXX=g++
ENV OMPI_CC=gcc
ENV OMPI_CXX=g++
ENV MPICC=mpicc
ENV MPICXX=mpicxx

ADD setup.sh /home/geosx/setup.sh
ADD buildtpl.sh /home/geosx/buildtpl.sh

RUN /home/geosx/setup.sh
RUN /home/geosx/buildtpl.sh

USER geosx
WORKDIR /home/geosx