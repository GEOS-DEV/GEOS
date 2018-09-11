FROM geosx/compiler:ubuntu18

LABEL maintainer="Randolph Settgast <settgast1@llnl.gov>"

ENV gccver=8
ENV CC=gcc-8
ENV CXX=g++-8
ENV OMPI_CC=gcc-8
ENV OMPI_CXX=g++-8
ENV MPICC=mpicc
ENV MPICXX=mpicxx

ADD setup.sh /home/geosx/setup.sh
ADD buildtpl.sh /home/geosx/buildtpl.sh

RUN /home/geosx/setup.sh
RUN /home/geosx/buildtpl.sh

USER geosx
WORKDIR /home/geosx