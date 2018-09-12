FROM geosx/compiler:ubuntu18

LABEL maintainer="Randolph Settgast <settgast1@llnl.gov>"

ENV CC=clang-6.0
ENV CXX=clang++-6.0
ENV OMPI_CC=clang-6.0
ENV OMPI_CXX=clang++-6.0
ENV MPICC=mpicc
ENV MPICXX=mpicxx

ADD setup.sh /home/geosx/setup.sh
ADD buildtpl.sh /home/geosx/buildtpl.sh

RUN /home/geosx/setup.sh
RUN /home/geosx/buildtpl.sh

USER geosx
WORKDIR /home/geosx