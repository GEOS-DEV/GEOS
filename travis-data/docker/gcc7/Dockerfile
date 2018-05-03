FROM geosx/compiler:ubuntu-clang-base

LABEL maintainer="Randolph Settgast <settgast1@llnl.gov>"

ENV gccver=7

ADD setup.sh /home/geosx/setup.sh

RUN /home/geosx/setup.sh

USER geosx
WORKDIR /home/geosx