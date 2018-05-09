FROM ubuntu:16.04

LABEL maintainer="Randolph Settgast <settgast1@llnl.gov>"

ADD generic-setup.sh /root/generic-setup.sh

RUN /root/generic-setup.sh

USER geosx
WORKDIR /home/geosx
