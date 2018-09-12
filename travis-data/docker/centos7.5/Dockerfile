FROM centos:7.5.1804

LABEL maintainer="Randolph Settgast <settgast1@llnl.gov>"

ADD generic-setup.sh /root/generic-setup.sh

RUN /root/generic-setup.sh

USER geosx
WORKDIR /home/geosx
