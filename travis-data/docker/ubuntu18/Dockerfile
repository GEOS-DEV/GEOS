FROM ubuntu:18.10

LABEL maintainer="Randolph Settgast <settgast1@llnl.gov>"

#RUN apt-get update && apt-get install -y vim nano zsh curl git sudo wget

ADD generic-setup.sh /root/generic-setup.sh

RUN /root/generic-setup.sh

USER geosx
WORKDIR /home/geosx
