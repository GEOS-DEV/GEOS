#!/usr/bin/env zsh

docker build --tag geosx/compiler:ubuntu18 ubuntu18
docker build --tag geosx/compiler:gcc7-ubuntu18 gcc7-ubuntu18
docker build --tag geosx/compiler:gcc8-ubuntu18 gcc8-ubuntu18
docker build --tag geosx/compiler:clang5-ubuntu18 clang5-ubuntu18
docker build --tag geosx/compiler:clang6-ubuntu18 clang6-ubuntu18
docker build --tag geosx/compiler:centos7.5 centos7.5
docker build --tag geosx/compiler:clang5-centos clang5-centos


docker push geosx/compiler:ubuntu18
docker push geosx/compiler:gcc7-ubuntu18
docker push geosx/compiler:gcc8-ubuntu18
docker push geosx/compiler:clang5-ubuntu18
docker push geosx/compiler:clang6-ubuntu18
docker push geosx/compiler:centos7.5
docker push geosx/compiler:clang5-centos
