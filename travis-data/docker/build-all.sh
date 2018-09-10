#!/usr/bin/env zsh

if (( $# == 0 )) ; then
  dockerfiles=(*ubuntu18*/Dockerfile)
else
  dockerfiles=($@)
fi

all_images="${compiler_images} ubuntu18"

function build-tag () {
  echo building $2 and tagging with $1
  docker build --tag $1 $2
}

build-tag geosx/compiler:ubuntu18 ubuntu18
docker push geosx/compiler:ubuntu18

echo ${dockerfiles}
for df in ${dockerfiles} ; do
  imgname=${df:h}
  imgpath=geosx/compiler:$imgname
  echo $imgname
  echo $imgpath
  [[ ${imgname} == 'ubuntu18' ]] && continue
  build-tag $imgpath $imgname

  echo pushing $imgpath
  docker push $imgpath
done

#for img in $all_images ; do
  # docker push geosx/compiler:$img
#done 
