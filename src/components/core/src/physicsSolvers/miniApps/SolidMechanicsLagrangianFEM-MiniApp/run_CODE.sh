#!/bin/bash

#Note to users:
#
# Please allocate an exclusive node prior to running this script
#

#Ray
#getNode='bsub -x -n1 -G guests -Is /bin/bash'
#run=''

##Quartz
#getNode='salloc -N 1 -n 1 -c 36 --exclusive'
#run='--cpu_bind=cores'

#CORI
#getNode='salloc -N 1 -n 1 -C knl,quad,cache -q debug -t 00:30:00 -L SCRATCH'
#run=''

#getNode
echo '-----[Array Of Objects Kernel]------' 
 ./ArrayOfObjects 10 30 
 ./ArrayOfObjects 20 30 
 ./ArrayOfObjects 40 30 
 ./ArrayOfObjects 60 30
 ./ArrayOfObjects 80 30 
 ./ArrayOfObjects 100 30
 ./ArrayOfObjects 120 30
 ./ArrayOfObjects 140 30
 ./ArrayOfObjects 160 30
 ./ArrayOfObjects 180 30
 ./ArrayOfObjects 200 30
 ./ArrayOfObjects 220 30
echo '\n'

echo '-----[Object Of Arrays Kernel]------' 
 ./ObjectOfArrays 10 30
 ./ObjectOfArrays 20 30
 ./ObjectOfArrays 40 30
 ./ObjectOfArrays 60 30
 ./ObjectOfArrays 80 30
 ./ObjectOfArrays 100 30
 ./ObjectOfArrays 120 30
 ./ObjectOfArrays 140 30
 ./ObjectOfArrays 160 30
 ./ObjectOfArrays 180 30
 ./ObjectOfArrays 200 30
 ./ObjectOfArrays 220 30
echo '\n'
