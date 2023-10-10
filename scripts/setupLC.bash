#!/bin/bash

## Builds GEOS using the TPLs on all LC systems. Must be run from the top level GEOS directory.

## Trap the interupt signal and kill all children.
trap 'killall' INT

killall() {
    trap '' INT TERM     # ignore INT and TERM while shutting down
    echo "**** Shutting down. Killing chid processes ****"     # added double quotes
    kill -TERM 0         # fixed order, send TERM not INT
    wait
    echo DONE
}

mkdir toBeDeleted
mv build-* toBeDeleted
rm -rf toBeDeleted &

echo "Building all LC host-configs from GEOS."
./scripts/setupLC-helper.bash quartz clang-14 "srun -N 1 -t 90 -n 1 -A geosecp" $@ &
./scripts/setupLC-helper.bash quartz gcc-12 "srun -N 1 -t 90 -n 1 -A geosecp" $@ &
./scripts/setupLC-helper.bash lassen gcc-8-cuda-11       "lalloc 1 -qpdebug" $@ &
./scripts/setupLC-helper.bash lassen clang-13-cuda-11       "lalloc 1 -qpdebug" $@ &
./scripts/setupLC-helper.bash lassen clang-10-cuda-11       "lalloc 1 -qpdebug" $@ &

wait
echo "Complete"
