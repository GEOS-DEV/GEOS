#!/bin/bash

# these come from CMake
MAKE=$1
TARGET=$2

# let the user see warnings, but also write to file
$MAKE $TARGET 2> >(tee doxygen.err) >/dev/null
err=$(cat doxygen.err | wc -l)
exit $err
