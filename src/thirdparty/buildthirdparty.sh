#!/bin/bash

BUILD_DIR=$1
INSTALL_DIR=$2
C_COMPILER=$3
CXX_COMPILER=$4


###############################################################
#
# CHAI
#
###############################################################
echo "************** Starting CHAI **************"
CHAI_BUILD_DIR=$BUILD_DIR/chai
CHAI_INSTALL_DIR=$INSTALL_DIR/chai

echo $CHAI_INSTALL_DIR

CHAI_BUILD_TYPE=cpu-no-rm

# Setup build and install directories

pushd $CHAI_BUILD_DIR/src
    echo "**** Building CHAI ****"
    make $CHAI_BUILD_TYPE CPP=$CXX_COMPILER $CHAI_ARGS CALIPER_DIR=$CALIPER_INSTALL C_COMPILER=$C_COMPILER
    if [ $? -ne 0 ]; then
        echo Error: CHAI make failed
        exit 1
    fi

    echo "**** Installing CHAI ****"
    mkdir -p $CHAI_INSTALL_DIR/lib && cp *.a $CHAI_INSTALL_DIR/lib
    if [ $? -ne 0 ]; then
        echo Error: CHAI install libraries failed
        exit 1
    fi

    mkdir -p $CHAI_INSTALL_DIR/include && cp *.h* $CHAI_INSTALL_DIR/include
    mkdir -p $CHAI_INSTALL_DIR/include/yallutil && cp yallutil/optional_locking_mechanism.cpp $CHAI_INSTALL_DIR/include/yallutil
    if [ $? -ne 0 ]; then
        echo Error: CHAI install headers failed
        exit 1
    fi
popd
echo "************** Finished CHAI **************"

###############################################################
#
# RAJA
#
###############################################################
echo "************** Starting RAJA **************"
RAJA_BUILD_DIR=$BUILD_DIR/raja
RAJA_INSTALL_DIR=$INSTALL_DIR/raja
RAJA_SOURCE_DIR=$BUILD_DIR/../../src/thirdparty/raja

echo RAJA_BUILD_DIR = $RAJA_BUILD_DIR
echo RAJA_INSTALL_DIR = $RAJA_INSTALL_DIR
echo RAJA_SOURCE_DIR = $RAJA_SOURCE_DIR

CUDA_ENABLED=OFF

# Setup build and install directories
rm -rf $RAJA_BUILD_DIR
rm -rf $RAJA_INSTALL_DIR
mkdir -p $RAJA_BUILD_DIR

pushd $RAJA_BUILD_DIR
    echo "**** Configuring RAJA ****"
    cmake $RAJA_SOURCE_DIR \
    -DCMAKE_C_COMPILER=$C_COMPILER \
    -DCMAKE_CXX_COMPILER=$CXX_COMPILER \
    -DRAJA_ENABLE_CUDA=$CUDA_ENABLED \
    -DRAJA_ENABLE_TESTS=OFF \
    -DCMAKE_INSTALL_PREFIX=$RAJA_INSTALL_DIR
    if [ $? -ne 0 ]; then
        echo Error: RAJA CMake failed
        exit 1
    fi

    echo "**** Building RAJA ****"
    make VERBOSE=1 all
    if [ $? -ne 0 ]; then
        echo Error: RAJA 'make all' failed
        exit 1
    fi

    echo "**** Installing RAJA ****"
    make install
    if [ $? -ne 0 ]; then
        echo Error: RAJA 'make install' failed
        exit 1
    fi
popd
echo "************** Finished RAJA **************"
