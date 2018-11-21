#-----[Build-type]------
Build-type = CUDA
#Build-type  = CPU

HOSTNAME := $(shell hostname)

#-----[RAJA and CUDA directories]----
#RAJA_DIR ?= /global/homes/v/vargas45/Git-Repos/RAJA/develop/build
#RAJA_DIR ?= /g/g17/vargas45/Git-Repo/RAJA/develop/build
#RAJA_DIR ?= /g/g17/vargas45/Git-Repo/MINI_GEOSX/RAJA/build
RAJA_DIR ?= /g/g17/vargas45/Git-Repo/MINI_GEOSX/RAJA/build
CHAI_DIR ?= /g/g17/vargas45/Git-Repo/MINI_GEOSX/CHAI/build
CXX_UTIL_DIR ?= /g/g17/vargas45/Git-Repo/MINI_GEOSX/GEOSX/src/coreComponents/cxx-utilities/build-blueos_3_ppc64le_ib-clang@coral-release
CUDA_DIR ?= /usr/tce/packages/cuda/cuda-9.2.64

rajaInc = -I$(RAJA_DIR)/include
rajaLib = $(RAJA_DIR)/lib/libRAJA.a
chaiInc = -I$(CHAI_DIR)/include
chaiLib = $(CHAI_DIR)/lib/libchai.a
cxxInc = -I$(CXX_UTIL_DIR)/include
cxxInc2 = -I/g/g17/vargas45/Git-Repo/MINI_GEOSX/GEOSX/src/coreComponents/cxx-utilities/src/src
cxxLib = $(CXX_UTIL_DIR)/lib/libcxx-utilities.a
cudaLib = -Wl,-rpath -Wl,$(CUDA_DIR)/lib64 -L$(CUDA_DIR)/lib64 -lcuda -lcudart -lcudadevrt -lnvToolsExt
#===================================

#---[Host compiler]-----
ifeq ($(findstring quartz,$(HOSTNAME)),quartz)
host-compiler = /usr/tce/packages/gcc/gcc-7.1.0/bin/g++
host-compilerFlags = '-std=c++14 -fopenmp -O3 -ffast-math -unroll-agressive -march=native -finline-functions -mavx  -Wall -Wextra  -fopenmp  -Wall  -Wno-abi  -Wextra  -Wfloat-equal  -Wcast-align  -Wpointer-arith  -Wwrite-strings  -Wcast-qual  -Wswitch-default  -Wno-vla  -Wno-switch-default  -Wno-unused-parameter  -Wno-unused-variable  -Wno-unused-function  -rdynamic '
compilerFlags = -std=c++14 -fopenmp -O3 -ffast-math -unroll-agressive -march=native -finline-functions -mavx  -Wall -Wextra -fopenmp -Wall  -Wno-abi  -Wextra -Wfloat-equal  -Wcast-align  -Wpointer-arith  -Wwrite-strings  -Wcast-qual  -Wswitch-default  -Wno-vla  -Wno-switch-default  -Wno-unused-parameter  -Wno-unused-variable  -Wno-unused-function  -rdynamic 

#host-compiler = /usr/tce/packages/clang/clang-4.0.0/bin/clang++
#host-compilerFlags = '-std=c++14 -fopenmp -O3 -ffast-math -finline-functions  -Wall -Wextra  -fopenmp  -Wall  -Wno-abi  -Wextra  -Wfloat-equal  -Wcast-align  -Wpointer-arith  -Wwrite-strings  -Wcast-qual  -Wswitch-default  -Wno-vla  -Wno-switch-default  -Wno-unused-parameter  -Wno-unused-variable  -Wno-unused-function '
#compilerFlags = -std=c++14 -fopenmp -O3 -ffast-math -finline-functions  -Wall -Wextra  -fopenmp  -Wall  -Wno-abi  -Wextra  -Wfloat-equal  -Wcast-align  -Wpointer-arith  -Wwrite-strings  -Wcast-qual  -Wswitch-default  -Wno-vla  -Wno-switch-default  -Wno-unused-parameter  -Wno-unused-variable  -Wno-unused-function 


else ifeq ($(findstring ray,$(HOSTNAME)),ray)
#host-compiler = /usr/tce/packages/xl/xl-beta-2018.06.01/bin/xlC
#host-compilerFlags = '-std=c++14 -fopenmp -O3 -qinline'
#compilerFlags = -std=c++14 -fopenmp -O3 -qinline

#Clang CORAL
host-compiler = /usr/tce/packages/clang/clang-coral-2018.05.23/bin/clang++
host-compilerFlags = '-std=c++14 -fopenmp -O3 -ffast-math -finline-functions  -Wall -Wextra  -fopenmp  -Wall  -Wno-abi  -Wextra  -Wfloat-equal  -Wcast-align  -Wpointer-arith  -Wwrite-strings  -Wcast-qual  -Wswitch-default  -Wno-vla  -Wno-switch-default  -Wno-unused-parameter  -Wno-unused-variable  -Wno-unused-function -mpower8-vector -ffast-math'
compilerFlags = -std=c++14 -fopenmp -O3 -ffast-math -finline-functions  -Wall -Wextra  -fopenmp  -Wall  -Wno-abi  -Wextra  -Wfloat-equal  -Wcast-align  -Wpointer-arith  -Wwrite-strings  -Wcast-qual  -Wswitch-default  -Wno-vla  -Wno-switch-default  -Wno-unused-parameter  -Wno-unused-variable  -Wno-unused-function -mpower8-vector -ffast-math

else ifeq ($(findstring cori,$(HOSTNAME)),cori)
host-compiler = cc
host-compilerFlags = '-O3 -std=c++14  -qopenmp  -xMIC-AVX512'
compilerFlags = -O3 -std=c++14  -qopenmp  -xMIC-AVX512

else
host-compiler = g++-6
host-compilerFlags = '-O3 -g -std=c++14 -m64 -fopenmp'
compilerFlags = -O3 -g -std=c++14 -m64 -fopenmp
endif
paths = -I ./$(iPath)
paths += $(rajaInc) $(chaiInc) $(cxxInc) $(cxxInc2)
linker = $(host-compiler)
#======================

#----[device compiler]----
device-compiler=nvcc
device-flags   = -g -std=c++14 -Xptxas=-v -lineinfo --expt-extended-lambda --restrict -O3
device-flags  += -ccbin=$(linker) -Xcompiler $(host-compilerFlags) -x=cu -arch=sm_60
#======================

#----[Cuda - Compilation]---------
ifeq ($(Build-type),CUDA)

all: ArrayOfObjects

ArrayOfObjects:
	@echo Compiling for CUDA - start
	$(device-compiler) $(device-flags) $(paths) -g -c -o ArrayOfObjects.o  ArrayOfObjects.cpp
	$(linker) -o ArrayOfObjects ArrayOfObjects.o $(cudaLib) -fopenmp $(rajaLib) $(chaiLib) $(cxxLib)
	@echo Compiling for CUDA - end

#ObjectOfArrays: ObjectOfArrays
#	@echo Compiling for CUDA - start
#	$(device-compiler) $(device-flags) $(paths) -g -c -o ObjectOfArrays.o ObjectOfArrays.cpp
#	$(linker) -o ObjectOfArrays ObjectOfArrays.o $(cudaLib) -fopenmp $(rajaLib)
#	@echo Compiling for CUDA - end

else
#----[CPU - Compilation]---------
all: ArrayOfObjects

ArrayOfObjects:
	@echo Compiling for CPU - start
	$(host-compiler) $(compilerFlags) $(paths) -g -o ArrayOfObjects ArrayOfObjects.cpp
	@echo Compiling for CPU - end
	@echo on $(HOSTNAME)

#ObjectOfArrays: ObjectOfArrays
#	@echo Compiling for CPU - start
#	$(host-compiler) $(compilerFlags) $(paths) -g -o ObjectOfArrays ObjectOfArrays.cpp
#	@echo Compiling for CPU - end
#	@echo on $(HOSTNAME)
endif
#======================



#-----[Clean up]-------
clean:
	rm -rf ArrayOfObjects ArrayOfObjects.o ObjectOfArrays ObjectOfArrays.o
