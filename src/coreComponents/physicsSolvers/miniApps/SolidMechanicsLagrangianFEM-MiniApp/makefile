#-----[BuildType]------
#BuildType ?= CUDA
BuildType ?= CPU

HOSTNAME := $(shell hostname)

#OPTIONS
#OPT ?= -DUSE_GEOSX_ARRAY
#OPT ?= -DUSE_RAJA_VIEW

#-----[RAJA and CUDA directories]----
ifeq ($(BuildType), CUDA)

RAJA_DIR ?= /g/g17/vargas45/Git-Repo/CLEAN_MINI_GEOSX/RAJA/build-gpu
CHAI_DIR ?= /g/g17/vargas45/Git-Repo/CLEAN_MINI_GEOSX/gpu-chai/build
UMPIRE_DIR ?= /g/g17/vargas45/Git-Repo/CLEAN_MINI_GEOSX/gpu-chai/src/tpl/umpire/build
CUDA_DIR ?=/usr/tce/packages/cuda/cuda-9.2.148
OPT += -DUSE_GPU
CXX_UTIL_DIR ?= /g/g17/vargas45/Git-Repo/CLEAN_MINI_GEOSX/GEOSX/src/coreComponents/cxx-utilities
else

RAJA_DIR ?= /g/g17/vargas45/Git-Repo/CLEAN_MINI_GEOSX/RAJA/build-cpu
CHAI_DIR ?= /g/g17/vargas45/Git-Repo/CLEAN_MINI_GEOSX/cpu-chai/build
UMPIRE_DIR ?= /g/g17/vargas45/Git-Repo/CLEAN_MINI_GEOSX/cpu-chai/src/tpl/umpire/build
CXX_UTIL_DIR ?= /g/g17/vargas45/Git-Repo/CLEAN_MINI_GEOSX/cpu-cxx-utilities
endif



rajaInc = -I$(RAJA_DIR)/include
rajaLib = $(RAJA_DIR)/lib/libRAJA.a
chaiInc= -I$(CHAI_DIR)/include
umpireInc= -I$(UMPIRE_DIR)/include

chaiLib = $(CHAI_DIR)/lib/*.a
umpireLib = $(UMPIRE_DIR)/lib/*.a

cxxInc = -I$(CXX_UTIL_DIR)/build-blueos_3_ppc64le_ib-clang@coral-release/include
cxxInc2 = -I$(CXX_UTIL_DIR)/src/src
cxxLib = $(CXX_UTIL_DIR)/build-blueos_3_ppc64le_ib-clang@coral-release/lib/libcxx-utilities.a
cudaLib = -Wl,-rpath -Wl,$(CUDA_DIR)/lib64 -L$(CUDA_DIR)/lib64 -lcuda -lcudart -lcudadevrt -lnvToolsExt
#===================================

#---[Host compiler]-----
ifeq ($(findstring quartz,$(HOSTNAME)),quartz)
host-compiler = /usr/tce/packages/gcc/gcc-7.1.0/bin/g++
host-compilerFlags = '-std=c++14 -fopenmp -O3 -ffast-math -finline-functions -unroll-agressive -march=native -mavx'
compilerFlags = -std=c++14 -fopenmp -O3 -ffast-math -finline-functions -unroll-agressive -march=native -mavx

#host-compiler = /usr/tce/packages/clang/clang-4.0.0/bin/clang++
host-compilerFlags = '-std=c++14 -fopenmp -O3 -ffast-math -finline-functions -unroll-agressive -march=native -mavx'
compilerFlags = -std=c++14 -fopenmp -O3 -ffast-math -finline-functions -unroll-agressive -march=native -mavx

else ifeq ($(findstring ray,$(HOSTNAME)),ray)
host-compiler = /usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-clang-coral-2018.08.08/bin/mpic++
#host-compiler = /usr/tce/packages/xl/xl-beta-2018.06.01/bin/xlC
#host-compilerFlags = '-std=c++14 -fopenmp -O3 -qinline'
#compilerFlags = -std=c++14 -fopenmp -O3 -qinline

#Clang CORAL
host-compiler = /usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-clang-coral-2018.08.08/bin/mpic++
host-compilerFlags = '-std=c++14 -fopenmp -O3 -ffast-math -finline-functions -mpower8-vector -ffast-math'
compilerFlags = -std=c++14 -fopenmp -O3 -ffast-math -finline-functions -mpower8-vector -ffast-math


else ifeq ($(findstring cori,$(HOSTNAME)),cori)
host-compiler = cc
host-compilerFlags = '-O3 -std=c++14  -qopenmp  -xMIC-AVX512'

else
host-compiler = g++-6
host-compilerFlags = '-O3 -g -std=c++14 -m64 -fopenmp'

endif
paths = -I ./$(iPath)
paths += $(rajaInc) $(umpireInc) $(chaiInc) $(cxxInc) $(cxxInc2)
linker = $(host-compiler)
#======================

#----[device compiler]----
device-compiler=nvcc
device-flags   = -g -std=c++14 -Xptxas=-v -lineinfo --expt-extended-lambda --restrict -O3 --expt-relaxed-constexpr -use_fast_math
device-flags  += -ccbin=$(linker) -Xcompiler $(host-compilerFlags) -x=cu -arch=sm_60
#======================

#----[Cuda - Compilation]---------
ifeq ($(BuildType),CUDA)

all: ArrayOfObjects ObjectOfArrays

#Note order
ArrayOfObjects:
	@echo Compiling for CUDA - start
	$(device-compiler) $(device-flags) $(OPT) $(paths) -g -c -o ArrayOfObjects.o ArrayOfObjects.cpp
	$(linker) -o ArrayOfObjects ArrayOfObjects.o $(cudaLib) -fopenmp $(rajaLib) $(chaiLib) $(chaiLibLoc) $(cxxLib) $(umpireLib) 
	@echo Compiling for CUDA - end

ObjectOfArrays:
	@echo Compiling for CUDA - start
	$(device-compiler) $(device-flags) $(OPT) $(paths) -g -c -o ObjectOfArrays.o ObjectOfArrays.cpp
	$(linker) -o ObjectOfArrays ObjectOfArrays.o $(cudaLib) -fopenmp $(rajaLib) $(chaiLib) $(chaiLibLoc) $(cxxLib) $(umpireLib) 
	@echo Compiling for CUDA - end

else
#----[CPU - Compilation]---------
all: ArrayOfObjects ObjectOfArrays

ArrayOfObjects:
	@echo Compiling for CPU - start
	$(host-compiler) $(compilerFlags) $(OPT) $(paths) -g -o ArrayOfObjects ArrayOfObjects.cpp $(rajaLib) $(chaiLib) $(chaiLibLoc) $(cxxLib) $(umpireLib)
	@echo Compiling for CPU - end
	@echo on $(HOSTNAME)

ObjectOfArrays: 
	@echo Compiling for CPU - start
	$(host-compiler) $(compilerFlags) $(OPT) $(paths) -g -o ObjectOfArrays ObjectOfArrays.cpp $(rajaLib) $(chaiLib) $(chaiLibLoc) $(cxxLib) $(umpireLib)
	@echo Compiling for CPU - end
	@echo on $(HOSTNAME)
endif
#======================



#-----[Clean up]-------
clean:
	rm -rf ArrayOfObjects ArrayOfObjects.o ObjectOfArrays ObjectOfArrays.o
