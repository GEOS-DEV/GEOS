#

PETSC_DIR=/Users/hannahmorgan/work/geos/thirdPartyLibs/build-default-release/petsc/src/petsc
PETSC_ARCH=arch-darwin-c-debug

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
include ${PETSC_DIR}/lib/petsc/conf/test

CLINKER = g++
CFLAGS =
CPPFLAGS = 
LIBFILES =
TARGET = testPetsc
OBJFILES = testPetsc.o
CLEANFILES = $(TARGET)

all: $(TARGET)

$(TARGET): $(OBJFILES)
	${CLINKER} -o $(TARGET) $(OBJFILES) ${PETSC_LIB}
