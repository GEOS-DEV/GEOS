###########################################################################
#
# Override default makefile settings
#
## Create or edit one of the following files:
##   Make.local.defs
##   Make.local.defs.(username) [e.g. Make.local.defs.benliu]
##   Make.local.defs.(host) [e.g. Make.local.defs.mcr]
##   Make.local.defs.(username).(host) [e.g. Make.local.defs.benliu.mcr]
## in the build directory to change the default compilation behavior
## `Make.local.defs.defaults' is provided as a template.

-include ../../Make.local.defs
-include ../../Make.local.defs.$(shell whoami)
#-include ../Make.local.defs.$(shell hostname|sed 's/[0-9]//g')
#-include ../Make.local.defs.$(shell whoami).$(shell hostname|sed 's/[0-9]//g')

###########################################################################


MAIN_OBJ ?=\
	main.o

3D_OBJ ?=\
	../Interpolation/VolumeKernel.o \
	FractalVolume.o \
	main_volume.o

LIB_OBJ ?=\
	FractalSurfaceLibrary.o

COMMON_OBJ ?=\
	FractalBaseT.o \
	StatisticalDistributionBaseT.o

2D_OBJ ?=\
	../Interpolation/SurfaceKernel.o \
	FractalSurface.o

OBJS = $(COMMON_OBJ) $(LIB_OBJ) $(MAIN_OBJ) $(3D_OBJ) $(2D_OBJ)

LIB_OBJS = $(COMMON_OBJ) $(2D_OBJ) $(LIB_OBJ)
MAIN_OBJS = $(COMMON_OBJ) $(2D_OBJ) $(MAIN_OBJ)
VMAIN_OBJS = $(COMMON_OBJ) $(3D_OBJ)

ifeq ($(HAVE_MPI), 1)
	CPPC=$(CPP)
else
	CPPC=$(SERIAL_CPP)
endif

CPP_FLAGS+=-I. -I../ -I../Interpolation/

ifeq ($(DEBUG), 1)
	CPP_FLAGS += -DRANGE_CHECKING=1 
else
	CPP_FLAGS += -DRANGE_CHECKING=0
endif


ifeq ($(COMPILER), gcc)
	CPP_FLAGS +=  -Wall -fmessage-length=0
	ifeq ($(BUILD_TYPE), debug)
		CPP_FLAGS +=
	else
		CPP_FLAGS +=
	endif
endif
ifeq ($(COMPILER), icc)
	ifeq ($(BUILD_TYPE), debug)
		CPP_FLAGS +=
	else
		CPP_FLAGS +=
	endif
endif

CPP_FLAGS += -ICommon

GEOS.x: $(OBJS)
	ar -cvq libfractalsurface.a $(LIB_OBJS)
	$(CPPC) $(MAIN_OBJS) -o apgen && echo && echo SUCCESS;
	$(CPPC) $(VMAIN_OBJS) -o volgen && echo && echo SUCCESS;

Make.GEOS.config: $(OBJS) Makefile
	echo GEOS_OBJS=$(OBJS) > $@
	echo GEOS_LD_FLAGS=$(LNK_FLAGS) $(LD_FLAGS) >> $@


# pull in dependency info for *existing* .o files
-include $(OBJS:.o=.d)

# compile and generate dependency info;
# more complicated dependency computation, so all prereqs listed
# will also become command-less, prereq-less targets
#   sed:    strip the target (everything before colon)
#   sed:    remove any continuation backslashes
#   fmt -1: list words one per line
#   sed:    strip leading spaces
#   sed:    add trailing colons

%.o: %.cpp
	$(CPPC) -c $(CPP_FLAGS) $*.cpp -o $*.o
	$(CPPC) -MM $(CPP_FLAGS) $*.cpp > $*.d
	@mv -f $*.d $*.d.tmp
	@sed -e 's|.*:|$*.o:|' < $*.d.tmp > $*.d
	@sed -e 's/.*://' -e 's/\\$$//' < $*.d.tmp | fmt -1 | \
	  sed -e 's/^ *//' -e 's/$$/:/' >> $*.d
	@rm -f $*.d.tmp

clean :
	rm -f apgen volgen $(OBJS) *.d

