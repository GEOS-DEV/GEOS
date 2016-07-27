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
-include ../../Make.local.defs.$(shell hostname|sed 's/[0-9]//g')
-include ../../Make.local.defs.$(shell whoami).$(shell hostname|sed 's/[0-9]//g')

###########################################################################

PROGRAM = bifroesttest
#CPP_FLAGS += -O3 -DGPAC_MPI
CPP_FLAGS += -g -DGPAC_MPI
INC = -I../../

CSRC  = $(wildcard *Bifroest*.cpp)
COBJS = $(patsubst %.cpp, %.o, $(CSRC))

all:	
	make clean
	make $(COBJS) Communication.o main.o
#$(CC) $(CPP_FLAGS) $(INC) Communication.cpp -c
#$(CC) $(CPP_FLAGS) $(INC) main.cpp -c
	make $(PROGRAM)

$(PROGRAM):	$(COBJS)
	$(CC) -o $@ $(COBJS) Communication.o main.o

%.o: %.cpp
	$(CC) $(CPP_FLAGS) $(INC) -c $<
clean:
	-rm $(PROGRAM) *.o
