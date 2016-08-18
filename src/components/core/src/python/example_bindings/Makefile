ANACONDA = /usr/gapps/GEOS/python/anaconda2
CC = g++
CFLAGS = -Wall -Wextra -g
LFLAGS = 
LIBS = -L$(ANACONDA)/lib -lpython2.7
INCLUDES = -I$(ANACONDA)/include/python2.7 -I$(ANACONDA)/lib/python2.7/site-packages/numpy/core/include

SRCS = python_bind.cpp
OBJS = $(SRCS:.cpp=.o)
MAIN = pybind

#

all:    $(MAIN)
	@echo  "Finished building C++ python binding example!"
	@echo  "Writing a run script in ./run_pybind.sh"
	@echo  "!/bin/bash" > run_pybind.sh
	@echo  "" >> run_pybind.sh
	@echo  "export PYTHONHOME=$(ANACONDA)" >> run_pybind.sh
	@echo  "export LD_LIBRARY_PATH=$(ANACONDA)/lib" >> run_pybind.sh
	@echo  "" >> run_pybind.sh
	@echo  "# Arguments are path, module, x, y:" >> run_pybind.sh
	@echo  "$(MAIN) /g/g17/sherman/GEOS/parse_geos_xml/bindings test_mod 1.2 3.4" >> run_pybind.sh


$(MAIN): $(OBJS) 
	$(CC) $(CFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) $(LFLAGS) $(LIBS)

$(OBJS): %.o: %.cpp
	$(CC) $(CFLAGS) $(INCLUDES) -c $<  -o $@

clean:
	$(RM) *.o $(MAIN)

#
