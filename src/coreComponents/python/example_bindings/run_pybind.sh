!/bin/bash

export PYTHONHOME=/g/g17/sherman/anaconda2
export LD_LIBRARY_PATH=/g/g17/sherman/anaconda2/lib

# Arguments are path, module, x, y:
pybind /g/g17/sherman/GEOS/parse_geos_xml/bindings test_mod 1.2 3.4
