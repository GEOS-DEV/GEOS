#!/bin/bash

srun -n1 -p pdebug python runQueries.py /p/lscratche/geosadmn/Penny/run01 "plot_* database" run01&
srun -n1 -p pdebug python runQueries.py /p/lscratche/geosadmn/Penny/run02 "plot_* database" run02&
srun -n1 -p pdebug python runQueries.py /p/lscratche/geosadmn/Penny/run03 "plot_* database" run03&
srun -n1 -p pdebug python runQueries.py /p/lscratche/geosadmn/Penny/run04 "plot_* database" run04&
srun -n1 -p pdebug python runQueries.py /p/lscratche/geosadmn/Penny/run05 "plot_* database" run05
srun -n1 -p pdebug python runQueries.py /p/lscratche/geosadmn/Penny/run06 "plot_* database" run06&
srun -n1 -p pdebug python runQueries.py /p/lscratche/geosadmn/Penny/run07 "plot_* database" run07&
srun -n1 -p pdebug python runQueries.py /p/lscratche/geosadmn/Penny/run08 "plot_* database" run08&
srun -n1 -p pdebug python runQueries.py /p/lscratche/geosadmn/Penny/run09 "plot_* database" run09&
