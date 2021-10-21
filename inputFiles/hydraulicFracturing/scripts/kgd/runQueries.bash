#!/bin/bash

#visit -cli -nowin -l msub/srun -nn 1 -np 8 -p pbatch -b cbronze -t 10:00    -s kgdQueries.py $PWD/siloFiles "zeroToughness_* database" zeroToughness
visit -cli -nowin -l msub/srun -nn 1 -np 16 -p pdebug                       -s kgdQueries.py $PWD/siloFiles "zeroToughness_* database" zeroToughness
visit -cli -nowin -l msub/srun -nn 1 -np 16 -p pdebug                       -s kgdQueries.py $PWD/siloFiles "zeroViscosity_* database" zeroViscosity
