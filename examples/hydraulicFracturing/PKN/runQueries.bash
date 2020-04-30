#!/bin/bash

#visit -cli -nowin -l msub/srun -nn 1 -np 32 -p pbatch -b cbronze -t 10:00    -s pennyShapedQueries.py $PWD/siloFiles "pennyShaped_* database" pennyShaped
visit -cli -nowin -l msub/srun -nn 1 -np 32 -p pdebug                       -s PKNQueries.py $PWD/siloFiles "zeroToughness_* database" zeroToughness
