#!/bin/bash

visit -cli -nowin -l msub/srun -nn 1 -np 8 -p pbatch -b cbronze -t 10:00    -s pennyShapedQueries.py $PWD/siloFiles "pennyShaped_* database" pennyShaped
visit -cli -nowin -l msub/srun -nn 1 -np 16 -p pdebug                       -s pennyShapedQueries.py $PWD/siloFiles "zeroViscosity_* database" zeroViscosity
