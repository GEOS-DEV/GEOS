#!/bin/bash

visit -cli -nowin -l msub/srun -nn 1 -np 8 -p pbatch -b cbronze -t 30:00    -s kgdQueries.py $PWD/siloFiles "zeroToughness_* database" zeroToughness
