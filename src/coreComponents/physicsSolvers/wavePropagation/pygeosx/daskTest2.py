'''
Created on 9/02/2021

@author: macpro
'''

import pygeosx
import sys
from dask import delayed
from dask.distributed import Client, as_completed, fire_and_forget
from dask_jobqueue import SLURMCluster

from mesh import *
from acquisition import *
from shotFileManager import *
import glob

basePath = sys.argv[3] #Absolute path where the GEOSX folder is located

wavePropagationPath = os.path.join(basePath, "GEOSX/src/coreComponents/physicsSolvers/wavePropagation/")
pygeosxPath = os.path.join(wavePropagationPath, "pygeosx/")
segyPath = os.path.join(pygeosxPath, "segyAcquisition/")

xmlPath = sys.argv[2]


def multiProcessing(shot_list, nb_job = 1):
    futures = []
    nb_shot_m1 = len(shot_list)
    ind = 0

    #Loop over the process launch
    for i in range(nb_job):
        nb_shot = int(nb_shot_m1/(nb_job-i))

        shot_file = exportShotList(i, shot_list[ind:ind + nb_shot])
        
        tracePath = os.path.join(rootPath, "outputSismoTrace/")
        if os.path.exists(tracePath):
            pass
        else:
            os.mkdir(tracePath)

        traceProcPath = os.path.join(tracePath, "traceProc" + str(i) + "/")
        if os.path.exists(traceProcPath):
            pass
        else:
            os.mkdir(traceProcPath)
        
        fsh = open("bash"+str(i)+".sh", 'w+')
        fsh.write("#!/bin/bash \n")
        fsh.write("cd /beegfs/jbesset/codes/GEOSX/build-test_module-release/lib/PYGEOSX/bin/ \n")
        fsh.write("module load physics/pygeosx \n")                                                                                                                    
        fsh.write("module load physics/geosx_deps \n")
        fsh.write("source activate \n")
        fsh.write("cd /beegfs/jbesset/codes/GEOSX/src/coreComponents/physicsSolvers/wavePropagation/pygeosx/ \n")                                                      
        fsh.write("time mpirun -np $(($SLURM_CPUS_PER_TASK -2)) --oversubscribe python main.py -i " + str(sys.argv[2]) + " -x 6 -y 3 " + shot_file + " " + traceProcPath)
        fsh.close()            
        
        futures.append(client.submit(os.system, "bash bash"+str(i)+".sh", key = "proc"+str(i)))
        ind = ind + nb_shot
        nb_shot_m1 = nb_shot_m1 - nb_shot
        
    results = client.gather(futures)
    print(results)

def remove_tmp_files():
    for root, dir, files in os.walk(os.path.join(rootPath, "shots_lists")):
        for file in files:
            os.remove(os.path.join(root, file))

    os.rmdir(os.path.join(rootPath, "shots_lists"))

    bashFiles = glob.glob(os.path.join(rootPath, "bash*"))
    for file in bashFiles:
        os.remove(file)


def main():
    
    fsh = open("bash.sh", 'w+')
    fsh.write("#!/bin/bash \n")
    fsh.write("cd /beegfs/jbesset/codes/GEOSX/build-test_module-release/lib/PYGEOSX/bin/ \n")
    fsh.write("module load physics/pygeosx \n")
    fsh.write("module load physics/geosx_deps \n")
    fsh.write("source activate \n")
    fsh.write("cd /beegfs/jbesset/codes/GEOSX/src/coreComponents/physicsSolvers/wavePropagation/pygeosx/ \n")
    fsh.write("time mpirun -np 4 --oversubscribe python firstInit.py -i " + str(sys.argv[2]) + " -x 4")
    fsh.close()

    future= client.submit(os.system, "bash bash.sh", key = "firstInit")
    print(future.result())
    maxT, dt, boundary_box = readInitVariable()
    
    frequency = 5.0
    wavelet   = ricker(maxT, dt, frequency)


    shot_list = equispaced_acquisition(boundary_box,
                                       wavelet,
                                       dt,
                                       [1001, 1901, 4],
                                       [1751, 1751, 1],
                                       151,
                                       [21, 1981, 10],
                                       [1751, 1751, 1],
                                       101,
                                       export = 0
                                       )

    multiProcessing(shot_list, nb_job = 4)
    
    remove_tmp_files()

if __name__ == "__main__":

    cluster = SLURMCluster(cores = 20, 
                           memory = "200GB", 
                           processes = 1, 
                           job_cpu = 20,
                           job_extra=['-o jobDask.%j.%N.out', '-e jobDask.%j.%N.err'])
    cluster.scale(4)

    client = Client(cluster)

    main()
    
