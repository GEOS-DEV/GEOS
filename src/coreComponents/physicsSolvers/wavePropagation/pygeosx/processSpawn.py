from fileManager import *
import multiprocessing as mp
from dask.distributed import Client, as_completed
from dask_jobqueue import SLURMCluster


def multiProcessing(shot_list, GEOSX, xml,
                    processes = 1,
                    mpiranks = None,
                    x = 1, y = 1, z = 1):
    p = []
    nb_shot_m1 = len(shot_list)
    ind = 0

    #Loop over the process launch
    for i in range(processes):
        nb_shot = int(nb_shot_m1/(processes - i))

        shot_file = exportShotList(i, shot_list[ind:ind + nb_shot])

        p.append( mp.Process(target = callProcess,
                             args   = (shot_file, GEOSX, xml, mpiranks, x, y, z)) )

        ind = ind + nb_shot
        nb_shot_m1 = nb_shot_m1 - nb_shot

        #Start process
        p[i].start()


    for i in range(number_of_processes):
        p[i].join()

    remove_shots_files()
def callProcess(shot_file, GEOSX, xml, mpiranks, x, y, z):

    wavePropagationPath = os.path.join(GEOSX, "GEOSX/src/coreComponents/physicsSolvers/wavePropagation/")
    pygeosxPath = os.path.join(wavePropagationPath, "pygeosx/")

    tracePath = os.path.join(rootPath, "outputSismoTrace/")
    if os.path.exists(tracePath):
        pass
    else:
        os.mkdir(tracePath)

    traceProcPath = os.path.join(tracePath, "traceProc" + str(mp.current_process()._identity)[1] + "/")
    if os.path.exists(traceProcPath):
        pass
    else:
        os.mkdir(traceProcPath)

    if mpiranks == None:
        os.system("python " + pygeosxPath + "geosxProcess.py -i " + xml + " " + shot_file + " " + traceProcPath)
    else:
        os.system("mpirun -np " + str(mpiranks) + " python " + pygeosxPath + "/geosxProcess.py -i " + xml + " -x " + str(x) + " -y " + str(y) + " -z " + str(z) + shot_file + " " + traceProcPath)




def daskProcessing(shot_list, GEOSX, xml,
                   memory = None,
                   cores = 1,
                   processes = 1,
                   nodes = 1,
                   x = 1, y = 1, z = 1,
                   cluster_type = "local"):
    futures = []
    nb_shot_m1 = len(shot_list)
    ind = 0

    wavePropagationPath = os.path.join(GEOSX, "GEOSX/src/coreComponents/physicsSolvers/wavePropagation/")
    pygeosxPath = os.path.join(wavePropagationPath, "pygeosx/")

    if cluster_type == "SLURM":
        cluster = SLURMcluster(cores = cores, memory = memory, processes = processes, job_cpu = cores)
        cluster.scale(nodes)
    elif cluster_type == "local":
        cluster = None

    client = Client(cluster)

    for i in range(nodes * processes):
        nb_shot = int(nb_shot_m1/(nodes*processes - i))

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
        fsh.write("cd " + GEOSX + "/GEOSX/build-test_module-release/lib/PYGEOSX/bin/ \n")
        if cluster_type == "SLURM":
            fsh.write("module load physics/pygeosx \n")
            fsh.write("module load physics/geosx_deps \n")
        fsh.write("source activate \n")
        fsh.write("cd " + pygeosxPath + " \n")
        if cores == 1:
            fsh.write("python geosxProcess.py -i " + xml + " " + shot_file + " " + traceProcPath)
        else:
            fsh.write("mpirun -np cores --oversubscribe python geosxProcess.py -i " + xml + " -x " + str(x) + " -y " + str(y) + " -z " + str(z) + " " + shot_file + " " + traceProcPath)
        fsh.close()

        futures.append(client.submit(os.system, "bash bash"+str(i)+".sh", key = "proc"+str(i)))
        ind = ind + nb_shot
        nb_shot_m1 = nb_shot_m1 - nb_shot

    for future, result in as_completed(futures, with_results=True):
        print(future, result)

    remove_shots_files()
    remove_bash_files()
    remove_dask_files()
    client.close()



def firstInit(GEOSX, xml,
              mpiranks = None,
              x = 1, y = 1, z = 1):
    wavePropagationPath = os.path.join(GEOSX, "GEOSX/src/coreComponents/physicsSolvers/wavePropagation/")
    pygeosxPath = os.path.join(wavePropagationPath, "pygeosx/")

    if mpiranks == None:
        os.system("python " + pygeosxPath + "firstInit.py -i " + xml)
    else:
        os.system("mpirun -np " + str(mpiranks) + " python firstInit.py -i " + xml + " -x " + str(x) + " -y " + str(y) + " -z " + str(z))

    maxT, dt = readInitVariable()

    return maxT, dt
