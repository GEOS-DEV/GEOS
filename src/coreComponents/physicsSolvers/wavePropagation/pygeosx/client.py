import asyncio
import uuid
import copy
import subprocess
import numpy as np
from time import sleep
import inspect
from utils import *

class SLURMCluster:
    def __init__(self,
                 name=None,
                 job_name=None,
                 nodes=None,
                 cores=None,
                 memory=None,
                 walltime=None,
                 node_name=None,
                 job_cpu=None,
                 directory=None,
                 node_list=None,
                 extra=None,
                 env_extra=None,
                 python=None):

        self.name=name
        self.nodes=nodes
        self.cores=cores
        self.job_cpu=job_cpu
        self.job_name=job_name
        self.walltime=walltime
        self.directory=directory
        self.node_name=node_name
        self.python=python
        self.err=None
        self.out=None
        self.state="Free"

        header_lines=["#!/usr/bin/bash"]
        if job_name is not None:
            header_lines.append("#SBATCH -J %s" % self.job_name)

        if nodes is not None:
            header_lines.append("#SBATCH -N %d" % self.nodes)

        if cores is not None:
            header_lines.append("#SBATCH -n %d" % self.cores)
        else:
            raise ValueError("You must specify how much cores you want to use")

        if walltime is not None:
            header_lines.append("#SBATCH -t %s" % self.walltime)

        if directory is not None:
            header_lines.append("#SBATCH -D  %s" % self.directory)

        if node_name is not None:
            header_lines.append("#SBATCH -C %s" %self.node_name)

        header_lines.append("#SBATCH -o slurm-%s-%%A.out" %self.job_name)
        header_lines.append("#SBATCH -e slurm-%s-%%A.err" %self.job_name)

        if extra is not None:
            header_lines.extend(["SBATCH %s" %arg for arg in extra])
        if env_extra is not None:
            header_lines.extend(["module load %s" %arg for arg in env_extra])

        self.header=header_lines

        if python is not None:
            if python=="$Python3_EXECUTABLE":
                self.header.append("source /beegfs/jbesset/codes/GEOSX/build-test_module-release/lib/PYGEOSX/bin/activate")
                self.python="python"
            else:
                self.python = python
        else:
            self.python = "python"


    def job_file(self):
        handle, filename = mkstemp(suffix=".sh", dir=os.getcwd())
        with open(filename, "w") as f:
            for line in self.header:
                f.write(line + "\n")
        
        os.close(handle)
        return filename


    def _add_job_to_script(self, cores):
        if cores > 1:
            self.run = "mpirun -np %d --map-by node %s startJob.py " %(cores, self.python)
        else:
            self.run = "%s startJob.py " %self.python


    def _add_xml_to_cmd(self, xml):
        self.run += "-i %s" %xml + " "


    def _add_partition_to_cmd(self, x, y, z):
        self.run += "-x %d -y %d -z %d" %(x, y, z) + " "


    def _add_args_to_cmd(self, args):
        for i in range(len(args)):
            self.run += str(args[i]) + " "


    def finalize_script(self):
        self.header.append(self.run)


class Client:
    def __init__(self,
                 cluster=None):

        self.cluster = cluster
        self.output = self.create_output()


    def submit(self,
               func,
               parameters,
              # xml=None,
               cores=None,
               x_partition=1,
               y_partition=1,
               z_partition=1):

        key = func.__name__ + "-" + str(uuid.uuid4())
        module   = inspect.getmodule(func).__name__
        jsonfile = obj_to_json(parameters)
        args     = [module, func.__name__, jsonfile, self.output, key]

        if cores is not None:
            if cores > self.cluster.cores:
                raise ValueError("Number of cores requested exceeds the number of cores available on cluster")
            else:
                self.cluster._add_job_to_script(cores)
               # self.cluster._add_xml_to_cmd(xml)
                if x_partition*y_partition*z_partition != cores:
                    raise ValueError("Partition x y z must match x*y*z=cores")
                else:
                    self.cluster._add_partition_to_cmd(x_partition, y_partition, z_partition)
        else:
            raise ValueError("You must specify the number of cores you want to use")

        self.cluster._add_args_to_cmd(args)
        self.cluster.finalize_script()
        
        bashfile = self.cluster.job_file()

        if isinstance(self.cluster, SLURMCluster):
            output = subprocess.check_output("/usr/bin/sbatch " + bashfile, shell=True)

        job_id = output.split()[-1].decode()

        print("Job : " + job_id+ " has been submited")
        os.remove(bashfile)

        return Future(key, output=self.output, cluster=self.cluster, job_id=job_id)


    def _submit(self,
                ind,
                func,
                parameters,
                #xml=None,
                cores=None,
                x_partition=1,
                y_partition=1,
                z_partition=1):

        key = func.__name__ + "-" + str(uuid.uuid4())
        module   = inspect.getmodule(func).__name__
        jsonfile = obj_to_json(parameters)
        args     = [module, func.__name__, jsonfile, self.output, key]

        if cores is not None:
            if cores > self.cluster[ind].cores:
                raise ValueError("Number of cores requested exceeds the number of cores available on cluster")
            else:
                self.cluster[ind]._add_job_to_script(cores)
                #self.cluster[ind]._add_xml_to_cmd(xml)
                if x_partition*y_partition*z_partition != cores:
                    raise ValueError("Partition x y z must match x*y*z=cores")
                else:
                    self.cluster[ind]._add_partition_to_cmd(x_partition, y_partition, z_partition)
        else:
            raise ValueError("You must specify the number of cores you want to use")

        self.cluster[ind]._add_args_to_cmd(args)
        self.cluster[ind].finalize_script()
        
        bashfile = self.cluster[ind].job_file()
        output = subprocess.check_output("/usr/bin/sbatch " + bashfile, shell=True)
        job_id = output.split()[-1].decode()
        
        print("Job : " + job_id + " has been submited")
        os.remove(bashfile)
        
        return Future(key, output=self.output, cluster=self.cluster[ind], job_id=job_id)



    def scale(self, n=None):
        cluster_list = [self.cluster]
        
        if n is not None:
            for i in range(n-1):
                cluster_list.append(copy.deepcopy(self.cluster))
        if len(cluster_list) > 1:
            self.cluster=cluster_list



    def map(self,
            func,
            parameters,
            #xml=None,
            cores=None,
            x_partition=1,
            y_partition=1,
            z_partition=1):
        
        futures = []
        nb_param = len(parameters)
        nb_cluster = len(self.cluster)

        for i in range(nb_param):
            ind = i % nb_cluster
            while True:
                if self.cluster[ind].state == "Free":
                    futures.append(self._submit(ind, 
                                                func, 
                                                parameters[i], 
                                                #xml,
                                                cores,
                                                x_partition,
                                                y_partition,
                                                z_partition))
                    self.cluster[ind].state = "Working"
                    break
                else:
                    sleep(1)
                    futures[ind].check_state()
            
        return futures

        
    def create_output(self):
        output = os.path.join(os.getcwd(),"output")
        if os.path.exists(output):
            pass
        else:
            os.mkdir(output)
        
        return output

    
    def gather(self, futures):
        results = []
        for future in futures:
            results.append(future.result())

        return results



class Future:
    def __init__(self,
                 key,
                 output=None,
                 cluster=None,
                 job_id=None):

        self.key=key
        self.output=os.path.join(output, key+".txt")
        self.cluster=cluster
        self.state="PENDING"
        self.job_id = job_id
        self.err="slurm-%s-%s.err"%(self.cluster.job_name, self.job_id)
        self.out="slurm-%s-%s.out"%(self.cluster.job_name, self.job_id)


    def result(self):
        while True:
            job_list = subprocess.check_output("squeue -o '%'A -h", shell=True).decode().split()
            if self.job_id in job_list:
                status = subprocess.check_output("squeue -j %s -h -o %%t" %self.job_id, shell=True).decode().strip()
                if status == "PD":
                    sleep(1)
                elif status == "R":
                    self.state = "RUNNING"
                    sleep(1)
                elif status in ["CG","CD"]:
                    self.state = "COMPLETED"
                    result = self.getResult()
                    return result
                else:
                    self.state = "ERROR"
                    return 1
            else:
                with open(self.err) as err:
                    if os.stat(self.err).st_size == 0:
                        self.state="COMPLETED"
                        result = self.getResult()
                    elif (np.array([firstChar[0] for firstChar in err.readlines()]) == "[").all() != True:
                        self.state="COMPLETED"
                        result = self.getResult()
                    else:
                        self.state="ERROR"
                        result = 1
                err.close()
                return result


    
    def check_state(self):
        job_list = subprocess.check_output("squeue -o '%'A -h", shell=True).decode().split()
        if self.job_id not in job_list:
            self.cluster.state = "Free"
            self.cluster.header = self.cluster.header[:-1]
            self.run =  ""

    
    def getResult(self):
        with open(self.output, 'r') as f:
            result = f.readlines()
        
        return result
        
