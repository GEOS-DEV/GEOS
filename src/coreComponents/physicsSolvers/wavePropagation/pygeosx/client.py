import asyncssh
import asyncio
from tempfile import mkstemp
from munch import munchify
import uuid
import copy
import os
import subprocess
from time import sleep


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
        self.state="Not used"

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

        if extra is not None:
            header_lines.extend(["SBATCH %s" %arg for arg in extra])
        if env_extra is not None:
            header_lines.extend(["module load %s" %arg for arg in env_extra])

        self.header=header_lines

        if python is not None:
            self.python = python
        else:
            self.python = "python"


    def job_file(self):
        handle, filename = mkstemp(suffix=".sh")
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
            self.run += args[i] + " "


    def finalize_script(self):
        self.header.append(self.run)



class Client:
    def __init__(self,
                 cluster=None):

        self.cluster=cluster


    def submit(self,
               func,
               parameters,
               xml=None,
               cores=None,
               x_partition=1,
               y_partition=1,
               z_partition=1):

        module   = inspect.getmodule(func)
        filename = obj_to_json(parameters, filename = uuid.randomuuid().toString())
        args     = [module, func.__name__, filename]

        if cores is not None:
            if cores > self.cluster.cores:
                raise ValueError("Number of cores requested exceeds the number of cores available on cluster")
            else:
                self.cluster._add_job_to_cmd(cores)
                self.cluster._add_xml_to_cmd(xml)
                if x_partition*y_partition*z_partition != cores:
                    raise ValueError("Partition x y z must match x*y*z=cores")
                else:
                    self.cluster._add_partition_to_cmd(x_partition, y_partition, z_partition)
        else:
            raise ValueError("You must specify the number of cores you want to use")

        self.cluster._add_args_to_cmd(args)
        self.cluster.finalize_script()

        output = subprocess.check_output("usr/bin/sbatch " + self.cluster.job_file() + "sh")
        job_id = output.split()[-1]

        key=func.__name__ + "-" + str(uuid.uuid4())
        return Future(key, cluster=self.cluster, jobid=job_id)


    async def submit(self,
                     ind,
                     func,
                     parameters,
                     cores=None,
                     x_partition=1,
                     y_partition=1,
                     z_partition=1):
        while True:
            if self.cluster[ind].state == "Not used":
                module   = inspect.getmodule(func)
                filename = obj_to_json(parameters, filename = uuid.randomuuid().toString())
                args     = [module, func.__name__, filename]

                if cores is not None:
                    if cores > self.cluster[ind].cores:
                        raise ValueError("Number of cores requested exceeds the number of cores available on cluster")
                    else:
                        self.cluster[ind]._add_job_to_cmd(cores)
                        self.cluster[ind]._add_xml_to_cmd(xml)
                        if x_partition*y_partition*z_partition != cores:
                            raise ValueError("Partition x y z must match x*y*z=cores")
                        else:
                            self.cluster[ind]._add_partition_to_cmd(x_partition, y_partition, z_partition)
                else:
                    raise ValueError("You must specify the number of cores you want to use")

                self.cluster[ind]._add_args_to_cmd(args)
                self.cluster[ind].finalize_script()

                output = subprocess.check_output("usr/bin/sbatch " + self.cluster[ind].job_file() + "sh")
                job_id = output.split()[-1]

                key=func.__name__ + "-" + str(uuid.uuid4())
                return Future(key, cluster=self.cluster[ind], jobid=job_id)
            else:
                sleep(1)


    def scale(self, n=None):
        cluster_list=[self.cluster]
        conn_list=[self.conn]

        if n is not None:
            for i in range(n-1):
                cluster_list.append(copy.deepcopy(self.cluster))
                conn_list.append(self._connection())
        if len(cluster_list) > 1:
            self.cluster=cluster_list
            self.conn=conn_list


    def map(self,
            func,
            parameters,
            cores=None,
            x_partition=1,
            y_partition=1,
            z_partition=1):
        futures=[]
        task=[]
        nb_param=len(parameters)
        nb_cluster=len(self.cluster)
        for i in range(nb_param):
            ind = i%nb_cluster
            task.append(asyncio.create(self.submit(ind, func, parameters[i])))
        for i in range(nb_param):
            futures.append(asyncio.run(task[i]))
        return futures




class Future:
    def __init__(self,
                 key,
                 client=None):

        self.key=key
        self.cluster=cluster
        self.state="PENDING"
        self.job_id=None


    def result():
        while True:
            job_list = subprocess.check_output("squeue -o '%'A -h", shell=True).split()
            if self.job_id in job_list:
                status = subprocess.check_output("squeue -j %s -h -o %%t" %self.job_id, shell=True).strip()
                if status == "PD":
                    sleep(1)
                elif status == "R":
                    self.state = "RUNNING"
                    sleep(1)
                elif status in ["CG","CD"]:
                    self.state = "COMPLETED"
                    return 0
                else:
                    self.state = "ERROR"
                    return 1
            else:
                if os.stat(self.client.err).st_size == 0:
                    self.state="COMPLETED"
                    return 0
                else:
                    self.state="ERROR"
                    return 1





#use isinstance(obj, list) to check whether or not an obj is a list or not
