import sys
import os
sys.path.append(os.getcwd() + '/..')
from threading import Thread
from queue import Queue
import uuid
import copy
import subprocess
from time import sleep
import inspect
from utils import args_to_json
from tempfile import mkstemp

class Cluster():
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
                 python=None,
                 launch=None,
                 gpu=False):

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
        self.work=None
        self.launch=launch
        self.gpu=gpu

    def job_file(self):
        handle, filename = mkstemp(suffix=".sh", dir=os.getcwd())
        with open(filename, "w") as f:
            for line in self.header:
                f.write(line + "\n")

        os.close(handle)

        return filename


    def _add_job_to_script(self, cores):
        if self.launch=="mpirun":
            self.run = "%s -np %d " %(self.launch, cores)
        else:
            self.run = "%s -n %d " %(self.launch, cores)

        if self.gpu==True:
            self.run += "-g 1 "
        self.run += "%s seismicUtilities/wrapper.py " % self.python

    def _add_cmd_line(self, cmd_line):
        if isinstance(cmd_line, list):
            for cmd in cmd_line:
                self.run += "%s " %cmd
        else:
            self.run += "%s " %cmd_line


    def _add_partition_to_cmd(self, x, y, z):
        self.run += "-x %d -y %d -z %d" %(x, y, z) + " "


    def _add_args_to_cmd(self, args):
        for i in range(len(args)):
            self.run += str(args[i]) + " "


    def finalize_script(self):
        self.header.append(self.run)

    def setErrOut(self, job_id):
        self.err = "%s-%s-%s.err"%(self.type, self.job_name, job_id)
        self.out = "%s-%s-%s.out"%(self.type, self.job_name, job_id)



class SLURMCluster(Cluster):
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
                 python=None,
                 launch="srun",
                 gpu=False) :

        super().__init__(name,
                         job_name,
                         nodes,
                         cores,
                         memory,
                         walltime,
                         node_name,
                         job_cpu,
                         directory,
                         node_list,
                         extra,
                         env_extra,
                         python,
                         launch,
                         gpu)

        self.type="slurm"

        header_lines=["#!/usr/bin/bash"]
        if job_name is not None:
            header_lines.append("#SBATCH -J %s" % self.job_name)

        if nodes is not None:
            header_lines.append("#SBATCH -N %d" % self.nodes)

        if launch == "mpirun":
            header_lines.append("#SBATCH -n %d" % self.cores)
       # if cores is not None:
       #     header_lines.append("#SBATCH -n %d" % self.cores)
       # else:
       #raise ValueError("You must specify how much cores you want to use")

        if walltime is not None:
            header_lines.append("#SBATCH -t %s" % self.walltime)

        if directory is not None:
            header_lines.append("#SBATCH -D  %s" % self.directory)

        if node_name is not None:
            header_lines.append("#SBATCH -C %s" %self.node_name)

        header_lines.append("#SBATCH -o slurm-%s-%%A.out" %self.job_name)
        header_lines.append("#SBATCH -e slurm-%s-%%A.err" %self.job_name)

        if extra is not None:
            header_lines.extend(["#SBATCH %s" %arg for arg in extra])
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


    def free(self):
        if self.work is not None:
            job_list = subprocess.check_output("squeue -o %%A -h", shell=True).decode().split()
            if self.work.job_id not in job_list:
                self.state = "Free"
                self.header = self.header[:-1]
                self.run =  ""
                return True
            else:
                return False
        else:
            return True


    def getJobsList(self):
        job_list = subprocess.check_output("squeue -o %%A -h", shell=True).decode().split()
        return job_list


    def getJobState(self):
        status = subprocess.check_output("squeue -j %s -h -o %%t" %self.work.job_id, shell=True).decode().strip()
        return status


    def runJob(self, bashfile):
        output = subprocess.check_output("/usr/bin/sbatch "+ bashfile, shell=True)
        return output


class LSFCluster(Cluster):
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
                 python=None,
                 launch="lrun",
                 gpu=False):

        super().__init__(name,
                         job_name,
                         nodes,
                         cores,
                         memory,
                         walltime,
                         node_name,
                         job_cpu,
                         directory,
                         node_list,
                         extra,
                         env_extra,
                         python,
                         launch,
                         gpu)

        self.type="lsf"

        header_lines=["#!/usr/bin/tcsh"]
        if job_name is not None:
            header_lines.append("#BSUB -J %s" % self.job_name)

        if nodes is not None:
            header_lines.append("#BSUB -nnodes %d" % self.nodes)

        #if cores is not None:
        #    header_lines.append("#BSUB -n %d" % self.cores)
        #else:
        #    raise ValueError("You must specify how much cores you want to use")

        if walltime is not None:
            header_lines.append("#BSUB -W %s" % self.walltime)

        if directory is not None:
            header_lines.append("#BSUB -cwd  %s" % self.directory)

        if node_name is not None:
            header_lines.append("#BSUB -clusters %s" %self.node_name)

        header_lines.append("#BSUB -o lsf-%s-%%J.out" %self.job_name)
        header_lines.append("#BSUB -e lsf-%s-%%J.err" %self.job_name)

        if extra is not None:
            header_lines.extend(["#BSUB %s" %arg for arg in extra])
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


    def free(self):
        if self.work is not None:
            job_list = self.getJobsList()
            if self.work.job_id not in job_list:
                self.state = "Free"
                self.header = self.header[:-1]
                self.run =  ""
                return True
            else:
                return False
        else:
            return True


    def getJobsList(self):
        job_list = subprocess.check_output('bjobs -noheader -o JOBID', shell=True).decode().split()
        return job_list


    def getJobState(self):
        status = subprocess.check_output('bjobs -J %s -noheader -o "stat:"' %self.work.job_id, shell=True).decode().strip()
        return status


    def runJob(self, bashfile):
        output = subprocess.check_output("bsub " + bashfile, shell=True)
        return output


class Client:
    def __init__(self,
                 cluster=None):

        self.cluster = [cluster]
        self.output = self.create_output()



    def submit(self,
               func,
               args,
               cmd_line=None,
               cores=None,
               daemon=False,
               x_partition=1,
               y_partition=1,
               z_partition=1):

        future = ["In queue"]

        thread = Thread(target = self._submit,
                        args = (func,
                                args,
                                cmd_line,
                                cores,
                                x_partition,
                                y_partition,
                                z_partition,
                                future,),
                        daemon=daemon,
                        )
        thread.start()
        return future



    def map(self,
            func,
            args,
            cmd_line=None,
            parallel_tool="srun",
            cores=None,
            daemon=False,
            x_partition=1,
            y_partition=1,
            z_partition=1):

        futures = ["In queue"]*len(args)

        thread = Thread(target = self._map,
                        args = (func,
                                args,
                                cmd_line,
                                cores,
                                x_partition,
                                y_partition,
                                z_partition,
                                futures,),
                        daemon=daemon,
                    )
        thread.start()
        return futures



    def _submit(self,
                func,
                args,
                cmd_line,
                cores,
                x_partition,
                y_partition,
                z_partition,
                future):

        while True:
            for cluster in self.cluster:
                if cluster.free():
                    future[0] = self.run(func,
                                         args,
                                         cmd_line,
                                         cores,
                                         x_partition,
                                         y_partition,
                                         z_partition,
                                         cluster = cluster)
                    break
            else:
                sleep(1)
                continue
            break

        return future



    def _map(self,
             func,
             args,
             cmd_line,
             cores,
             x_partition,
             y_partition,
             z_partition,
             futures):

        work_queue = Queue()

        for arg in args:
            work_queue.put(arg)

            i=0
        while not work_queue.empty():
            work = work_queue.get()
            while True:
                for cluster in self.cluster:
                    if cluster.free():
                        futures[i] = self.run(func,
                                              work,
                                              cmd_line,
                                              cores,
                                              x_partition,
                                              y_partition,
                                              z_partition,
                                              cluster = cluster)
                        i+=1
                        break
                else:
                    sleep(1)
                    continue
                break



    def run(self,
            func,
            args,
            cmd_line,
            cores,
            x_partition,
            y_partition,
            z_partition,
            cluster):

        module   = inspect.getmodule(func).__name__
        key      = module + "-" + func.__name__ + "-" +str(x_partition) + "-" + str(y_partition) + "-" + str(z_partition) + "-" + str(uuid.uuid4())
        jsonfile = args_to_json(args)
        cmd_args = [module, func.__name__, jsonfile, self.output, key]

        if cores is not None:
            if cores > cluster.cores:
                raise ValueError("Number of cores requested exceeds the number of cores available on cluster")
            else:
                cluster._add_job_to_script(cores)
                if cmd_line is not None:
                    cluster._add_cmd_line(cmd_line)
                if x_partition*y_partition*z_partition != cores:
                    raise ValueError("Partition x y z must match x*y*z=cores")
                else:
                    cluster._add_partition_to_cmd(x_partition, y_partition, z_partition)
        else:
            raise ValueError("You must specify the number of cores you want to use")

        cluster._add_args_to_cmd(cmd_args)
        cluster.finalize_script()

        bashfile = cluster.job_file()
        output = cluster.runJob(bashfile)

        if cluster.type == "slurm":
            job_id = output.split()[-1].decode()
        else:
            job_id = output.split()[1][1:-1].decode()

        cluster.setErrOut(job_id)
        print("Job : " + job_id+ " has been submited")
        os.remove(bashfile)

        future = Future(key, output=self.output, cluster=cluster, args=args, job_id=job_id)

        cluster.work = future
        cluster.state = "Working"

        return future



    def scale(self, n=None):
        cluster = self.cluster[0]
        if n is not None:
            for i in range(n-1):
                self.cluster.append(copy.deepcopy(cluster))




    def create_output(self):
        output = os.path.join(os.getcwd(),"output")
        if os.path.exists(output):
            pass
        else:
            os.mkdir(output)

        return output



    def gather(self, futures, retry=False):
        n = len(futures)
        results = [0]*n
        while any([result == 0] for result in results):
            for i in range(n):
                if future[i].state in ["RUNNING", "PENDING"]:
                    results[i] = custom_future._result()
                elif future[i].state == "ERROR" and retry == True:
                    futures[i] = self.retry(custom_future)
                    print("Job", custom_future.job_id, "has failed, relaunch...")
                elif future[i].state == "ERROR" and retry == False:
                    results[i] = custom_future._result()

        return results




class Future:
    def __init__(self,
                 key,
                 output=None,
                 cluster=None,
                 args=None,
                 job_id=None):

        self.key=key
        self.args=args
        self.output=os.path.join(output, key+".txt")
        self.cluster=cluster
        self.state="PENDING"
        self.job_id = job_id
        self.err="%s-%s-%s.err"%(self.cluster.type, self.cluster.job_name, self.job_id)
        self.out="%s-%s-%s.out"%(self.cluster.type, self.cluster.job_name, self.job_id)


    def result(self):
        while self.state in ["RUNNING", "PENDING"]:
            result = self._result()
        return result


    def _result(self):
        job_list = self.cluster.getJobsList()
        if self.job_id in job_list:
            status = self.cluster.getJobState()
            if status in ["PD","PEND"]:
                sleep(1)
                return 0
            elif status in ["R","RUN"]:
                self.state = "RUNNING"
                sleep(1)
                return 0
            elif status in ["CD","DONE"]:
                self.state = "COMPLETED"
                result = self.getResult()
                print("Job " + self.job_id + " has been completed\n")
                return result
            else:
                self.state = "ERROR"
                print("Error in job " + self.job_id +" :\n")
                with open(self.err) as err:
                    for line in err:
                        if line[0] != "[":
                            print(line)
                return "err"
        else:
            with open(self.err) as err:
                if os.stat(self.err).st_size == 0:
                    self.state="COMPLETED"
                    result = self.getResult()
                elif (all([firstChar[0] for firstChar in err.readlines()]) == "[") == True:
                    self.state="COMPLETED"
                    result = self.getResult()
                    print("Job " + self.job_id + " has been completed\n")
                else:
                    self.state="ERROR"
                    err.seek(0)
                    print("Error in job " + self.job_id +" :\n")
                    for line in err.readlines():
                        if line[0] != "[":
                            print(line)

                    result = "err"
            err.close()
            return result


    def getResult(self):
        with open(self.output, 'r') as f:
            result = f.readlines()

        return result
