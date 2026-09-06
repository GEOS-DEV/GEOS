#!/bin/bash
#
# Example SLURM batch script for GEOS using OpenMPI .
#
# Chevron Technology Center, Innovation and HPC R&D 
#
# ** You have to modify the relevant SBATCH options to match the resource
#     allocation requirements of your HPC job
#
# Michael E Thomadakis

# The general syntax of a Slurm job script is 
#SLURM_Keyword SLURM_Options ...

# Account
#SBATCH --account=rso

# SLURM Partition (queue) for GEOS: rso and geos
#SBATCH -p rso 

# Azure SKU types: hbv2, hbv3, hbv4 or hc44 (hbv3 are AMD Zen3, hbv2 are older AMD Zen2 and
#   hc44 are the ancient Intel Skylake
#SBATCH --constraint=hbv3 

# Number of nodes of type --constraint=
##SBATCH -N N_nodes 
#SBATCH -N 4 

# N_tasks : number of "tasks", i.e., MPI ranks; can be ommited 
##SBATCH -n N_tasks 

# Request exclusive use of the nodes; need to do for perf tests
#SBATCH --exclusive 


# number of cores per tasks if running OMP use -c >= OMP_NUM_THREADS 
#SBATCH -c 2 

# Time limits; if the job exheeds this amont. SLURM will kill it
#SBATCH -t 0:30:00 # time (D-HH:MM)

# Where to save STDOU and STDERR 
#SBATCH -o geos.%N.%j.out # STDOUT
#SBATCH -e geos.%N.%j.err # STDERR

# Miscaleneous 
# memory pool for all cores; not really needed
##SBATCH --mem 100 
 
# Example GEOS run for OPTO2/OpenMPI/OMP support, CPU-only binary
#
# Source proper initialization .rc file 
#
# 
# To see all Open MPI rc files in public GEOS installations use : 
#   ls -l /chap/geos/.init/GEOS/*CPU*ompi* 
# or 
#   ls -l /data/saet/software/.init/GEOS/*CPU*ompi* 
# 
# For GEOS built by you use 
#   ls -l ${HOME}/.init/GEOS/*CPU*ompi*

# Example
 
# Using INITRC to initialize the SLURM environment
: ${INITRC:="/data/saet/software/.init/GEOS/Init-x86_64-RHEL8-GEOS-0.2.0-CPU-OPTO3-Hypre-GCC_13.2.0-ompi_hpcx-OMP-relwithdebinfo.rc"}; 
export INITRC 
echo "Slurm $0 : Initializing GEOS using INITRC=$INITRC"
source $INITRC 


# Show loaded modules and location of geosx binary
module list 
which geosx mpirun

# cd to /where/the/model-xml/is
cd /where/the/model-xml/is

# Show current directory
pwd

# Example with 2 OMP thread 120 ranks total, 2 cores per rank and 30 ranks / node 
: ${OMP_NUM_THREADS:=2}
export OMP_NUM_THREADS

mpirun -np 120 --map-by ppr:30:node:PE=2  \
    $GEOSX_DIR/bin/geosx -i ./xxxx.xml -x X -y Y -z Z


## sbatch Documentation 
## Usage: sbatch [OPTIONS(0)...] [ : [OPTIONS(N)...]] script(0) [args(0)...]

# Parallel run options:
#   -a, --array=indexes         job array index values
#   -A, --account=name          charge job to specified account
#       --bb=<spec>             burst buffer specifications
#       --bbf=<file_name>       burst buffer specification file
#   -b, --begin=time            defer job until HH:MM MM/DD/YY
#       --comment=name          arbitrary comment
#       --cpu-freq=min[-max[:gov]] requested cpu frequency (and governor)
#   -c, --cpus-per-task=ncpus   number of cpus required per task
#   -d, --dependency=type:jobid[:time] defer job until condition on jobid is satisfied
#       --deadline=time         remove the job if no ending possible before
#                               this deadline (start > (deadline - time[-min]))
#       --delay-boot=mins       delay boot for desired node features
#   -D, --chdir=directory       set working directory for batch script
#   -e, --error=err             file for batch script's standard error
#       --export[=names]        specify environment variables to export
#       --export-file=file|fd   specify environment variables file or file
#                               descriptor to export
#       --get-user-env          load environment from local cluster
#       --gid=group_id          group ID to run job as (user root only)
#       --gres=list             required generic resources
#       --gres-flags=opts       flags related to GRES management
#   -H, --hold                  submit job in held state
#       --ignore-pbs            Ignore #PBS and #BSUB options in the batch script
#   -i, --input=in              file for batch script's standard input
#   -J, --job-name=jobname      name of job
#   -k, --no-kill               do not kill job on node failure
#   -L, --licenses=names        required license, comma separated
#   -M, --clusters=names        Comma separated list of clusters to issue
#                               commands to.  Default is current cluster.
#                               Name of 'all' will submit to run on all clusters.
#                               NOTE: SlurmDBD must up.
#       --container             Path to OCI container bundle
#   -m, --distribution=type     distribution method for processes to nodes
#                               (type = block|cyclic|arbitrary)
#       --mail-type=type        notify on state change: BEGIN, END, FAIL or ALL
#       --mail-user=user        who to send email notification for job state
#                               changes
#       --mcs-label=mcs         mcs label if mcs plugin mcs/group is used
#   -n, --ntasks=ntasks         number of tasks to run
#       --nice[=value]          decrease scheduling priority by value
#       --no-requeue            if set, do not permit the job to be requeued
#       --ntasks-per-node=n     number of tasks to invoke on each node
#   -N, --nodes=N               number of nodes on which to run (N = min[-max])
#   -o, --output=out            file for batch script's standard output
#   -O, --overcommit            overcommit resources
#   -p, --partition=partition   partition requested
#       --parsable              outputs only the jobid and cluster name (if present),
#                               separated by semicolon, only on successful submission.
#       --power=flags           power management options
#       --priority=value        set the priority of the job to value
#       --profile=value         enable acct_gather_profile for detailed data
#                               value is all or none or any combination of
#                               energy, lustre, network or task
#       --propagate[=rlimits]   propagate all [or specific list of] rlimits
#   -q, --qos=qos               quality of service
#   -Q, --quiet                 quiet mode (suppress informational messages)
#       --reboot                reboot compute nodes before starting job
#       --requeue               if set, permit the job to be requeued
#   -s, --oversubscribe         over subscribe resources with other jobs
#   -S, --core-spec=cores       count of reserved cores
#       --signal=[[R][B]:]num[@time] send signal when time limit within time seconds
#       --spread-job            spread job across as many nodes as possible
#       --switches=max-switches{@max-time-to-wait}
#                               Optimum switches and max time to wait for optimum
#       --thread-spec=threads   count of reserved threads
#   -t, --time=minutes          time limit
#       --time-min=minutes      minimum time limit (if distinct)
#       --uid=user_id           user ID to run job as (user root only)
#       --use-min-nodes         if a range of node counts is given, prefer the
#                               smaller count
#   -v, --verbose               verbose mode (multiple -v's increase verbosity)
#   -W, --wait                  wait for completion of submitted job
#       --wckey=wckey           wckey to run job under
#       --wrap[=command string] wrap command string in a sh script and submit

# Constraint options:
#       --cluster-constraint=[!]list specify a list of cluster constraints
#       --contiguous            demand a contiguous range of nodes
#   -C, --constraint=list       specify a list of constraints
#   -F, --nodefile=filename     request a specific list of hosts
#       --mem=MB                minimum amount of real memory
#       --mincpus=n             minimum number of logical processors (threads)
#                               per node
#       --reservation=name      allocate resources from named reservation
#       --tmp=MB                minimum amount of temporary disk
#   -w, --nodelist=hosts...     request a specific list of hosts
#   -x, --exclude=hosts...      exclude a specific list of hosts

# Consumable resources related options:
#       --exclusive[=user]      allocate nodes in exclusive mode when
#                               cpu consumable resource is enabled
#       --exclusive[=mcs]       allocate nodes in exclusive mode when
#                               cpu consumable resource is enabled
#                               and mcs plugin is enabled
#       --mem-per-cpu=MB        maximum amount of real memory per allocated
#                               cpu required by the job.
#                               --mem >= --mem-per-cpu if --mem is specified.

# Affinity/Multi-core options: (when the task/affinity plugin is enabled)
#                               For the following 4 options, you are
#                               specifying the minimum resources available for
#                               the node(s) allocated to the job.
#       --sockets-per-node=S    number of sockets per node to allocate
#       --cores-per-socket=C    number of cores per socket to allocate
#       --threads-per-core=T    number of threads per core to allocate
#   -B  --extra-node-info=S[:C[:T]]  combine request of sockets per node,
#                               cores per socket and threads per core.
#                               Specify an asterisk (*) as a placeholder,
#                               a minimum value, or a min-max range.

#       --ntasks-per-core=n     number of tasks to invoke on each core
#       --ntasks-per-socket=n   number of tasks to invoke on each socket
#       --hint=                 Bind tasks according to application hints
#                               (see "--hint=help" for options)
#       --mem-bind=             Bind memory to locality domains (ldom)
#                               (see "--mem-bind=help" for options)

# GPU scheduling options:
#       --cpus-per-gpu=n        number of CPUs required per allocated GPU
#   -G, --gpus=n                count of GPUs required for the job
#       --gpu-bind=...          task to gpu binding options
#       --gpu-freq=...          frequency and voltage of GPUs
#       --gpus-per-node=n       number of GPUs required per allocated node
#       --gpus-per-socket=n     number of GPUs required per allocated socket
#       --gpus-per-task=n       number of GPUs required per spawned task
#       --mem-per-gpu=n         real memory required per allocated GPU

# Help options:
#   -h, --help                  show this help message
#       --usage                 display brief usage message

# Other options:
#   -V, --version               output version information and exit
