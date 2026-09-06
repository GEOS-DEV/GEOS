#!/bin/bash 
# Michael E. Thomadakis Innovation and HPC R&D
# Prepare GPU environment on local node for MPI codes 
# Maps GPU numbers to ranks

# Assumes MPI has already been initialized!

: ${V:="0"};
: ${N_gpus:="0"};

H=$(hostname) 
HID="${H}:${$}"

if [ -z $N_gpus ]; then
    echo "## ${HID} $0 Error : N_gpus is undefined ! Exiting "
    exit 100;
elif [ $N_gpus -lt 1 ]; then
    echo "## ${HID} $0 Error : N_gpus is undefined ! Exiting "
    exit 100;
fi

if [ -z $MPI ]; then 
    export MPIRUN=$(which mpirun) 
    export MPI=$($MPIRUN -version | gawk '
      { if ($1=="Intel(R)" && $2=="MPI") {printf "impi"}
        else if ($2=="(Open" && $3=="MPI)") {printf "ompi"} }' )
fi

if [ "${MPI}" == "impi" ]; then
    # export MPIVER=$($MPIRUN -version | gawk '$1=="Intel(R)" && $2=="MPI" {printf "impi-%d.%d",$8, $10}')
    export GEOS_GLOBAL_RANK=$MPI_LOCALRANKID
    export GEOS_LOCAL_RANK=$PMI_RANK
else
    # export MPIVER=$($MPIRUN -version | gawk '$2=="(Open" && $3=="MPI)" {printf "ompi-%s",$4}')
    export GEOS_GLOBAL_RANK=$OMPI_COMM_WORLD_RANK
    export GEOS_LOCAL_RANK=$OMPI_COMM_WORLD_LOCAL_RANK
fi

export CUDA_VISIBLE_DEVICES=$(echo "${GEOS_LOCAL_RANK} % ${N_gpus}" | bc)
export GEOS_RANK_AFFINITY=$(taskset -c -p $$ | gawk '{print $6}')

export FHID="${HID}:${GEOS_GLOBAL_RANK}-${GEOS_LOCAL_RANK}"


if [ $V -ge 1 ]; then
    echo "## ${FHID} ................................................. "
    echo "#  ${FHID}  N_gpus= $N_gpus ; MPI= $MPI, MPIVER= $MPIVER, CUDA_VISIBLE_DEVICES= { $CUDA_VISIBLE_DEVICES }"

    echo "#  ${FHID}  GPUMPICLI= $GPUMPICLI"
    echo "#  ${FHID}  GEOS_RANK_AFFINITY=$GEOS_RANK_AFFINITY"
    if [ $V -ge 2 ]; then
	echo "## ${FHID} ................................................. "
	print "## ${FHID} {  "
	env | sort | gawk '{printf " %s,  ", $0; }'
	print "} ${FHID} \n"
	echo "## ${FHID} ................................................. "
    fi

    echo "## ${FHID} .................................................................................................... "
    echo "## ${FHID}  $GPUMPICLI $@"
    echo "## ${FHID} .................................................................................................... "

fi 

if [ $# -eq 0 ]; then
    $GPUMPICLI
else
   "$@"
fi 

## 
# CUDA_VISIBLE_DEVICES

# A comma-separated sequence of GPU identifiers MIG support: MIG-<GPU-UUID>/<GPU instance
# ID>/<compute instance ID>

# GPU identifiers are given as integer indices or as UUID strings. GPU UUID strings should
# follow the same format as given by nvidia-smi, such as
# GPU-8932f937-d72c-4106-c12f-20bd9faed9f6. However, for convenience, abbreviated forms are
# allowed; simply specify enough digits from the beginning of the GPU UUID to uniquely identify
# that GPU in the target system. For example, CUDA_VISIBLE_DEVICES=GPU-8932f937 may be a valid
# way to refer to the above GPU UUID, assuming no other GPU in the system shares this
# prefix. Only the devices whose index is present in the sequence are visible to CUDA
# applications and they are enumerated in the order of the sequence. If one of the indices is
# invalid, only the devices whose index precedes the invalid index are visible to CUDA
# applications. For example, setting CUDA_VISIBLE_DEVICES to 2,1 causes device 0 to be invisible
# and device 2 to be enumerated before device 1. Setting CUDA_VISIBLE_DEVICES to 0,2,-1,1 causes
# devices 0 and 2 to be visible and device 1 to be invisible. MIG format starts with MIG keyword
# and GPU UUID should follow the same format as given by nvidia-smi. For example,
# MIG-GPU-8932f937-d72c-4106-c12f-20bd9faed9f6/1/2. Only single MIG instance enumeration is
# supported.

## RANKS
## nvidia-smi top -m : GPU-2-GPU matrix 
#                 -c : CPU_id GPU list with affinity to core CPU_id

# IntelMPI ranks
# mpirun -prepend-rank -hostfile ~/hosts/A100_2.hosts -np 5 -ppn 3 bash -c "env | egrep -i 'rank' | sort "
# [4] MPI_LOCALNRANKS=2
# [4] MPI_LOCALRANKID=1
# [4] PMI_RANK=4


# OpenMPI ranks
# mpirun --tag-output --hostfile ~/hosts/A100_2.hosts --oversubscribe -np 5 --map-by ppr:3:node bash -c "env | egrep -i rank | sort "
# [1,4]<stdout>:OMPI_ARGV=-c env | egrep -i rank | sort 
# [1,4]<stdout>:OMPI_COMM_WORLD_LOCAL_RANK=1
# [1,4]<stdout>:OMPI_COMM_WORLD_NODE_RANK=1
# [1,4]<stdout>:OMPI_COMM_WORLD_RANK=4
# [1,4]<stdout>:OMPI_FIRST_RANKS=0
# [1,4]<stdout>:OMPI_MCA_orte_ess_node_rank=1
# [1,4]<stdout>:PMIX_RANK=4


