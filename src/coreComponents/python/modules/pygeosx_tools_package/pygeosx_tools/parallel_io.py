
from mpi4py import MPI
import numpy as np


# Get the MPI rank
comm = MPI.COMM_WORLD
rank = comm.Get_rank()


def get_global_array_range(local_values):
    # 1D arrays will return a scalar, ND arrays an array
    N = np.shape(local_values)
    local_min = 1e100
    local_max = -1e100
    if (len(N) > 1):
        local_min = np.zeros(N[1]) + 1e100
        local_max = np.zeros(N[1]) - 1e100

    # For >1D arrays, keep the last dimension
    query_axis = 0
    if (len(N) > 2):
        query_axis = tuple([ii for ii in range(0, len(N)-1)])

    # Ignore zero-length results
    if len(local_values):
        local_min = np.amin(local_values, axis=query_axis)
        local_max = np.amax(local_values, axis=query_axis)

    # Gather the results onto rank 0
    all_min = comm.gather(local_min, root=0)
    all_max = comm.gather(local_max, root=0)
    global_min = 1e100
    global_max = -1e100
    if (rank == 0):
        global_min = np.amin(np.array(all_min), axis=0)
        global_max = np.amax(np.array(all_max), axis=0)
    return global_min, global_max


def gather_array(local_values, allgather=False, concatenate=True):
    # Find buffer size
    N = np.shape(local_values)
    M = np.prod(N)
    all_M = []
    max_M = 0
    if allgather:
        all_M = comm.allgather(M)
        max_M = np.amax(all_M)
    else:
        all_M = comm.gather(M, root=0)
        if (rank == 0):
            max_M = np.amax(all_M)
        max_M = comm.bcast(max_M, root=0)

    # Pack the array into a buffer
    send_buff = np.zeros(max_M)
    send_buff[:M] = np.reshape(local_values, (-1))
    receive_buff = np.zeros((comm.size, max_M))

    # Gather the buffers
    if allgather:
        comm.Allgather([send_buff, MPI.DOUBLE],
                       [receive_buff, MPI.DOUBLE])
    else:
        comm.Gather([send_buff, MPI.DOUBLE],
                    [receive_buff, MPI.DOUBLE],
                    root=0)

    # Unpack the buffers
    all_values = []
    R = list(N)
    R[0] = -1
    if ((rank == 0) | allgather):
        # Reshape each rank's contribution
        for ii in range(comm.size):
            if (all_M[ii] > 0):
                tmp = np.reshape(receive_buff[ii, :all_M[ii]], R)
                all_values.append(tmp)

        # Concatenate into a single array
        if concatenate:
            all_values = np.concatenate(all_values, axis=0)

    return all_values
