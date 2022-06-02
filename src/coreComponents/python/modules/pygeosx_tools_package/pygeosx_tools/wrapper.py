import sys
import numpy as np
from mpi4py import MPI
import matplotlib.pyplot as plt
import pylvarray
import pygeosx


# Get the MPI rank
comm = MPI.COMM_WORLD
rank = comm.Get_rank()


def get_wrapper(problem, target_key, write_flag=False):
    """
    Get a local copy of a wrapper as a numpy ndarray

    Args:
        filename (str): Catalog file name
        problem (pygeosx.Group): GEOSX problem handle
        target_key (str): Key for the target wrapper
        write_flag (bool): Sets write mode (default=False)

    Returns:
        np.ndarray: The wrapper as a numpy ndarray
    """
    local_values = problem.get_wrapper(target_key).value()

    if hasattr(local_values, "set_access_level"):
        # Array types will have the set_access_level method
        # These require additional manipulation before use
        if write_flag:
            local_values.set_access_level(pylvarray.MODIFIABLE, pylvarray.CPU)
        else:
            local_values.set_access_level(pylvarray.CPU)

        if hasattr(local_values, "to_numpy"):
            local_values = local_values.to_numpy()
    return local_values


def get_wrapper_par(problem, target_key, allgather=False, ghost_key=""):
    """
    Get a global copy of a wrapper as a numpy ndarray.
    Note: if ghost_key is set, it will try to remove any ghost elements

    Args:
        problem (pygeosx.Group): GEOSX problem handle
        target_key (str): Key for the target wrapper
        allgather (bool): Flag to trigger allgather across ranks (False)
        ghost_key (str): Key for the corresponding ghost wrapper (default='')

    Returns:
        np.ndarray: The wrapper as a numpy ndarray
    """
    if comm.size == 1:
        # This is a serial problem
        return get_wrapper(problem, target_key)

    else:
        # This is a parallel problem
        # Get the local wrapper size, shape
        local_values = get_wrapper(problem, target_key)

        # Filter out ghost ranks if requested
        if ghost_key:
            ghost_values = get_wrapper(problem, ghost_key)
            local_values = local_values[ghost_values < -0.5]

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
            if rank == 0:
                max_M = np.amax(all_M)
            max_M = comm.bcast(max_M, root=0)

        # Pack the array into a buffer
        send_buff = np.zeros(max_M)
        send_buff[:M] = np.reshape(local_values, (-1))
        receive_buff = np.zeros((comm.size, max_M))

        # Gather the buffers
        if allgather:
            comm.Allgather([send_buff, MPI.DOUBLE], [receive_buff, MPI.DOUBLE])
        else:
            comm.Gather([send_buff, MPI.DOUBLE], [receive_buff, MPI.DOUBLE], root=0)

        # Unpack the buffers
        all_values = []
        R = list(N)
        R[0] = -1
        if (rank == 0) | allgather:
            # Reshape each rank's contribution
            for ii in range(comm.size):
                if all_M[ii] > 0:
                    tmp = np.reshape(receive_buff[ii, : all_M[ii]], R)
                    all_values.append(tmp)

            # Concatenate into a single array
            all_values = np.concatenate(all_values, axis=0)
        return all_values


def gather_wrapper(problem, key, ghost_key=""):
    """
    Get a global copy of a wrapper as a numpy ndarray on rank 0

    Args:
        problem (pygeosx.Group): GEOSX problem handle
        target_key (str): Key for the target wrapper

    Returns:
        np.ndarray: The wrapper as a numpy ndarray
    """
    return get_wrapper_par(problem, key, ghost_key=ghost_key)


def allgather_wrapper(problem, key, ghost_key=""):
    """
    Get a global copy of a wrapper as a numpy ndarray on all ranks

    Args:
        problem (pygeosx.Group): GEOSX problem handle
        target_key (str): Key for the target wrapper

    Returns:
        np.ndarray: The wrapper as a numpy ndarray
    """
    return get_wrapper_par(problem, key, allgather=True, ghost_key=ghost_key)


def get_global_value_range(problem, key):
    """
    Get the range of a target value across all processes

    Args:
        problem (pygeosx.Group): GEOSX problem handle
        target_key (str): Key for the target wrapper

    Returns:
        tuple: The global min/max of the target
    """
    local_values = get_wrapper(problem, key)

    # 1D arrays will return a scalar, ND arrays an array
    N = np.shape(local_values)
    local_min = 1e100
    local_max = -1e100
    if len(N) > 1:
        local_min = np.zeros(N[1]) + 1e100
        local_max = np.zeros(N[1]) - 1e100

    # For >1D arrays, keep the last dimension
    query_axis = 0
    if len(N) > 2:
        query_axis = tuple([ii for ii in range(0, len(N) - 1)])

    # Ignore zero-length results
    if len(local_values):
        local_min = np.amin(local_values, axis=query_axis)
        local_max = np.amax(local_values, axis=query_axis)

    # Gather the results onto rank 0
    all_min = comm.gather(local_min, root=0)
    all_max = comm.gather(local_max, root=0)
    global_min = 1e100
    global_max = -1e100
    if rank == 0:
        global_min = np.amin(np.array(all_min), axis=0)
        global_max = np.amax(np.array(all_max), axis=0)
    return global_min, global_max


def print_global_value_range(problem, key, header, scale=1.0, precision="%1.4f"):
    """
    Print the range of a target value across all processes

    Args:
        problem (pygeosx.Group): GEOSX problem handle
        target_key (str): Key for the target wrapper
        header (str): Header to print with the range
        scale (float): Multiply the range with this value before printing (default = 1.0)
        precision (str): Format for printing the range (default = %1.4f)

    Returns:
        tuple: The global min/max of the target
    """
    global_min, global_max = get_global_value_range(problem, key)
    global_min *= scale
    global_max *= scale

    if rank == 0:
        if isinstance(global_min, np.ndarray):
            min_str = ", ".join([precision % (x) for x in global_min])
            max_str = ", ".join([precision % (x) for x in global_max])
            print("%s: min=[%s], max=[%s]" % (header, min_str, max_str))
        else:
            min_str = precision % (global_min)
            max_str = precision % (global_max)
            print("%s: min=%s, max=%s" % (header, min_str, max_str))

    # Return a copy of the min/max in case we want to use them
    return global_min, global_max


def set_wrapper_to_value(problem, key, value):
    """
    Set the value of a wrapper

    Args:
        problem (pygeosx.Group): GEOSX problem handle
        target_key (str): Key for the target wrapper
        value (float): Value to set the wrapper
    """
    local_values = get_wrapper(problem, key, write_flag=True)
    local_values[...] = value


def set_wrapper_with_function(problem, target_key, input_keys, fn, target_index=-1):
    """
    Set the value of a wrapper using a function

    Args:
        problem (pygeosx.Group): GEOSX problem handle
        target_key (str): Key for the target wrapper
        input_keys (str, list): The input key(s)
        fn (function): Vectorized function used to calculate target values
        target_index (int): Target index to write the output (default = all)
    """
    if isinstance(input_keys, str):
        input_keys = [input_keys]
    local_target = get_wrapper(problem, target_key, write_flag=True)
    local_inputs = [get_wrapper(problem, k) for k in input_keys]

    # Run the function, check the shape of outputs/target
    fn_output = fn(*local_inputs)
    N = np.shape(local_target)
    M = np.shape(fn_output)

    if target_index < 0:
        if N == M:
            # Function output, target shapes are the same
            local_target[...] = fn_output

        elif len(M) == 1:
            # Apply the function output across each of the target dims
            local_target[...] = np.tile(np.expand_dims(fn_output, axis=1), (1, N[1]))

        else:
            raise Exception("Shape of function output %s is not compatible with target %s" % (str(M), str(N)))
    elif len(M) == 1:
        if len(N) == 2:
            # 2D target, with 1D output applied to a given index
            local_target[:, target_index] = fn_output

        else:
            # ND target, with 1D output tiled across intermediate indices
            expand_axes = tuple([ii for ii in range(1, len(N) - 1)])
            tile_axes = tuple([1] + [ii for ii in N[1:-1]])
            local_target[..., target_index] = np.tile(np.expand_dims(fn_output, axis=expand_axes), tile_axes)

    else:
        raise Exception(
            "Shape of function output %s is not compatible with target %s (target axis=%i)"
            % (str(M), str(N), target_index)
        )


def search_datastructure_wrappers_recursive(group, filters, matching_paths, level=0, group_path=[]):
    """
    Recursively search the group and its children for wrappers that match the filters

    Args:
        problem (pygeosx.Group): GEOSX problem handle
        filters (list): a list of strings
        matching_paths (list): a list of matching values
    """
    for wrapper in group.wrappers():
        wrapper_path = str(wrapper).split()[0]
        wrapper_test = group_path + [wrapper_path[wrapper_path.rfind("/") + 1 :]]
        if all(f in wrapper_test for f in filters):
            matching_paths.append("/".join(wrapper_test))

    for sub_group in group.groups():
        sub_group_name = str(sub_group).split()[0].split("/")[-1]
        search_datastructure_wrappers_recursive(
            sub_group, filters, matching_paths, level=level + 1, group_path=group_path + [sub_group_name]
        )


def get_matching_wrapper_path(problem, filters):
    """
    Recursively search the group and its children for wrappers that match the filters
    A successful match is identified if the wrapper path contains all of the
    strings in the filter.  Note: multiple matches are possible
        For example, if filters=['a', 'b', 'c'], the following would match
            - a/b/c
            - c/b/a
            - d/e/c/f/b/a/a

    Args:
        problem (pygeosx.Group): GEOSX problem handle
        filters (list): a list of strings

    Returns:
        str: Key of the matching wrapper

    """
    matching_paths = []
    search_datastructure_wrappers_recursive(problem, filters, matching_paths)

    if len(matching_paths) == 1:
        if rank == 0:
            print("Found matching wrapper: %s" % (matching_paths[0]))
        return matching_paths[0]

    else:
        if rank == 0:
            print("Error occured while looking for wrappers:")
            print("Filters: [%s]" % (", ".join(filters)))
            print("Matching wrappers: [%s]" % (", ".join(matching_paths)))
        raise Exception("Search resulted in 0 or >1 wrappers mathching filters")


def run_queries(problem, records):
    """
    Query the current GEOSX datastructure
    Note: The expected record request format is as follows.
    For now, the only supported query is to find the min/max values of the target
        record = {'path/of/wrapper': {'label': 'aperture (m)', # A label to include with plots
                  'scale': 1.0,  # Value to scale results by
                  'history: [],  # A list to store values over time
                  'fhandle': plt.figure()  # A figure handle }}

    Args:
        problem (pygeosx.Group): GEOSX problem handle
        records (dict): A dict of dicts that specifies the queries to run
    """
    for k in records.keys():
        if k == "time":
            current_time = get_wrapper(problem, "Events/time")
            records[k]["history"].append(current_time * records[k]["scale"])
        else:
            tmp = print_global_value_range(problem, k, records[k]["label"], scale=records[k]["scale"])
            records[k]["history"].append(tmp)
    sys.stdout.flush()


def plot_history(records, output_root=".", save_figures=True, show_figures=True):
    """
    Plot the time-histories for the records structure.
    Note: If figures are shown, the GEOSX process will be blocked until they are closed

    Args:
        records (dict): A dict of dicts containing the queries
        output_root (str): Path to save figures (default = './')
        save_figures (bool): Flag to indicate whether figures should be saved (default = True)
        show_figures (bool): Flag to indicate whether figures should be drawn (default = False)
    """
    if rank == 0:
        for k in records.keys():
            if k != "time":
                # Set the active figure
                fa = plt.figure(records[k]["fhandle"].number)

                # Assemble values to plot
                t = np.array(records["time"]["history"])
                x = np.array(records[k]["history"])
                N = np.shape(x)  # (time, min/max, dimensions)

                # Add plots
                if len(N) == 2:
                    # This is a 1D field
                    plt.gca().cla()
                    plt.plot(t, x[:, 0], label="min")
                    plt.plot(t, x[:, 1], label="max")
                    plt.xlabel(records["time"]["label"])
                    plt.ylabel(records[k]["label"])

                else:
                    # This is a 2D field
                    columns = 2
                    rows = int(np.ceil(N[2] / float(columns)))

                    # Setup axes
                    if ("axes" not in records[k]) or (len(fa.axes) == 0):
                        records[k]["axes"] = [plt.subplot(rows, columns, ii + 1) for ii in range(0, N[2])]

                    for ii in range(0, N[2]):
                        ax = records[k]["axes"][ii]
                        ax.cla()
                        ax.plot(t, x[:, 0, ii], label="min")
                        ax.plot(t, x[:, 1, ii], label="max")
                        ax.set_xlabel(records["time"]["label"])
                        ax.set_ylabel("%s (dim=%i)" % (records[k]["label"], ii))
                plt.legend(loc=2)
                records[k]["fhandle"].tight_layout(pad=1.5)

                if save_figures:
                    fname = k[k.rfind("/") + 1 :]
                    plt.savefig("%s/%s.png" % (output_root, fname), format="png")
        if show_figures:
            plt.show()
