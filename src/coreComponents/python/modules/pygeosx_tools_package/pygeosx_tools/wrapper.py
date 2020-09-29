
import sys
import numpy as np
from mpi4py import MPI
import matplotlib.pyplot as plt


# Get the MPI rank
comm = MPI.COMM_WORLD
rank = comm.Get_rank()


def get_wrapper(problem, key, write_flag=False):
  local_values = problem.get_wrapper(key).value(write_flag)
  if not isinstance(local_values, np.ndarray):
    local_values = local_values.to_numpy()
  return local_values


def get_global_value_range(problem, key):
  """
  @brief get the range of a target value across all processes
  @param problem the GEOSX problem handle
  @param key the path of the target value
  @return global_min the minimum value of the target
  @return global_max the maximum value of the target
  """
  local_values = get_wrapper(problem, key)

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


def print_global_value_range(problem, key, header, scale=1.0, precision='%1.4f'):
  """
  @brief print the range of a target value across all processes
  @param problem the GEOSX problem handle
  @param key the path of the target value
  @param header the short name to print with the range
  @param scale multiply the range with this value before printing (default = 1.0)
  @param precision the format for printing the range (default = %1.4f)
  @return global_min the scaled minimum value of the target
  @return global_max the scaled maximum value of the target
  """
  global_min, global_max = get_global_value_range(problem, key)
  global_min *= scale
  global_max *= scale

  if (rank == 0):
    if isinstance(global_min, np.ndarray):
      min_str = ', '.join([precision % (x) for x in global_min])
      max_str = ', '.join([precision % (x) for x in global_max])
      print('%s: min=[%s], max=[%s]' % (header, min_str, max_str))
    else:
      min_str = precision % (global_min)
      max_str = precision % (global_max)
      print('%s: min=%s, max=%s' % (header, min_str, max_str))

  # Return a copy of the min/max in case we want to use them
  return global_min, global_max


def set_wrapper_to_value(problem, key, value):
  """
  @brief set the value of a wrapper
  @param problem the GEOSX problem handle
  @param key the path of the target wrapper
  @param value set the values in the wrapper to this number
  """
  local_values = get_wrapper(problem, key, write_flag=True)
  local_values[...] = value


def set_wrapper_with_function(problem, target_key, input_keys, fn, target_index=-1):
  """
  @brief set the value of a wrapper using a function
  @param problem the GEOSX problem handle
  @param target_key the path of the target wrapper
  @param input_keys string (or list of strings) containing the input path(s)
  @param fn the vectorized function used to calculate target values
  @param target_index apply the function output to this target index (default = all)
  """
  if isinstance(input_keys, str):
    input_keys = [input_keys]
  local_target = get_wrapper(problem, target_key, write_flag=True)
  local_inputs = [get_wrapper(problem, k) for k in input_keys]

  # Run the function, check the shape of outputs/target
  fn_output = fn(*local_inputs)
  N = np.shape(local_target)
  M = np.shape(fn_output)

  if (target_index < 0):
    if (N == M):
      # Function output, target shapes are the same
      local_target[...] = fn_output

    elif (len(M) == 1):
      # Apply the function output across each of the target dims
      local_target[...] = np.tile(np.expand_dims(fn_output, axis=1), (1, N[1]))

    else:
      raise Exception('Shape of function output %s is not compatible with target %s' % (str(M), str(N)))
  elif (len(M) == 1):
    if (len(N) == 2):
      # 2D target, with 1D output applied to a given index
      local_target[:, target_index] = fn_output

    else:
      # ND target, with 1D output tiled across intermediate indices
      expand_axes = tuple([ii for ii in range(1, len(N) - 1)])
      tile_axes = tuple([1] + [ii for ii in N[1:-1]])
      local_target[..., target_index] = np.tile(np.expand_dims(fn_output, axis=expand_axes), tile_axes)

  else:
    raise Exception('Shape of function output %s is not compatible with target %s (target axis=%i)' % (str(M), str(N), target_index))


def search_datastructure_wrappers_recursive(group, filters, matching_paths, level=0, group_path=[]):
  """
  @brief recursively search the group and its children for wrappers that match the filters
  @param group the current GEOSX group
  @param filters a list of strings
  @param matching_paths a list of matching values
  @details a successful match is identified if the wrapper path contains all of the
           strings in the filter.  Note: multiple matches are possible

           For example, if filters=['a', 'b', 'c'], the following would match
             - a/b/c
             - c/b/a
             - d/e/c/f/b/a/a
  """
  for wrapper in group.wrappers():
    wrapper_path = str(wrapper).split()[0]
    wrapper_test = group_path + [wrapper_path[wrapper_path.rfind('/')+1:]]
    if all(f in wrapper_test for f in filters):
      matching_paths.append('/'.join(wrapper_test))

  for sub_group in group.groups():
    sub_group_name = str(sub_group).split()[0].split('/')[-1]
    search_datastructure_wrappers_recursive(sub_group,
                                            filters,
                                            matching_paths,
                                            level=level+1,
                                            group_path=group_path+[sub_group_name])


def get_matching_wrapper_path(problem, filters):
  """
  @brief recursively search the GEOSX data structure for a wrapper that matches the filter
  @param problem the GEOSX problem handle
  @param filters a list of strings
  @returns matching_path the path to the matching wrapper
  @details a successful match is identified if the wrapper path contains all of the
           strings in the filter.  Note: multiple matches are possible.  This function
           will throw an error is 0 or >1 wrappers are found

           For example, if filters=['a', 'b', 'c'], the following would match
             - a/b/c
             - c/b/a
             - d/e/c/f/b/a/a
  """
  matching_paths = []
  search_datastructure_wrappers_recursive(problem, filters, matching_paths)

  if (len(matching_paths) == 1):
    if (rank == 0):
      print('Found matching wrapper: %s' % (matching_paths[0]))
    return matching_paths[0]

  else:
    if (rank == 0):
      print('Error occured while looking for wrappers:')
      print('Filters: [%s]' % (', '.join(filters)))
      print('Matching wrappers: [%s]' % (', '.join(matching_paths)))
    raise Exception('Search resulted in 0 or >1 wrappers mathching filters')


def run_queries(problem, records):
  """
  @brief query the current GEOSX datastructure
  @param problem the GEOSX problem handle
  @param records a dict of dicts that specifies the queries to run
         and holds the results.  The expected format is as follows:
         record = {'path/of/wrapper': {'label': 'aperture (m)', # A label to include with plots
                                       'scale': 1.0,  # Value to scale results by
                                       'history: [],  # A list to store values over time
                                       'fhandle': plt.figure()  # A figure handle
                                       }}
  @details for now, the only supported query is to find the min/max values of the target
  """
  for k in records.keys():
    if (k == 'time'):
      current_time = get_wrapper(problem, "Events/time")
      records[k]['history'].append(current_time * records[k]['scale'])
    else:
      tmp = print_global_value_range(problem, k, records[k]['label'], scale=records[k]['scale'])
      records[k]['history'].append(tmp)
  sys.stdout.flush()


def plot_history(records,
                 output_root='.',
                 save_figures=True,
                 show_figures=True):
  """
  @brief plot the time-histories for the records structure
  @param records a dict of dicts that specifies the queries to run
         and holds the results.  The expected format is as follows:
         record = {'path/of/wrapper': {'label': 'aperture (m)', # A label to include with plots
                                       'scale': 1.0,  # Value to scale results by
                                       'history: [],  # A list to store values over time
                                       'fhandle': plt.figure()  # A figure handle
                                       }}
  @param output_root the path to save figures (default = './')
  @param save_figures flag to indicate whether figures should be saved (default = True)
  @param show_figures flag to indicate whether figures should be drawn (default = False)
  @details Note: if figures are shown, the GEOSX process will be blocked until they are closed
  """
  if (rank == 0):
    for k in records.keys():
      if (k != 'time'):
        # Set the active figure
        fa = plt.figure(records[k]['fhandle'].number)

        # Assemble values to plot
        t = np.array(records['time']['history'])
        x = np.array(records[k]['history'])
        N = np.shape(x)   # (time, min/max, dimensions)

        # Add plots
        if (len(N) == 2):
          # This is a 1D field
          plt.gca().cla()
          plt.plot(t, x[:, 0], label='min')
          plt.plot(t, x[:, 1], label='max')
          plt.xlabel(records['time']['label'])
          plt.ylabel(records[k]['label'])

        else:
          # This is a 2D field
          columns = 2
          rows = int(np.ceil(N[2] / float(columns)))

          # Setup axes
          if (('axes' not in records[k]) or (len(fa.axes) == 0)):
            records[k]['axes'] = [plt.subplot(rows, columns, ii+1) for ii in range(0, N[2])]

          for ii in range(0, N[2]):
            ax = records[k]['axes'][ii]
            ax.cla()
            ax.plot(t, x[:, 0, ii], label='min')
            ax.plot(t, x[:, 1, ii], label='max')
            ax.set_xlabel(records['time']['label'])
            ax.set_ylabel('%s (dim=%i)' % (records[k]['label'], ii))
        plt.legend(loc=2)
        records[k]['fhandle'].tight_layout(pad=1.5)

        if save_figures:
          fname = k[k.rfind('/')+1:]
          plt.savefig('%s/%s.png' % (output_root, fname), format='png')
    if show_figures:
      plt.show()

