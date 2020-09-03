
import sys
import numpy as np
from mpi4py import MPI
import matplotlib.pyplot as plt
import pygeosx
from geosx_xml_tools import xml_processor


# Get the MPI rank
comm = MPI.COMM_WORLD
rank = comm.Get_rank()


def apply_xml_preprocessor():
  """
  @brief modify the input file before handing it to GEOSX
  """
  new_input_file = ''
  file_index = sys.argv.index('-i') + 1
  if (rank == 0):
    print('Applying xml preprocessor...')
    new_input_file = xml_processor.process(sys.argv[file_index], '%s.processed' % (sys.argv[file_index]))
    print('  the compiled filename is: %s' % (new_input_file))

  # Broadcast and set the new input file name
  sys.argv[file_index] = comm.bcast(new_input_file, root=0)


def get_global_value_range(problem, key):
  """
  @brief get the range of a target value across all processes
  @param problem the GEOSX problem handle
  @param key the path of the target value
  @return global_min the minimum value of the target
  @return global_max the maximum value of the target
  """
  local_values = problem.getWrapper(key).value(True).toNumPy()

  # 1D arrays will return a scalar, ND arrays an array
  N = np.shape(local_values)
  local_min = 1e100
  local_max = -1e100
  if (len(N) > 1):
    local_min = np.zeros(N[1]) + 1e100
    local_max = np.zeros(N[1]) - 1e100

  # Ignore zero-length results
  if len(local_values):
    local_min = np.amin(local_values, axis=0)
    local_max = np.amax(local_values, axis=0)

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


def search_datastructure_wrappers_recursive(group, filters, matching_paths):
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
    wrapper_path = str(wrapper).split()[0].split('/')
    if all(f in wrapper_path for f in filters):
      matching_paths.append('/'.join(wrapper_path[1:]))

  for sub_group in group.groups():
    search_datastructure_wrappers_recursive(sub_group, filters, matching_paths)


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
      current_time = problem.getWrapper("Events/time").value(False)
      records[k]['history'].append(current_time * records[k]['scale'])
    else:
      tmp = print_global_value_range(problem, k, records[k]['label'], scale=records[k]['scale'])
      records[k]['history'].append(tmp)
  sys.stdout.flush()


def plot_history(records,
                 output_root='.',
                 save_figures=True,
                 show_figures=False):
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
        plt.figure(records[k]['fhandle'].number)

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
          for ii in range(0, N[2]):
            ax = plt.subplot(rows, columns, ii+1)
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


def run_problem():
  """
  @brief Run the GEOSX problem
  """
  # Initialize the code
  apply_xml_preprocessor()
  problem = pygeosx.initialize(rank, sys.argv)

  # Apply initial conditions
  pygeosx.applyInitialConditions()

  # Rather than specifying the wrapper paths explicitly,
  # search for them using a set of filters
  fracture_location_key = get_matching_wrapper_path(problem, ['Fracture', 'elementCenter'])
  fracture_aperture_key = get_matching_wrapper_path(problem, ['Fracture', 'effectiveAperture'])

  # Note: we will be plotting our results.
  # This will modify the fonts so that they a bit easier to read
  plot_font = {'weight': 'normal', 'size': 8}
  plt.rc('font', **plot_font)

  # Setup values to record
  records = {fracture_location_key: {'label': 'Fracture Extents (m)',
                                     'scale': 1.0,
                                     'history': [],
                                     'fhandle': plt.figure()},
             fracture_aperture_key: {'label': 'Aperture (mm)',
                                     'scale': 1e3,
                                     'history': [],
                                     'fhandle': plt.figure()},
             'time': {'label': 'Time (min)',
                      'scale': 1.0 / 60.0,
                      'history': []}}

  # Run the code
  while pygeosx.run() != pygeosx.COMPLETED:
    run_queries(problem, records)
    plot_history(records)


if __name__ == '__main__':
  run_problem()


