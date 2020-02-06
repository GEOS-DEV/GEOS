

from silo_parser import simple_silo_parser
from functools import partial
import multiprocessing
import glob
import signal
import os
import hdf5_wrapper
import numpy as np


def sig_handler(signum, frame):
  print('Processing error!')


signal.signal(signal.SIGSEGV, sig_handler)


def parse_single(field_type, field_names, parallel_folder, f):
    # Parse the silo file
    print('  %s' % (f))
    data = simple_silo_parser.parse_file(f, field_type, field_names, parallel_folder=parallel_folder)

    # Write to an hdf5 database
    cycle = data['cycle']
    sub_name = './sub_hdf5/%s_%08d.hdf5' % (field_type, cycle)
    with hdf5_wrapper.hdf5_wrapper(sub_name, mode='w') as sub_database:
      for k in data.keys():
        sub_database[k] = np.squeeze(np.array(data[k]))


def parse_multiple(plot_root, field_type, field_names, parallel_folder='data', npar_proc=1):
  print('Processing files...')
  os.makedirs('sub_hdf5', exist_ok=True)
  fnames = sorted(glob.glob(plot_root + '*'))
  if (npar_proc > 1):
    pfunc = partial(parse_single, field_type, field_names, parallel_folder)
    pool = multiprocessing.Pool(processes=npar_proc)
    pool.map(pfunc, fnames)
    pool.close()
    pool.join()
  else:
    for f in fnames:
      parse_single(field_type, field_names, parallel_folder, f)

  print('Linking files...')
  with hdf5_wrapper.hdf5_wrapper('post_database.hdf5', mode='w') as database:
    # Link in sub-databases
    for k in glob.glob('./sub_hdf5/*.hdf5'):
      sub_name = k[k.rfind('/')+1:k.rfind('.')]
      # database.link(sub_name, os.path.abspath(k))
      database.link(sub_name, k)

    # Build the metadata
    build_metadata_hdf5(database)


def build_metadata_hdf5(data):
  # Build metadata
  metadata = {'cycles': [], 'times': []}
  keymap = {}

  # Find cycles, timesteps
  for ka in sorted(data.keys()):
    metadata['cycles'].append(data[ka]['cycle'])
    metadata['times'].append(data[ka]['time'])

    for kb in data[ka].keys():
      if kb not in keymap.keys():
        keymap[kb] = []
      keymap[kb].append(ka)

  # Write to the database
  data['cycles'] = np.array(sorted(np.unique(metadata['cycles'])))
  data['times'] = np.array(sorted(np.unique(metadata['times'])))
  data['keymap'] = keymap


