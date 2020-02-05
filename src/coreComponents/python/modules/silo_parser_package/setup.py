import os
from distutils.core import setup
from distutils.extension import Extension


def parse_environment_variables():
  options = {'geosx_build': '',
             'silo': '/usr/gapps/GEOSX/thirdPartyLibs/2019-12-09/install-quartz-clang@9.0.0-release/silo',
             'hdf5': '/usr/gapps/GEOSX/thirdPartyLibs/2019-12-09/install-quartz-clang@9.0.0-release/hdf5',
             'mpi': '/usr/tce/packages/mvapich2/mvapich2-2.3-clang-9.0.0',
             'libz': '/usr/lib64/',
             'libdl': '/usr/lib64/',
             'libm': '/usr/lib64/'}

  # Override library paths using os environment variables
  for k in options.keys():
    if k in os.environ.keys():
      options[k] = os.environ[k]

  # If specified, find these using a cmake cache file in a GEOSX build dir
  if (options['geosx_build']):
    with open(os.path.expanduser('%s/CMakeCache.txt' % (options['geosx_build']))) as f:
      for line in f:
        if ('SILO_INCLUDE_DIRS:PATH' in line):
          tmp = line.split('=')[1]
          options['silo'] = tmp[:tmp.rfind('/')]
        elif ('BLT_HDF5_INCLUDES:STRING' in line):
          tmp = line.split('=')[1]
          options['hdf5'] = tmp[:tmp.rfind('/')]
        elif ('BLT_MPI_INCLUDES:STRING' in line):
          tmp = line.split('=')[1]
          options['mpi'] = tmp[:tmp.rfind('/')]

  # Assemble setup arguments
  extra_includes = ['%s/include' % (options['silo']),
                    '%s/include' % (options['hdf5'])]

  extra_libs = ['%s/lib' % (options['silo']),
                '%s/lib' % (options['hdf5'])]

  extra_files = ['%s/lib/libsiloh5.a' % (options['silo']),
                 '%s/lib/libhdf5.a' % (options['hdf5']),
                 '%s/lib/libmpi.so' % (options['mpi']),
                 '%s/lib/libmpicxx.so' % (options['mpi']),
                 '%s/libz.so' % (options['libz']),
                 '%s/libdl.so' % (options['libdl']),
                 '%s/libm.so' % (options['libm'])]

  return extra_includes, extra_libs, extra_files


extra_includes, extra_libs, extra_files = parse_environment_variables()


setup(name='silo_parser_wrapper',
      version='0.0.1',
      python_requires='>=3.0',
      description='Simple parser for silo files',
      author='Chris Sherman',
      author_email='sherman27@llnl.gov',
      packages=['silo_parser'],
      entry_points={'console_scripts': ['silo_parser = silo_parser.__main__:main']},
      ext_modules=[Extension('silo_parser.simple_silo_parser',
                             sources=['./silo_parser/simple_silo_parser.cpp'],
                             include_dirs=extra_includes,
                             library_dirs=extra_libs,
                             extra_objects=extra_files)])

