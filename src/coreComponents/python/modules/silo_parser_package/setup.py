from distutils.core import setup
from distutils.extension import Extension

silo_root = '/usr/gapps/GEOSX/thirdPartyLibs/2019-12-09/install-quartz-clang@9.0.0-release/silo'
hdf5_root = '/usr/gapps/GEOSX/thirdPartyLibs/2019-12-09/install-quartz-clang@9.0.0-release/hdf5'
mpi_root = '/usr/tce/packages/mvapich2/mvapich2-2.3-clang-9.0.0'
cmake_cache = '/usr/WS2/sherman/GEOSX/build-quartz-clang@9.0.0-release/CMakeCache.txt'

# setup(name='silo_parser_wrapper',
#       version='0.0.1',
#       description='Simple parser for silo files',
#       author='Chris Sherman',
#       author_email='sherman27@llnl.gov',
#       packages=['silo_parser'],
#       entry_points={'console_scripts': ['silo_parser = silo_parser.__main__:main']},
#       ext_modules=[Extension('silo_parser.simple_silo_parser',
#                              sources=['./silo_parser/simple_silo_parser.cpp'],
#                              extra_compile_args=['-I%s/include' % (silo_root)],
#                              library_dirs=['%s/lib' % (silo_root),
#                                            '%s/lib' % (hdf5_root)],
#                              runtime_library_dirs=['%s/lib' % (silo_root),
#                                                    '%s/lib' % (hdf5_root)],
#                              libraries=['siloh5', 'hdf5'])])


extra_includes = ['%s/include' % (silo_root),
                  '%s/include' % (hdf5_root)]

extra_libs = ['%s/lib' % (silo_root),
              '%s/lib' % (hdf5_root)]

extra_files = ['%s/lib/libsiloh5.a' % (silo_root),
               '%s/lib/libhdf5.a' % (hdf5_root),
               '/usr/lib64/libz.so',
               '/usr/lib64/libdl.so',
               '/usr/lib64/libm.so',
               '%s/lib/libmpi.so' % (mpi_root),
               '%s/lib/libmpicxx.so' % (mpi_root)]


setup(name='silo_parser_wrapper',
      version='0.0.1',
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

