from distutils.core import setup

setup(name='time_history_plotting',
      version='0.1.0',
      description='Scripts to plot time-series data from GEOSX time-history output files.',
      author='William Tobin',
      author_email='tobin6@llnl.gov',
      packages=['plot_time_history'],
      install_requires=['matplotlib', 'hdf5_wrapper', 'h5py', 'numpy'])
