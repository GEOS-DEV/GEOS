from distutils.core import setup

setup(
    name="hdf5_wrapper",
    version="0.1.0",
    description="Simple wrapper for h5py objects",
    author="Chris Sherman",
    author_email="sherman27@llnl.gov",
    packages=["hdf5_wrapper"],
    install_requires=["h5py>=2.10.0", "numpy>=1.16.2"],
)
