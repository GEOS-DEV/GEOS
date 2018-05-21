###############################################################################
GEOSX Getting Started Guide
###############################################################################

Getting Ready
=================================
GEOSX resides in a git repository hosted at https://github.com/GEOSX/GEOSX.
It is suggested that you setup ssh keys, and use ssh for your clones as discussed `here <https://help.github.com/articles/adding-a-new-ssh-key-to-your-github-account/>`_.
It is suggested that a directory be created to host the various clones that one may require for an effective development workflow, so the first step is to make a directory

1. mkdir geosx; cd geosx

Downloading the Code
=================================
There are currently two separate repositories that should be downloaded.
The first is the main repository, which may be cloned and initialized by the following steps: 

2. git clone git@github.com:GEOSX/GEOSX.git
   
   cd GEOSX
   
   git submodule init
   
   git submodule update
   
   cd ..

next, you must clone the third-party libraries, which may be cloned and initialized by the following steps: 

3. git clone git@github.com:GEOSX/thirdPartyLibs.git

   cd thirdPartyLibs
   
   git submodule init
   
   git submodule update
   
   cd ..


Compiling the Code
=================================

GEOSX compilations are typically driven by a hostconfig file, which reside in GEOSX/host-configs.
If your platform does not have a host-config in the repository, you are encouraged to maintain one.

The first step in compiling GEOSX is to run cmake and generate the makefiles.
Starting with the third-party libraries, the config-build.script will run cmake for you.
Note that the 'make' step should be run serially, as the indiviudal package builds are run in parallel by default.

4. cd thirdPartyLibs

   python scripts/config-build.py -hc ../GEOSX/host-configs/your-platform.cmake -bt Release
   
   cd build-your-platform-release
   
   make -j1

The next step is to compile the main code. 
Again, the config-build sets up cmake for you.

5. cd geosx/GEOSX

   python scripts/config-build.py -hc host-configs/your-platform.cmake
   
   cd build-your-platform-release
   
   make -j4
