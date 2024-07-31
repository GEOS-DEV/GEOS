.. _AppleMacOS:

Building Apple MacOS
==============

Install homebrew
----------------
.. code-block::

/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
(echo; echo 'eval "$(/opt/homebrew/bin/brew shellenv)"') >> /Users/settgast1/.zprofile
eval "$(/opt/homebrew/bin/brew shellenv)"

Install packages using homebrew
-------------------------------
.. code-block::

brew install bison cmake gfortran git-lfs open-mpi lapack python3 ninja m4
echo 'export PATH="/opt/homebrew/opt/bison/bin:$PATH"' >> ~/.zshrc
echo 'export PATH="/opt/homebrew/opt/m4/bin:$PATH"' >> ~/.zshrc
git lfs install

Clone GEOS
----------
.. code-block::

git clone git@github.com:GEOS-DEV/GEOS.git
cd GEOS
git submodule init
git submodule update
cd ..

Clone thirdPartyLibs
--------------------
.. code-block::

git clone git@github.com:GEOS-DEV/thirdPartyLibs.git
cd thirdPartyLibs
git submodule init 
git submodule update
git lfs pull


Configure and build thirdPartyLibs
---------------------------------
.. code-block::

python3 scripts/config-build.py -hc ../GEOS/host-configs/apple/macOS_arm.cmake -bt Release

You will get a warning you can ignore
.. code-block::

CMake Warning at /Users/settgast1/Codes/geos/GEOS/host-configs/tpls.cmake:10 (message):
  'GEOS_TPL_DIR' does not exist.


.. code-block::
  cd build-macOS_arm-release
  make

  You will get an error at the end, ignore
  [100%] Linking CXX executable ../../../tests/blt_mpi_smoke
  ld: warning: -commons use_dylibs is no longer supported, using error treatment instead
  ld: file not found: @rpath/libquadmath.0.dylib for architecture arm64
  clang: error: linker command failed with exit code 1 (use -v to see invocation)
  make[2]: *** [tests/blt_mpi_smoke] Error 1
  make[1]: *** [blt/tests/smoke/CMakeFiles/blt_mpi_smoke.dir/all] Error 2
  make: *** [all] Error 2


build GEOS
.. code-block::

  cd ../../GEOS
  python3 scripts/config-build.py -hc host-configs/apple/macOS_arm.cmake -bt Release --ninja
  cd build-macOS_arm-release
  ninja geosx
