.. _AppleMacOS:

Building on Apple macOS
=======================

The overall process is the one described in the :ref:`QuickStart`. This page lists what is specific to macOS
(tested on Apple Silicon with Homebrew).

Install the prerequisites with Homebrew
---------------------------------------

Install `Homebrew <https://brew.sh>`_ if needed:

.. code-block:: sh

  /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

At the end of the installation, Homebrew prints the command that adds it to your shell environment
(for ``zsh``, the default shell, it appends ``eval "$(/opt/homebrew/bin/brew shellenv)"`` to ``~/.zprofile``). Run it.

Then install the build dependencies. ``bison`` and ``m4`` must be the Homebrew versions, ahead of the system ones in your ``PATH``:

.. code-block:: sh

  brew install bison cmake gfortran open-mpi lapack python3 ninja m4
  echo 'export PATH="/opt/homebrew/opt/bison/bin:$PATH"' >> ~/.zshrc
  echo 'export PATH="/opt/homebrew/opt/m4/bin:$PATH"' >> ~/.zshrc

The C/C++ compiler is Apple Clang from the Xcode command line tools (``xcode-select --install``).
OpenMP is disabled with Apple Clang.

Host-config
-----------

Use ``host-configs/apple/macOS_arm.cmake`` (Apple Silicon) or ``host-configs/apple/macOS_x86.cmake`` (Intel).
Both include ``macOS_base.cmake``, which assumes Homebrew is installed in ``/opt/homebrew``.
You do not need to write your own host-config.

Build
-----

Clone the repositories as described in :ref:`QuickStartDownload`, then:

.. code-block:: sh

  cd thirdPartyLibs
  python3 scripts/config-build.py -hc ../GEOS/host-configs/apple/macOS_arm.cmake -bt Release
  cd build-macOS_arm-release
  make

  cd ../../GEOS
  python3 scripts/config-build.py -hc host-configs/apple/macOS_arm.cmake -bt Release --ninja
  cd build-macOS_arm-release
  ninja install

The executable is ``GEOS/install-macOS_arm-release/bin/geos``. Continue with :ref:`QuickStartFirstRun`.

Known messages during the TPL build
-----------------------------------

* While configuring the TPLs, CMake prints ``'GEOS_TPL_DIR' does not exist``.
  This is expected: the TPLs have not been installed yet. The same message when configuring **GEOS** means the
  TPLs were not found (check that ``thirdPartyLibs/install-macOS_arm-release`` exists).

* The TPL build may end with a link error in BLT's MPI smoke test:

  .. code-block:: none

    ld: file not found: @rpath/libquadmath.0.dylib for architecture arm64
    make[2]: *** [tests/blt_mpi_smoke] Error 1

  This test belongs to the BLT build system and is not used by GEOS. The error can be ignored as long as
  the libraries (``hypre``, ``vtk``, ``raja``, ...) are present in ``thirdPartyLibs/install-macOS_arm-release``;
  if some are missing, run ``make`` again to finish the build.
