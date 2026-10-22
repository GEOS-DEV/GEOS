.. _BuildGuide:

###############################################################################
Build Guide
###############################################################################

This guide is the reference for building GEOS and its third-party libraries (TPLs).
If you are building GEOS for the first time, follow the :ref:`QuickStart` instead, which walks through the recommended route
(build the TPLs with the ``thirdPartyLibs`` repository, then GEOS) and the faster Docker route.

* :ref:`Prerequisites`: compilers and system tools required.
* :ref:`Dependencies`: list of third-party libraries, their CMake variables, and how to build them.
* :ref:`BuildProcess`: ``config-build.py`` options and GEOS CMake options.
* :ref:`SpackUberenv`: building the TPLs with Spack (HPC systems, advanced).
* :ref:`Continuous_Integration_process`: Docker images used by the continuous integration.
* :ref:`AppleMacOS`: macOS-specific instructions.

.. toctree::
   :maxdepth: 1

   Prerequisites

   Dependencies

   BuildProcess

   SpackUberenv

   ContinuousIntegration

   AppleMacOS
