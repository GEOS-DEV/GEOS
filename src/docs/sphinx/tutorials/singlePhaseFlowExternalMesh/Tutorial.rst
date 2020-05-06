.. _TutorialSinglePhaseFlowExternalMesh:

#########################################
Using an external mesh
#########################################



**Context**

In this tutorial, we use a simple single-phase flow solver (see :ref:`SinglePhaseFlow`)
from GEOSX to solve for pressure propagation on a mesh that we import into GEOSX.
The main goal of this tutorial is to work on importing external meshes,
an important feature to use GEOSX on meshes representing realistic models.

**Objectives**

At the end of this tutorial you will know:

  - the syntax and format of input meshes,
  - how to input external files into a GEOSX input XML file,
  - how to use and visualize hexahedral and tetrahedral meshes.


**Input file**

This tutorial uses external input files. You will need a GEOSX XML file and a mesh file.
The xml input file for this test case is located at:

.. code-block:: console

  src/coreComponents/physicsSolvers/fluidFlow/integratedTests/singlePhaseFlow/3D_10x10x10_compressible.xml



.. _ExternalHexahedral:
.. _2_ImportingExternalMesh:

All results are written in a format compatible with `VisIt
<https://wci.llnl.gov/simulation/computer-codes/visit/>`_.
