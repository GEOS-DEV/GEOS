.. _TriaxialDriver:


Triaxial Driver
===============

.. contents:: Table of Contents
    :depth: 3

Introduction
------------

When calibrating solid material parameters to experimental data, it can be a hassle to launch a full finite element simulation to mimic experimental loading conditions. Instead, GEOS provides a ``TriaxialDriver`` allowing the user to run loading tests on a single material point. This makes it easy to understand the material response and fit it to lab data. The driver itself is launched like any other GEOS simulation, but with a particular XML structure:

.. code-block:: sh

   ./bin/geosx -i myTest.xml


XML Structure
-------------
A typical XML file to run the triaxial driver will have the following key elements. We present the whole file first, before digging into the individual blocks.

.. literalinclude:: ../../../../inputFiles/constitutiveDriver/testTriaxial_druckerPragerExtended.xml
  :language: xml

The first thing to note is that the XML structure is identical to a standard GEOS input deck. In fact, once the constitutive block is calibrated, one could start adding solver and discretization blocks to the same file to create a proper field simulation. This makes it easy to go back and forth between calibration and simulation.

The ``TriaxialDriver`` is added as a ``Task``, a particular type of executable event often used for simple actions. It is added as a ``SoloEvent`` to the event queue. This leads to a trivial event queue, since all we do is launch the driver and then quit.

The driver itself takes as input the name of the solid model ``material``.
The ``steps`` parameter controls how many steps are taken along the parametric path.
Results will be written in a simple ASCII table format (described below) to the file ``output``. If the ``output`` is not specified, output will be written to the standard log (screen).
The ``logLevel`` parameter controls the verbosity of log output during execution. 
The ``precision`` parameter determines the precision used to write the output data.
 
Internally, the triaxial driver uses a simple form of time-stepping to advance through the loading steps, allowing for both rate-dependent and rate-independent models to be tested. This timestepping is handled independently from the more complicated time-stepping pattern used by physics ``Solvers`` and coordinated by the ``EventManager``. In particular, in the XML file above, the ``maxTime`` parameter in the ``Events`` block is an event manager control, controlling when/if certain events occur. Once launched, the triaxial driver internally determines its own max time and timestep size using a combination of the strain function's time coordinates and the requested number of loadsteps. It is therefore helpful to think of the driver as an instantaneous *event* (from the event manager's point of view), but one which has a separate, internal clock.

The key parameters for the TriaxialDriver are:

.. include:: /docs/sphinx/datastructure/TriaxialDriver.rst

.. note::

   GEOS uses the *engineering* sign convention where compressive stresses and strains are *negative*.
   This is one of the most frequent issues users make when calibrating material parameters, as
   stress- and strain-like quantities often need to be negative to make physical sense. You may note in the
   XML above, for example, that ``stressFunction`` and ``strainFunction`` have negative values for
   a compressive test.

Test Modes
----------
The most complicated part of the driver is understanding how the stress and strain functions are applied in different testing modes. The driver mimics laboratory core tests, with loading controlled in the axial and radial directions. These conditions may be either strain-controlled or stress-controlled, with the user providing time-dependent functions to describe the loading. The following table describes the available test modes in detail:

+--------------------+-------------------------+--------------------------+---------------------------+
| **mode**           | **axial loading**       | **radial loading**       | **initial stress**        |
+--------------------+-------------------------+--------------------------+---------------------------+
| ``strainControl``  | axial strain controlled | radial strain controlled | isotropic stress using    |
|                    | with ``axialControl``   | with ``radialControl``   | ``initialStress``         |
+--------------------+-------------------------+--------------------------+---------------------------+
| ``stressControl``  | axial stress controlled | radial stress controlled | isotropic stress using    |
|                    | with ``axialControl``   | with ``radialControl``   | ``initialStress``         |
+--------------------+-------------------------+--------------------------+---------------------------+
| ``mixedControl``   | axial strain controlled | radial stress controlled | isotropic stress using    |
|                    | with ``axialControl``   | with ``radialControl``   | ``initialStress``         |
+--------------------+-------------------------+--------------------------+---------------------------+

Note that a classical triaxial test can be described using either the ``stressControl`` or ``mixedControl`` mode. We recommend using the ``mixedControl`` mode when possible, because this almost always leads to well-posed loading conditions. In a pure stress controlled test, it is possible for the user to request that the material sustain a load beyond its intrinsic strength envelope, in which case there is no feasible solution and the driver will fail to converge. Imagine, for example, a perfectly plastic material with a yield strength of 10 MPa, but the user attempts to load it to 11 MPa. 

A volumetric test can be created by setting the axial and radial control functions to the same time history function. Similarly, an oedometer test can be created by setting the radial strain to zero. 

The user should be careful to ensure that the initial stress set via the ``initialStress`` value is consistent any applied stresses set through axial or radial loading functions. Otherwise, the material may experience sudden and unexpected deformation at the first timestep because it is not in static equilibrium.

Output Format
-------------
The ``output`` key is used to identify a file to which the results of the simulation are written.  
If this key is omitted, file output will be suppressed and instead the resulting table will be output to the screen.  

When written to standard output, the data is written in a table format similar to the one below.

.. code:: sh

   ---------------------------------------------------------------------------------------------------------------------------------------------------
   |                                                            Output for triaxialDriver                                                            |
   |-------------------------------------------------------------------------------------------------------------------------------------------------|
   |     time      |                    strain                     |                    stress                     |  newton_iter  |  residual_norm  |
   |---------------|-----------------------------------------------|-----------------------------------------------|---------------|-----------------|
   |               |     axial     |   radial_1    |   radial_2    |     axial     |   radial_1    |   radial_2    |               |                 |
   |---------------|---------------|---------------|---------------|---------------|---------------|---------------|---------------|-----------------|
   |   0.0000e+00  |   0.0000e+00  |   0.0000e+00  |   0.0000e+00  |  -1.9600e+07  |  -1.9600e+07  |  -1.9600e+07  |   0.0000e+00  |     0.0000e+00  |
   |   2.0000e-02  |  -7.0000e-05  |   1.8837e-05  |   1.8837e-05  |  -2.0584e+07  |  -1.9600e+07  |  -1.9600e+07  |   1.0000e+00  |     0.0000e+00  |
   ...
   |   9.8000e-01  |  -3.4300e-03  |   8.6963e-04  |   8.6963e-04  |  -5.8735e+07  |  -1.9600e+07  |  -1.9600e+07  |   3.0000e+00  |     4.3034e-07  |
   |   1.0000e+00  |  -3.5000e-03  |   9.0207e-04  |   9.0207e-04  |  -5.8760e+07  |  -1.9600e+07  |  -1.9600e+07  |   3.0000e+00  |     5.9024e-07  |
   ---------------------------------------------------------------------------------------------------------------------------------------------------

When written to a file, the file is a simple ASCII format with a brief header followed by test data:

.. code:: sh

  # column 1 = time
  # column 2 = strain,axial
  # column 3 = strain,radial_1
  # column 4 = strain,radial_2
  # column 5 = stress,axial
  # column 6 = stress,radial_1
  # column 7 = stress,radial_2
  # column 8 = newton_iter
  # column 9 = residual_norm
   0.0000e+00 0.0000e+00 0.0000e+00 0.0000e+00-4.6000e+06-4.6000e+06-4.6000e+06 0.0000e+00 0.0000e+00
   2.0000e-02-8.0000e-05 2.1528e-05 2.1528e-05-5.7249e+06-4.6000e+06-4.6000e+06 1.0000e+00 0.0000e+00
   4.0000e-02-1.6000e-04 4.3056e-05 4.3056e-05-6.8499e+06-4.6000e+06-4.6000e+06 1.0000e+00 0.0000e+00
   6.0000e-02-2.4000e-04 6.4585e-05 6.4585e-05-7.9748e+06-4.6000e+06-4.6000e+06 1.0000e+00 0.0000e+00
   ...

This file can be readily plotted using any number of plotting tools. Each row corresponds to one timestep of the driver, starting from initial conditions in the first row.

We note that the file contains two columns for radial strain and two columns for radial stress. For an isotropic material, the stresses and strains along the two radial axes will usually be identical. We choose to output this way, however, to accommodate both anisotropic materials and true-triaxial loading conditions. In these cases, the stresses and strains in the radial directions could potentially differ.

These columns can be added and subtracted to produce other quantities of interest, like mean stress or deviatoric stress. For example, we can plot the output produce stress / strain curves (in this case for a plastic rather than simple elastic material):

.. figure:: TriaxialDriver.svg
   :width: 600px
   :align: center
   :alt: stress/strain figure
   :figclass: align-center

   Stress/strain behavior for a plastic material.

In this plot, we have reversed the sign convention to be consistent with typical experimental plots. Note also that the ``strainFunction`` includes two unloading cycles, allowing us to observe both plastic loading and elastic unloading.

Model Convergence
-----------------

The last two columns of the output file contain information about the convergence behavior of the material driver. In ``triaxial`` mode, the mixed nature of the stress/strain control requires using a Newton solver to converge the solution. This last column reports the number of Newton iterations and final residual norm. Large values here would be indicative of the material model struggling (or failing) to converge. Convergence failures can result from several reasons, including:

1. Inappropriate material parameter settings
2. Overly large timesteps
3. Infeasible loading conditions (i.e. trying to load a material to a physically-unreachable stress point)
4. Poor model implementation

We generally spend a lot of time vetting the material model implementations (#4). When you first encounter a problem, it is therefore good to explore the other three scenarios first. If you find something unusual in the model implementation or are just really stuck, please submit an issue on our issue tracker so we can help resolve any bugs.
