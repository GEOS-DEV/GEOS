.. _PVTDriver:

PVT Driver
===============

.. contents:: Table of Contents
    :depth: 3

Introduction
------------

When calibrating fluid material parameters to experimental or other reference data, it can be a hassle to launch a full flow simulation just to confirm density, viscosity, and other fluid properties are behaving as expected.  
Instead, GEOS provides a ``PVTDriver`` allowing the user to test fluid property models for a well defined set of pressure, temperature, and composition conditions.  
The driver itself is launched like any other GEOS simulation, but with a particular XML structure:

.. code-block:: sh

   ./bin/geosx -i myFluidTest.xml

This driver will work for any multi-phase fluid model (e.g. black-oil, co2-brine, compositional multiphase) enabled within GEOS.    

XML Structure
-------------
A typical XML file to run the driver will have several key elements. 
Here, we will walk through an example file included in the source tree at

.. code-block:: sh

   inputFiles/constitutiveDriver/testPVT_docExample.xml

The first thing to note is that the XML file structure is identical to a standard GEOS input deck.  
In fact, once the constitutive block is calibrated, one could start adding solver and discretization blocks to the same file to create a proper field simulation.  
This makes it easy to go back and forth between calibration and simulation.

The first step is to define a parameterized fluid model to test.
Here, we create a particular type of CO2-Brine mixture:

.. literalinclude:: ../../../../inputFiles/constitutiveDriver/testPVT_docExample.xml
  :language: xml
  :start-after: <!-- SPHINX_PVTDRIVER_CONSTITUTIVE_START --> 
  :end-before: <!-- SPHINX_PVTDRIVER_CONSTITUTIVE_END -->

We also define two time-history functions for the pressure (Pascal units) and temperature (Kelvin units) conditions we want to explore.

.. literalinclude:: ../../../../inputFiles/constitutiveDriver/testPVT_docExample.xml
  :language: xml
  :start-after: <!-- SPHINX_PVTDRIVER_FUNCTIONS_START --> 
  :end-before: <!-- SPHINX_PVTDRIVER_FUNCTIONS_END -->

Note that the time-axis here is just a pseudo-time, allowing us to `parameterize <https://en.wikipedia.org/wiki/Parametric_equation>`_ arbitrarily complicated paths through a (pressure,temperature) diagram.
The actual time values have no impact on the resulting fluid properties.
Here, we fix the temperature at 350K and simply ramp up pressure from 1 MPa to 50 MPa:

A ``PVTDriver`` is then added as a ``Task``, a particular type of executable event often used for simple actions. 

.. literalinclude:: ../../../../inputFiles/constitutiveDriver/testPVT_docExample.xml
  :language: xml
  :start-after: <!-- SPHINX_PVTDRIVER_TASKS_START --> 
  :end-before: <!-- SPHINX_PVTDRIVER_TASKS_END -->

The driver itself takes as input the fluid model, the pressure and temperature control functions, and a "feed composition". 
The latter is the mole fraction of each component in the mixture to be tested.
The ``steps`` parameter controls how many steps are taken along the parametric (P,T) path.
Results will be written in a simple ASCII table format (described below) to the file ``output``. If the ``output`` is not specified, output will be written to the standard log (screen).
The ``logLevel`` parameter controls the verbosity of log output during execution. 
The ``precision`` parameter determines the precision used to write the output data.
 
The driver task is added as a ``SoloEvent`` to the event queue.  
This leads to a trivial event queue, since all we do is launch the driver and then quit.

.. literalinclude:: ../../../../inputFiles/constitutiveDriver/testPVT_docExample.xml
  :language: xml
  :start-after: <!-- SPHINX_PVTDRIVER_EVENTS_START --> 
  :end-before: <!-- SPHINX_PVTDRIVER_EVENTS_END -->

Internally, the driver uses a simple form of time-stepping to advance through the (P,T) steps. 
This timestepping is handled independently of the more complicated time-stepping pattern used by physics ``Solvers`` and coordinated by the ``EventManager``.  
In particular, in the XML file above, the ``maxTime`` parameter in the ``Events`` block is an event manager control, controlling when/if certain events occur.  
Once launched, the PVTDriver internally determines its own max time and timestep size using a combination of the input functions' time coordinates and the requested number of loadsteps.  
It is therefore helpful to think of the driver as an instantaneous *event* (from the event manager's point of view), but one which has a separate, internal clock.

Parameters
----------
The key XML parameters for the PVTDriver are summarized in the following table:

.. include:: /docs/sphinx/datastructure/PVTDriver.rst

Output Format
-------------
The ``output`` key is used to identify a file to which the results of the simulation are written.  
If this key is omitted, file output will be suppressed and instead the resulting table will be output to the screen.  

When written to standard output, the data is written in a table format similar to the one below.

.. code:: sh

   -----------------------------------------------------------------------------------------------------------------------------------------------
   |                                                             Output for testCO2                                                              |
   |---------------------------------------------------------------------------------------------------------------------------------------------|
   |    time     |  pressure   |  temperature  |   density   |      phase fraction       |       phase density       |      phase viscosity      |
   |-------------|-------------|---------------|-------------|---------------------------|---------------------------|---------------------------|
   |             |             |               |             |     gas     |    water    |     gas     |    water    |     gas     |    water    |
   |-------------|-------------|---------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|
   |   0.00e+00  |   1.00e+06  |     3.50e+02  |   1.56e+01  |   1.00e+00  |   0.00e+00  |   1.56e+01  |   1.00e+03  |   1.75e-05  |   4.13e-04  |
   |   2.00e-02  |   1.98e+06  |     3.50e+02  |   3.18e+01  |   1.00e+00  |   0.00e+00  |   3.18e+01  |   1.00e+03  |   1.76e-05  |   4.13e-04  |
   |   4.00e-02  |   2.96e+06  |     3.50e+02  |   4.92e+01  |   1.00e+00  |   0.00e+00  |   4.92e+01  |   1.00e+03  |   1.78e-05  |   4.13e-04  |
   |   6.00e-02  |   3.94e+06  |     3.50e+02  |   6.78e+01  |   1.00e+00  |   0.00e+00  |   6.78e+01  |   1.00e+03  |   1.80e-05  |   4.13e-04  |
   ...
   |   9.80e-01  |   4.90e+07  |     3.50e+02  |   8.80e+02  |   1.00e+00  |   0.00e+00  |   8.80e+02  |   1.03e+03  |   8.74e-05  |   4.13e-04  |
   |   1.00e+00  |   5.00e+07  |     3.50e+02  |   8.85e+02  |   1.00e+00  |   0.00e+00  |   8.85e+02  |   1.03e+03  |   8.84e-05  |   4.13e-04  |
   -----------------------------------------------------------------------------------------------------------------------------------------------

When written to a file, the file is a simple ASCII format with a brief header followed by test data:

.. code:: sh

   # column 1 = time
   # column 2 = pressure
   # column 3 = temperature
   # column 4 = density
   # column 5 = phase fraction,gas
   # column 6 = phase fraction,water
   # column 7 = phase density,gas
   # column 8 = phase density,water
   # column 9 = phase viscosity,gas
   # column 10 = phase viscosity,water
    0.0000e+00  1.0000e+06  3.5000e+02  1.5581e+01  1.0000e+00  0.0000e+00  1.5581e+01  1.0022e+03  1.7476e-05  4.1330e-04 
    2.0000e-02  1.9800e+06  3.5000e+02  3.1833e+01  1.0000e+00  0.0000e+00  3.1833e+01  1.0028e+03  1.7598e-05  4.1330e-04 
    4.0000e-02  2.9600e+06  3.5000e+02  4.9192e+01  1.0000e+00  0.0000e+00  4.9192e+01  1.0034e+03  1.7771e-05  4.1330e-04 
    6.0000e-02  3.9400e+06  3.5000e+02  6.7831e+01  1.0000e+00  0.0000e+00  6.7831e+01  1.0040e+03  1.8004e-05  4.1330e-04 
   ...

Note that the number of columns will depend on how many phases and components are present and on whether the fluid is thermal or not.
In this case, we have a two-phase, two-component isothermal mixture.
The total density is reported in column 4, while phase fractions, phase densities, and phase viscosities are reported in subsequent columns.

If the ``outputCompressibility`` flag is activated, an extra column will be added for the total fluid compressibility after the density.
This is defined as :math:`c_t=\frac{1}{\rho_t}\left(\partial{\rho_t}/\partial P\right)` where :math:`\rho_t` is the total density.
If the ``outputMassDensity`` flag is activated, extra columns will be added for the mass density of each phase.
The number of columns will also depend on whether the ``outputPhaseComposition`` flag is activated or not. If it is activated, there will be an extra column for the mole fraction of each component in each phase.
The phase order will match the one defined in the input XML (here, the co2-rich phase followed by the water-rich phase).
This file can be readily plotted using any number of plotting tools.  Each row corresponds to one timestep of the driver, starting from initial conditions in the first row.

