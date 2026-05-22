Relperm Driver
===============

.. contents:: Table of Contents
    :depth: 3

Introduction
------------

Fitting relative permeability parameters to experimental or benchmark data often should not require a full flow simulation, when the only goal is to check that the permeability curves respond correctly to changes in phase saturation.  
As such, GEOS provides a ``RelpermDriver`` allowing the user to test relative permeability models for a well-defined sweep of saturation conditions.  
The driver itself is launched like any other GEOS simulation, but with a particular XML structure:

.. code-block:: sh

   ./bin/geosx -i myRelpermTest.xml

This driver will work for any 2-phase or 3-phase relative permeability model enabled within GEOS.    

XML Structure
-------------
A typical XML file to run the driver will have several key elements. 
Here, we will walk through an example file included in the source tree at

.. code-block:: sh

   inputFiles/constitutiveDriver/testRelperm_docExample.xml

A key point is that the XML file follows the same structure as a typical GEOS input deck. As with the other constitutive drivers, once the constitutive block has been calibrated, solver and discretization sections can be added directly to transform it into a full field simulation setup. This continuity allows for smooth transitions between calibration and simulation workflows.

The first step is to define a parameterized relative permeability model to test.
Here, we create a particular type of table-based two-phase relative permeability for sandstone:

.. literalinclude:: ../../../../inputFiles/constitutiveDriver/testRelperm_docExample.xml
  :language: xml
  :start-after: <!-- SPHINX_RELPERMDRIVER_CONSTITUTIVE_START --> 
  :end-before: <!-- SPHINX_RELPERMDRIVER_CONSTITUTIVE_END -->

We also define the underlying table functions used by the model.

.. literalinclude:: ../../../../inputFiles/constitutiveDriver/testRelperm_docExample.xml
  :language: xml
  :start-after: <!-- SPHINX_RELPERMDRIVER_FUNCTIONS_START --> 
  :end-before: <!-- SPHINX_RELPERMDRIVER_FUNCTIONS_END -->

A ``RelpermDriver`` is then added as a ``Task``, a particular type of executable event often used for simple actions. 

.. literalinclude:: ../../../../inputFiles/constitutiveDriver/testRelperm_docExample.xml
  :language: xml
  :start-after: <!-- SPHINX_RELPERMDRIVER_TASKS_START --> 
  :end-before: <!-- SPHINX_RELPERMDRIVER_TASKS_END -->

The driver itself takes as input the relative permeability model (referenced by the ``relperm`` parameter) and the ``steps`` parameter, which controls how many saturation increments are evaluated.
Results will be written in a simple ASCII table format (described below) to a specified file if an ``output`` parameter is provided. If it is not specified, output will be written to the standard log (screen).
The ``logLevel`` parameter controls the verbosity of log output during execution. 
 
The driver task is added as a ``SoloEvent`` to the event queue.  
This leads to a trivial event queue, since all we do is launch the driver and then quit.

.. literalinclude:: ../../../../inputFiles/constitutiveDriver/testRelperm_docExample.xml
  :language: xml
  :start-after: <!-- SPHINX_RELPERMDRIVER_EVENT_START --> 
  :end-before: <!-- SPHINX_RELPERMDRIVER_EVENT_END -->

Parameters
----------
The key XML parameters for the RelpermDriver are summarized in the following table:

.. include:: /docs/sphinx/datastructure/RelpermDriver.rst

Phase Count and Step Logic
--------------------------
The internal behavior of the driver depends on the number of phases present in the relative permeability model:

* **2-Phase Models:** The driver iterates exactly ``steps`` + 1 times, spanning from the minimum non-wetting saturation to its maximum valid bounds.
* **3-Phase Models:** The driver iterates through combinations of the wetting and non-wetting phase saturations based on the specified number of ``steps``. Only combinations where the sum of saturations is valid (i.e., less than or equal to 1.0) are evaluated, with the third intermediate phase making up the balance. Thus, the total number of output rows for a 3-phase model depends heavily on the saturation end-points and will be generally greater than the provided ``steps`` parameter.

Hysteresis Support
------------------
For relative permeability models featuring hysteresis (e.g., ``TableRelativePermeabilityHysteresis``), the driver requires historical phase saturations to determine the departure paths of the curves. This can be directly specified via the ``historicalSaturations`` array in the XML. If omitted, the driver will default to extracting the extremum phase volume fractions from the underlying curve definitions.

Output Format
-------------
The ``output`` key is used to identify a file to which the results of the simulation are written.  
If this key is omitted, file output will be suppressed and instead the resulting table will be output to the screen.  

When written to standard output, the data is written in a table format similar to the one below.

.. code:: sh

   ---------------------------------------------------------------------------------
   |                           Output for test_sandstone                           |
   |-------------------------------------------------------------------------------|
   |     index     |          saturation           |            relperm            |
   |---------------|-------------------------------|-------------------------------|
   |               |      gas      |     water     |      gas      |     water     |
   |---------------|---------------|---------------|---------------|---------------|
   |   0.0000e+00  |   0.0000e+00  |   1.0000e+00  |   0.0000e+00  |   1.0000e+00  |
   |   1.0000e+00  |   3.9525e-02  |   9.6047e-01  |   2.3250e-03  |   8.2337e-01  |
   |   2.0000e+00  |   7.9050e-02  |   9.2095e-01  |   7.0157e-03  |   6.6543e-01  |
   |   3.0000e+00  |   1.1857e-01  |   8.8143e-01  |   1.4124e-02  |   5.2659e-01  |
   |   4.0000e+00  |   1.5810e-01  |   8.4190e-01  |   2.5831e-02  |   4.1943e-01  |
   |   5.0000e+00  |   1.9762e-01  |   8.0238e-01  |   4.0028e-02  |   3.2628e-01  |
   |   6.0000e+00  |   2.3715e-01  |   7.6285e-01  |   5.6826e-02  |   2.4769e-01  |
   |   7.0000e+00  |   2.7668e-01  |   7.2332e-01  |   7.7434e-02  |   1.8665e-01  |
   |   8.0000e+00  |   3.1620e-01  |   6.8380e-01  |   1.0046e-01  |   1.3539e-01  |
   |   9.0000e+00  |   3.5572e-01  |   6.4428e-01  |   1.2703e-01  |   9.6594e-02  |
   |   1.0000e+01  |   3.9525e-01  |   6.0475e-01  |   1.5710e-01  |   6.7650e-02  |
   |   1.1000e+01  |   4.3477e-01  |   5.6523e-01  |   1.8949e-01  |   4.4341e-02  |
   |   1.2000e+01  |   4.7430e-01  |   5.2570e-01  |   2.2543e-01  |   2.8331e-02  |
   |   1.3000e+01  |   5.1382e-01  |   4.8618e-01  |   2.6486e-01  |   1.7448e-02  |
   |   1.4000e+01  |   5.5335e-01  |   4.4665e-01  |   3.0663e-01  |   9.4157e-03  |
   |   1.5000e+01  |   5.9287e-01  |   4.0713e-01  |   3.5200e-01  |   4.8321e-03  |
   |   1.6000e+01  |   6.3240e-01  |   3.6760e-01  |   4.0105e-01  |   2.3675e-03  |
   |   1.7000e+01  |   6.7192e-01  |   3.2808e-01  |   4.5150e-01  |   6.9571e-04  |
   |   1.8000e+01  |   7.1145e-01  |   2.8855e-01  |   5.0696e-01  |   2.9714e-04  |
   |   1.9000e+01  |   7.5097e-01  |   2.4903e-01  |   5.6478e-01  |   6.6429e-05  |
   |   2.0000e+01  |   7.9050e-01  |   2.0950e-01  |   6.2490e-01  |   4.4409e-20  |
   ---------------------------------------------------------------------------------

When written to a file, the file is a simple ASCII format with a brief header followed by the actual saturation data:

.. code:: sh

   # column 1 = index
   # column 2 = saturation,gas
   # column 3 = saturation,water
   # column 4 = relperm,gas
   # column 5 = relperm,water
    0.0000e+00 0.0000e+00 1.0000e+00 0.0000e+00 1.0000e+00
    1.0000e+00 3.9525e-02 9.6047e-01 2.3250e-03 8.2337e-01
    2.0000e+00 7.9050e-02 9.2095e-01 7.0157e-03 6.6543e-01
    3.0000e+00 1.1857e-01 8.8143e-01 1.4124e-02 5.2659e-01
    4.0000e+00 1.5810e-01 8.4190e-01 2.5831e-02 4.1943e-01
    5.0000e+00 1.9762e-01 8.0238e-01 4.0028e-02 3.2628e-01
   ...

The exact number of columns will depend on the phase count configured in the chosen model. In this 2-phase example, there are 5 total columns tracking the step index, saturations, and resulting relative permeabilities. If hysteresis is activated on the model, additional columns tracking historical saturations will be present between the saturations and the calculated relative permeabilities. 
This file format can be easily parsed using standard tools or plotting scripts to examine the generated curves.
