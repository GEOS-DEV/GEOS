Relperm Driver
===============

.. contents:: Table of Contents
    :depth: 3

Introduction
------------

Fitting relative permeability parameters to experimental or benchmark data often should not require a full flow simulation, when the only goal is to check that the relative permeability curves respond correctly to changes in phase saturation. As such, GEOS provides a ``RelpermDriver`` allowing the user to test relative permeability models for a well-defined sweep of saturation conditions. The driver itself is launched like any other GEOS simulation, but with a particular XML structure:

.. code-block:: sh

   ./bin/geosx -i myRelpermTest.xml

This driver will work for any 2-phase or 3-phase relative permeability model enabled within GEOS.

XML Structure
-------------
A typical XML file to run the driver will have several key elements. Here, we will walk through an example file included in the source tree at:

.. code-block:: sh

   inputFiles/constitutiveDriver/testRelperm_docExample.xml

A key point is that the XML file follows the same structure as a typical GEOS input deck. As with the other constitutive drivers, once the constitutive block has been calibrated, solver and discretization sections can be added directly to transform it into a full field simulation setup. This continuity allows for smooth transitions between calibration and simulation workflows.

The first step is to define a parameterized relative permeability model to test. Here, we create particular types of table-based two-phase relative permeabilities for sandstone, one standard and one with hysteresis:

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

The driver itself takes as input the relative permeability model (referenced by the ``relperm`` parameter). The model is evaluated over time increments defined by the ``steps`` parameter, while the saturations at each step are evaluated using functions provided in the ``saturationControls`` array. Results will be written in a simple ASCII table format to a specified file if an ``output`` parameter is provided; otherwise, it is written to the standard log. 

The driver task is added as a ``SoloEvent`` to the event queue.

.. literalinclude:: ../../../../inputFiles/constitutiveDriver/testRelperm_docExample.xml
  :language: xml
  :start-after: <!-- SPHINX_RELPERMDRIVER_EVENT_START -->
  :end-before: <!-- SPHINX_RELPERMDRIVER_EVENT_END -->

Parameters
----------
The key XML parameters for the RelpermDriver are summarized in the following table:

.. include:: /docs/sphinx/datastructure/RelpermDriver.rst

Saturation Assignment
---------------------
The phase saturations are explicitly driven by user-defined functions over time. 

To configure the driver, the user must specify a list of target phases via ``phaseNames`` and corresponding time-history functions via ``saturationControls``. For an $N$-phase model, exactly $N-1$ phases and saturation functions must be provided. The saturation of the remaining implicit phase is calculated automatically by the driver to ensure the phase volume fractions (saturations) always sum to exactly 1.0. 

During execution, the driver evaluates the provided functions at each step. To ensure physical validity:
* If a saturation function evaluates to a negative value at a given step, the saturation is clamped to zero.
* If the evaluated saturations sum to a value greater than 1.0, they are proportionally scaled down so that their collective sum equals 1.0, and the implicit phase is assigned a saturation of zero.

Examples and Hysteresis
-----------------------
The documentation example contains two distinct driver tests that showcase different saturation trajectories:

1. **Standard Relperm Test (``test_sandstone``):** This test evaluates a standard relative permeability model. The gas saturation is driven by the ``gasSaturation`` function, which simply increases linearly from 0.0 to 1.0 over the evaluation time.
2. **Hysteresis Relperm Test (``test_sandstone_hysteresis``):** This test evaluates a hysteresis-enabled model and requires a non-monotonic saturation path to test historical tracking. It uses the ``gasSaturationHysteresis`` function, which dynamically varies the gas saturation: rising from 0.0 to 0.4 (primary drainage), dropping to 0.3 (imbibition), rising again to 0.6 (secondary drainage), and finally dropping to 0.4 (secondary imbibition).

These saturation functions are explicitly defined in the XML as shown below:

.. literalinclude:: ../../../../inputFiles/constitutiveDriver/testRelperm_docExample.xml
  :language: xml
  :start-after: <!-- SPHINX_RELPERMDRIVER_SATURATION_FUNCTIONS_BEGIN -->
  :end-before: <!-- SPHINX_RELPERMDRIVER_SATURATION_FUNCTIONS_END -->

.. plot:: ../inputFiles/constitutiveDriver/testRelperm_data/plot_relperm.py

Output of the ``test_sandstone_hysteresis`` task. The plotted trajectory follows the defined ``gasSaturationHysteresis`` values ``{ 0.0, 0.4, 0.3, 0.6, 0.4 }``. The color scale indicates the internally tracked maximum historical gas saturation (Sghy), dictating the departure curves:
   
* **Primary Drainage (0.0 to 0.4):** Sghy updates alongside the current saturation (dark blue to yellow).
* **Primary Imbibition (0.4 to 0.3):** The curve drops from the main envelope; Sghy remains constant at 0.4 (yellow).
* **Secondary Drainage (0.3 to 0.6):** The curve climbs back up, with Sghy resuming updates past 0.4 until it hits 0.6 (dark red).
* **Secondary Imbibition (0.6 to 0.4):** The relative permeability drops along a new, lower path, while Sghy stays fixed at 0.6 (dark red).

Output Format
-------------
The ``output`` key is used to identify a file to which the results of the simulation are written. If this key is omitted, file output will be suppressed and instead the resulting table will be output to the screen. When written to standard output, the data is written in a table format similar to the one below.

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
   ---------------------------------------------------------------------------------

When written to a file, the file is a simple ASCII format with a brief header followed by the actual saturation data. An example is shown below.

.. code:: sh

   # column 1 = index
   # column 2 = saturation,gas
   # column 3 = saturation,water
   # column 4 = relperm,gas
   # column 5 = relperm,water
    0.0000e+00 1.0000e+00 0.0000e+00 1.0000e+00 0.0000e+00
    2.0000e-01 8.0000e-01 2.0000e-01 7.0006e-01 7.1429e-04
    4.0000e-01 6.0000e-01 4.0000e-01 4.4233e-01 1.5671e-02
    6.0000e-01 4.0000e-01 6.0000e-01 2.3163e-01 9.8057e-02
    8.0000e-01 2.0000e-01 8.0000e-01 7.6686e-02 3.6119e-01
    1.0000e+00 0.0000e+00 1.0000e+00 0.0000e+00 1.0000e+00

The exact number of columns will depend on the phase count configured in the chosen model. If hysteresis is activated on the model, additional columns tracking historical extremum saturations will automatically be present between the instantaneous saturations and the calculated relative permeabilities. This file can be readily plotted using any number of plotting tools.  Each row corresponds to one timestep of the driver, starting from initial conditions in the first row.

Exploring 3-phases data
------------------------

Another ability of ``RelpermDriver`` is to test 3-phases relperm interpolation models. GEOS offers an implementation for the ``Baker`` model and the ``StoneII`` model.
The example input is gathered here,

.. code-block:: sh

   inputFiles/constitutiveDriver/testRelperm_3_phase_example.xml

.. note::

  Note :math:`S_g = 1 - S_o - S_w` and hence is not an input of the ``RelpermDriver``


As this shares the same usage pattern, only the constitutive blocks is highlighted.

.. literalinclude:: ../../../../inputFiles/constitutiveDriver/testRelperm_3_phase_example.xml
  :language: xml
  :start-after: <!-- SPHINX_RELPERMDRIVER_CONSTITUTIVE_3P_BEGIN -->
  :end-before: <!-- SPHINX_RELPERMDRIVER_CONSTITUTIVE_3P_END -->


Eventually, to run it, the ``waterSaturation`` and ``oilSaturation`` controls are simply browsing the full range of accessible water and oil saturations. 
This results in a regular sampling that is further interpolated by the plotting tool.


When written to a file, the output only differs by the additional `saturation,oil` and `relperm,oil` columns as in the following extract,

.. code:: sh

  # column 1 = index
  # column 2 = saturation,gas
  # column 3 = saturation,oil
  # column 4 = saturation,water
  # column 5 = relperm,gas
  # column 6 = relperm,oil
  # column 7 = relperm,water
    0.0000e+00 1.0000e+00 0.0000e+00 0.0000e+00 7.5000e-01 0.0000e+00 0.0000e+00
    1.0000e-03 9.9913e-01 8.6667e-04 0.0000e+00 7.5000e-01 0.0000e+00 0.0000e+00
    2.0000e-03 9.9827e-01 1.7333e-03 0.0000e+00 7.5000e-01 0.0000e+00 0.0000e+00
    3.0000e-03 9.9740e-01 2.6000e-03 0.0000e+00 7.5000e-01 0.0000e+00 0.0000e+00
    4.0000e-03 9.9653e-01 3.4667e-03 0.0000e+00 7.5000e-01 0.0000e+00 0.0000e+00
    5.0000e-03 9.9567e-01 4.3333e-03 0.0000e+00 7.5000e-01 0.0000e+00 0.0000e+00
    6.0000e-03 9.9480e-01 5.2000e-03 0.0000e+00 7.5000e-01 0.0000e+00 0.0000e+00
    7.0000e-03 9.9393e-01 6.0667e-03 0.0000e+00 7.5000e-01 0.0000e+00 0.0000e+00
    8.0000e-03 9.9307e-01 6.9333e-03 0.0000e+00 7.5000e-01 0.0000e+00 0.0000e+00



Then a ternary diagram such as the following for ``StoneII`` model. 

.. plot:: ../inputFiles/constitutiveDriver/testRelperm_data/plot_relperm_ternary_table.py


While ``Baker`` constructs a weighted sum of linear interpolant :math:`kr_{ow}` and :math:`kr_{og}` (see :ref:`ThreePhaseRelativePermeability`) is an easier and more stable first approach, ``StoneII`` is more physically consistent in carbonates and sandstones. It exhibits the characteristic straight isoperm patterns (as shown above).
