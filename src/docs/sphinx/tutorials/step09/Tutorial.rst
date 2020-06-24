.. _TutorialHydraulicFractureWithAdvancedXML:

#####################################################
Tutorial 10: Hydraulic Fracturing 
#####################################################

**Context**

In this tutorial, we use a fully coupled hydrofracture solver from GEOSX to solve for the propagation of a single within a reservoir with hetrogeneous in-situ properties.
Advanced xml features will be used throughout the example.

**Objectives**

At the end of this tutorial you will know:

  - how to use multiple solvers for hydraulic fracturing problems,
  - how to specify pre-existing fractures and where new fractures can develop,
  - how to construct a mesh with bias
  - how to specify heterogeneous in-situ properties and initial conditions
  - how to use parameters, symbolic math, and units in xml files.


**Input files**

This tutorial uses a set of input files and table files located at:

.. code-block:: console

  examples/hydraulicFracturing/heterogeneousInSituProperties

Note: because these files use the advanced xml features, they must be preprocessed using pygeos.
To install pygeos, see :ref:`advanced_xml_features`

------------------------------------------------------------------
Description of the case
------------------------------------------------------------------

Here, our goal is to demonstrate how hydraulic fractures are modeled in an typical environment.
The in-situ properties and initial conditions are a based upon a randomly generated, fractal, 1D layer-cake model.


.. image:: hf_example.png



------------------------------------------------------------------
Preparing the input files
------------------------------------------------------------------

The inputs for this case are contained inside a case-specific (``heterogeneousInSitu_singleFracture.xml``) and base (``heterogeneousInSitu_base.xml``) XML file.
The ``tables`` directory contains the pre-constructed geologic model.
This tutorial will first focus on the case-specific input file, which contains the key parameter definitions, then consider the base-file.





Included: including external xml files
---------------------------------------------

At the head of the case-specific xml file is a block that will instruct GEOSX to include an external file.
In our case, this points to the base hydraulic fracturing input file.

.. literalinclude:: ../../../../../examples/hydraulicFracturing/heterogeneousInSituProperties/heterogeneousInSitu_singleFracture.xml
  :language: xml
  :start-after: <!-- SPHINX_HYDROFRACTURE_INCLUDED -->
  :end-before: <!-- SPHINX_HYDROFRACTURE_INCLUDED_END -->




Parameters: defining variables to be used throughout the file
--------------------------------------------------------------

The ``Parameters`` block defines a series of variables that can be used throughout the input file.
These variables allow a given input file to be easily understood and/or modified for a specific environment, even by non-expert users.

Explain types of parameters (string, real64, etc.), units, symbolic math, how parameters are used later in the input.

.. literalinclude:: ../../../../../examples/hydraulicFracturing/heterogeneousInSituProperties/heterogeneousInSitu_singleFracture.xml
  :language: xml
  :start-after: <!-- SPHINX_HYDROFRACTURE_PARAMETERS -->
  :end-before: <!-- SPHINX_HYDROFRACTURE_PARAMETERS_END -->


Mesh: building a mesh with biased boundaries
--------------------------------------------------------------

.. literalinclude:: ../../../../../examples/hydraulicFracturing/heterogeneousInSituProperties/heterogeneousInSitu_singleFracture.xml
  :language: xml
  :start-after: <!-- SPHINX_HYDROFRACTURE_MESH -->
  :end-before: <!-- SPHINX_HYDROFRACTURE_MESH_END -->


Geometry: defining a fracture nodeset
--------------------------------------------------------------

.. literalinclude:: ../../../../../examples/hydraulicFracturing/heterogeneousInSituProperties/heterogeneousInSitu_singleFracture.xml
  :language: xml
  :start-after: <!-- SPHINX_HYDROFRACTURE_GEOMETRY -->
  :end-before: <!-- SPHINX_HYDROFRACTURE_GEOMETRY_END -->


Boundary Conditions: defining boundary conditions
-----------------------------------------------------------------

.. literalinclude:: ../../../../../examples/hydraulicFracturing/heterogeneousInSituProperties/heterogeneousInSitu_singleFracture.xml
  :language: xml
  :start-after: <!-- SPHINX_HYDROFRACTURE_BC -->
  :end-before: <!-- SPHINX_HYDROFRACTURE_BC_END -->


.. literalinclude:: ../../../../../examples/hydraulicFracturing/heterogeneousInSituProperties/heterogeneousInSitu_base.xml
  :language: xml
  :start-after: <!-- SPHINX_HYDROFRACTURE_BC_BASE -->
  :end-before: <!-- SPHINX_HYDROFRACTURE_BC_BASE_END -->


Solvers: setting up the coupled hydraulic fracturing solver
-----------------------------------------------------------------


.. literalinclude:: ../../../../../examples/hydraulicFracturing/heterogeneousInSituProperties/heterogeneousInSitu_base.xml
  :language: xml
  :start-after: <!-- SPHINX_HYDROFRACTURE_SOLVERS -->
  :end-before: <!-- SPHINX_HYDROFRACTURE_SOLVERS_END -->


Events: setting up flexible events
-----------------------------------------------------------------


.. literalinclude:: ../../../../../examples/hydraulicFracturing/heterogeneousInSituProperties/heterogeneousInSitu_base.xml
  :language: xml
  :start-after: <!-- SPHINX_HYDROFRACTURE_EVENTS -->
  :end-before: <!-- SPHINX_HYDROFRACTURE_EVENTS_END -->


Functions: building functions to set in-situ properties
-----------------------------------------------------------------


.. literalinclude:: ../../../../../examples/hydraulicFracturing/heterogeneousInSituProperties/heterogeneousInSitu_base.xml
  :language: xml
  :start-after: <!-- SPHINX_HYDROFRACTURE_FUNCTIONS -->
  :end-before: <!-- SPHINX_HYDROFRACTURE_FUNCTIONS_END -->




------------------------------------------------------------------
Running the case and inspecting the results
------------------------------------------------------------------

Preprocessing the input file
---------------------------------

Because we are using advanced xml features in this example, the input file must be pre-processed using pygeos.
To build the final input file ``hydrofracture_processed.xml``, run the following:

``geosx_bin_dir/pygeos examples/hydraulicFracturing/heterogeneousInSituProperties/heterogeneousInSitu_singleFracture.xml -o hydrofracture_processed.xml``


Running the case
---------------------------------

To run the case, use the following command:

``geosx_bin_dir/geosx -i hydrofracture_processed.xml``

When it is finished, if successful, you should see something like this:

.. code-block:: sh

  Cleaning up events
  Rank 0: Writing out restart file at poroElastic_Terzaghi_restart_000000014/rank_0000000.hdf5

  init time = 0.015293s, run time = 0.44605s
  Umpire            HOST high water mark:  540.6 KB



Inspecting results
---------------------------------

We have requested VisIt formatted output in this file.
We can therefore import our Silo files into VisIt and visualize the outcome.





------------------------------------------------------------------
To go further
------------------------------------------------------------------


**Feedback on this tutorial**

This concludes the poroelastic tutorial.
For any feedback on this tutorial, please submit a `GitHub issue on the project's GitHub page <https://github.com/GEOSX/GEOSX/issues>`_.



**For more details**

  - More on advanced xml features, please see :ref:`advanced_xml_features`.
  - More on functions, please see :ref:`FunctionManager`.


