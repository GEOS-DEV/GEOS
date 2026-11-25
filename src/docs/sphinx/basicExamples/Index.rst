.. _BasicExamples:


Basic Examples
=================================

These examples show complete problems organized by physical process.
They assume that you have worked through the :ref:`Tutorials`, and focus on the features specific to each process.
Each example runs on a laptop or a small number of MPI processes.

A suggested order is:

 #. :ref:`TutorialDeadOilBottomLayersSPE10`: two-phase dead-oil flow in a heterogeneous reservoir (two layers of SPE10).
 #. :ref:`TutorialDeadOilEgg`: adds wells, coupled to the reservoir, with rate and pressure constraints (Egg model).
 #. :ref:`TutorialCO2FieldCase`: CO\ :sub:`2` injection with a CO\ :sub:`2`-brine fluid model, running in parallel, and restarting.
 #. :ref:`TutorialPoroelasticity`: coupled flow and geomechanics (Terzaghi's consolidation problem).
 #. :ref:`TutorialHydraulicFractureWithAdvancedXML`: hydraulic fracture propagation, using XML parameters and preprocessing.
 #. :ref:`triaxialDriverExample`: testing a constitutive model at a single material point, without a mesh.

More examples, including validation studies against analytical solutions and benchmarks, are found in the :ref:`AdvancedExamples`.

.. toctree::
   :maxdepth: 1

   multiphaseFlow/Example

   multiphaseFlowWithWells/Example

   co2Injection/Example

   poromechanics/Example

   hydraulicFracturing/Example

   triaxialDriver/Example
