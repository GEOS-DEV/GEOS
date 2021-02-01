.. _TutorialSneddon:


#######################################################
Tutorial X: Pressurized fracture in an infinite medium
#######################################################


**Context**

In this tutorial, we use a coupled solver to solve a poroelastic Terzaghi-type
problem, a classic benchmark in poroelasticity.
We do so by coupling a single phase flow solver with a small-strain Lagrangian mechanics solver.


**Objectives**

At the end of this tutorial you will know:

  - how to define embedded fractures in the porous domain,
  - how to use the SolidMechanicsEmbeddedFractures solver to solve mechanics problems with
  embedded fractures.


**Input file**

This tutorial uses no external input files and everything required is
contained within a single GEOSX input file.
The xml input file for this test case is located at:

.. code-block:: console

  src/coreComponents/physicsSolvers/solidMechanics/benchmarks/Sneddon-Validation.xml


------------------------------------------------------------------
Description of the case
------------------------------------------------------------------

We compute the displacement field induced by the presence of a pressurized fracture,
of length :math:`L_f`, in a porous medium.

.. _problemSketchFig:
.. figure:: sneddon_sketch.svg
   :align: center
   :width: 500
   :figclass: align-center

   Sketch of the the setup for Terzaghi's problem.

GEOSX will calculate the displacement field in the porous matrix and the displacement
jump at the fracture surface.
We will use the analytical solution for the fracture aperture, :math:`w_n` (normal component of the
jump) to, i.e.

.. math::
   w_n (s) = \frac{4(1 - \nu^2)p_f}{E} \, \sqrt{ \frac{L_f^2}{4} - s^2 }

where
- :math:`E` is the Young's modulus
- :math:`\nu` is the Poisson's ratio
- :math:`p_f` is the fracture pressure
- :math:`s` is the local fracture coordinate in :math:`[-\frac{L_f}{2}, \frac{L_f}{2}]`

------------------------------------------------------------------
Preparing the input files
------------------------------------------------------------------

All inputs for this case are contained inside a single XML file.
In this tutorial, we focus our attention on the ``Solvers`` tags,
the ``ElementRegions`` tags.

Solvers: setting up the embedded fractures mechanics solver
-----------------------------------------------------------

GEOSX is a multi-physics tool. Different combinations of
physics solvers available in the code can be applied
in different regions of the mesh at different moments of the simulation.
The XML ``Solvers`` tag is used to list and parameterize these solvers.

To specify a coupling between two solvers, as done here,
we define and characterize each single-physics solver separately.
Then, we define a *coupling solver* between these single-physics
solvers as another, separate, solver.
This approach allows for generality and flexibility in our multi-physics resolutions.
The order in which these solver specifications is done is not important.
It is important, though, to instantiate each single-physic solvers
with meaningful names. The names given to these single-physics solver instances
will be used to recognize them and create the coupling.

To define a poroelastic coupling, we will effectively define three solvers:

 - the single-physics flow solver, a solver of type ``SinglePhaseFVM`` called here ``SinglePhaseFlowSolver`` (more information on these solvers at :ref:`SinglePhaseFlow`),
 - the small-stress Lagrangian mechanics solver, a solver of type ``SolidMechanicsLagrangianSSLE`` called here ``LinearElasticitySolver`` (more information here: :ref:`SolidMechanicsLagrangianFEM`),
 - the coupling solver that will bind the two single-physics solvers above, an object of type ``Poroelastic`` called here ``PoroelasticitySolver`` (more information at :ref:`PoroelasticSolver`).

Note that the ``name`` attribute of these solvers is
chosen by the user and is not imposed by GEOSX.

The two single-physics solvers are parameterized as explained
in their respective documentation, each with their own tolerances,
verbosity levels, target regions,
and other solver-specific attributes.

Let us focus our attention on the coupling solver.
This solver (``PoroelasticitySolver``) uses a set of attributes that specifically describe the coupling for a poroelastic framework.
For instance, we must point this solver to the correct fluid solver (here: ``SinglePhaseFlowSolver``), the correct solid solver (here: ``LinearElasticitySolver``).
Now that these two solvers are tied together inside the coupling solver,
we have a coupled multiphysics problem defined.
More parameters are required to characterize a coupling.
Here, we specify the coupling type (``FIM``, fully implicit method; a choice among several possible options),
the discretization method (``FE1``, defined further in the input file),
and the target regions (here, we only have one, ``Region1``).

.. literalinclude:: ../../../../coreComponents/physicsSolvers/solidMechanics/benchmarks/Sneddon-Validation.xml
  :language: xml
  :start-after: <!-- SPHINX_SNEDDON_SOLVER -->
  :end-before: <!-- SPHINX_SNEDDON_SOLVER_END -->

  Setting up mesh, material properties and boundary conditions
  --------------------------------------------------------------------

  Last, let us take a closer look at the geometry of this simple problem.
  We use the internal mesh generator to create a beam-like mesh,
  with one single element along the Y and Z axes, and 21 elements along the X axis.
  All the elements are hexahedral elements (C3D8) of the same dimension (1x1x1 meters).

  .. literalinclude:: ../../../../coreComponents/physicsSolvers/solidMechanics/benchmarks/Sneddon-Validation.xml
    :language: xml
    :start-after: <!-- SPHINX_POROELASTIC_MESH -->
    :end-before: <!-- SPHINX_POROELASTIC_MESH_END -->

  The parameters used in the simulation are summarized in the following table.

  +----------------+-----------------------+------------------+-------------------+
  | Symbol         | Parameter             | Units            | Value             |
  +================+=======================+==================+===================+
  | :math:`E`      | Young's modulus       | [Pa]             | 1.0*10\ :sup:`4`  |
  +----------------+-----------------------+------------------+-------------------+
  | :math:`\nu`    | Poisson's ration      | [-]              | 0.2               |
  +----------------+-----------------------+------------------+-------------------+
  | :math:`L_f`    | Fracture length       | [m]              |                   |
  +----------------+-----------------------+------------------+-------------------+
  | :math:`p_f`    | fracture pressure     | [Pa]             | 10.0              |
  +----------------+-----------------------+------------------+-------------------+

  Material properties and boundary conditions are specified in the
  ``Constitutive`` and ``FieldSpecifications`` sections.

Adding an embedded fracture
---------------------------------

Numerical methods in multiphysics settings are similar to single physics numerical methods. All can be defined under the same ``NumericalMethods`` XML tag.
In this problem, we use finite volume for flow and finite elements for solid mechanics.
Both of these methods require additional parameterization attributes to be defined here.

As we have seen before, the coupling solver and the solid mechanics solver require the specification of a discretization method called ``FE1``.
This discretization method is defined here as a finite element method
using linear basis functions and Gaussian quadrature rules.
For more information on defining finite elements numerical schemes,
please see the dedicated :ref:`FiniteElementDiscretization` section.

The finite volume method requires the specification of a discretization scheme.
Here, we use a two-point flux approximation as described in the dedicated documentation (found here: :ref:`FiniteVolumeDiscretization`).

.. literalinclude:: ../../../../coreComponents/physicsSolvers/solidMechanics/benchmarks/Sneddon-Validation.xml
  :language: xml
  :start-after: <!-- SPHINX_SNEDDON_GEOMETRY -->
  :end-before: <!-- SPHINX_SNEDDON_GEOMETRY_END -->

------------------------------------------------------------------
Running the case and inspecting the results
------------------------------------------------------------------

Running the case
---------------------------------

To run the case, use the following command:

``path/to/geosx -i src/coreComponents/physicsSolvers/solidMechanics/benchmarks/Sneddon-Validation.xml``

Inspecting results
---------------------------------

This plot compares the analytical pressure solution (continuous lines) at selected
times with the numerical solution (markers).

.. plot::

   import matplotlib
   import matplotlib.pyplot as plt
   import numpy as np
   import h5py
   import xml.etree.ElementTree as ElementTree
   from mpmath import *
   import math


   class Sneddon:

       def __init__(self, mechanicalParameters, length, pressure):
           K = mechanicalParameters["bulkModulus"]
           G = mechanicalParameters["shearModulus"]
           E = (9 * K * G) / (3*K+G)
           nu = E / (2 * G) - 1

           self.scaling = ( 4 * (1 - nu**2) ) * pressure / E;
           self.halfLength = length;

       def computeAperture(self, x):
           return self.scaling * ( self.halfLength**2  - x**2 )**0.5;


       def getMechanicalParametersFromXML( xmlFilePath ):

           tree = ElementTree.parse(xmlFilePath)
           param = tree.find('Constitutive/LinearElasticIsotropic')

           mechanicalParameters = dict.fromkeys(["bulkModulus", "shearModulus"])
           mechanicalParameters["bulkModulus"] = float(param.get("defaultBulkModulus"))
           mechanicalParameters["shearModulus"] = float(param.get("defaultShearModulus"))
           return mechanicalParameters


  def getFracturePressureFromXML( xmlFilePath ):
      tree = ElementTree.parse(xmlFilePath)

      param = tree.findall('FieldSpecifications/FieldSpecification')

      found_traction = False
      for elem in param:
          if elem.get("fieldName") == "fractureTraction" and elem.get("component") == "0":
              pressure = float(elem.get("scale"))*(-1)
              found_traction = True
          if found_traction: break

      return pressure


   def getFractureLengthFromXML(xmlFilePath):
       tree = ElementTree.parse(xmlFilePath)
       boundedPlane = tree.find('Geometry/BoundedPlane')
       dimensions = boundedPlane.get("dimensions")
       dimensions = [float(i) for i in dimensions[1:-1].split(",")]
       length = dimensions[0] / 2
       origin = boundedPlane.get("origin")
       origin = [float(i) for i in origin.split(",")]

       return length, origin[0]


   def main():
       # File path
       hdf5File1Path = "displacemenJump_history.hdf5"
       hdf5File2Path = "cell_centers.hdf5"
       xmlFilePath = "../../../../coreComponents/physicsSolvers/solidMechanics/benchmarks/Sneddon-Validation.xml"

       # Read HDF5
       hf = h5py.File(hdf5File1Path, 'r')
       jump = hf.get('displacementJump')
       jump = np.array(jump)
       aperture = jump[0,:,0]

       hf = h5py.File(hdf5File2Path, 'r')
       x = hf.get('elementCenter')
       x = x[0,:,0]

       # Filter out extra entries in the hdf5 file. It is just to make the plot look nicer
       voidIndexes = np.asarray( np.where(x == 0) )
       if voidIndexes.size !=0:
         lastValue = voidIndexes[0][0]
         aperture = aperture[0:lastValue]
         x = x[0:lastValue]

      # Extract info from XML
      mechanicalParameters = getMechanicalParametersFromXML(xmlFilePath)
      appliedPressure = getFracturePressureFromXML(xmlFilePath)

      # Get length of the fracture
      length, origin = getFractureLengthFromXML(xmlFilePath)

      x = x - origin
      # Initialize Sneddon's analytical solution
      sneddonAnalyticalSolution = Sneddon(mechanicalParameters, length, appliedPressure)

      # Plot analytical (continuous line) and numerical (markers) aperture solution
      x_analytical = np.linspace(-length, length, 101, endpoint=True)
      aperture_analytical = np.empty(len(x_analytical))

      cmap = plt.get_cmap("tab10")
      i=0
      for xCell in x_analytical:
         aperture_analytical[i] = sneddonAnalyticalSolution.computeAperture( xCell )
         i += 1
      plt.plot(x_analytical, aperture_analytical, color=cmap(-1), label='analytical solution')
      plt.plot(x, aperture, 'o', color=cmap(2), label='numerical solution')

      plt.grid()
      plt.xlabel('length [m]')
      plt.ylabel('aperture [m]')
      plt.legend(bbox_to_anchor=(0.5, 0.2), loc='center', borderaxespad=0.)
      plt.show()

   if __name__ == "__main__":
       main()


------------------------------------------------------------------
To go further
------------------------------------------------------------------

**Feedback on this tutorial**

This concludes the Sneddon tutorial.
For any feedback on this tutorial, please submit a `GitHub issue on the project's GitHub page <https://github.com/GEOSX/GEOSX/issues>`_.
