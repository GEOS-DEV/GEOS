.. _TutorialPoroelasticity:


##############################################
Tutorial 8: Terzaghi's poroelastic problem
##############################################


**Context**

In this tutorial, we use a coupled solver to solve
a poroelastic Terzaghi-type problem, a classic benchmark in poroelasticity.
We do so by coupling a single phase flow solver
with a small-strain Lagrangian mechanics solver.


**Objectives**

At the end of this tutorial you will know:

  - how to use multiple solvers for poroelastic problems,
  - how to define finite elements and finite volume numerical methods.


**Input file**

This tutorial uses no external input files and everything required is
contained within a single GEOSX input file.
The xml input file for this test case is located at:

.. code-block:: console

  PoroElastic_Terzaghi_FIM.xml




------------------------------------------------------------------
Description of the case
------------------------------------------------------------------

We simulate the consolidation of a poroelastic
fluid-saturated column of height :math:`L` having unit cross-section.
The column is instantaneously loaded at time :math:`t` = 0 s with a constant compressive traction :math:`w` applied on the face highlighted in red in the figure below.

.. _problemSketchFig:
.. figure:: terzaghi_sketch.svg
   :align: center
   :width: 400
   :figclass: align-center

   Sketch of the the setup for Terzaghi's problem.

GEOSX will calculate the local deformation and pressure changes along the column as a function of time.
We will use the analytical solution for pressure to check the accuracy of the solution obtained with GEOSX, namely

.. math::
   p(x,t) = \frac{4}{\pi} p_0 \sum_{m=0}^{\infty}
            \frac{1}{2m + 1}
            \text{exp} \left[ -\frac{(2m + 1)^2 \pi^2 c_c t}{4 L^2} \right]
            \text{sin} \left[ \frac{(2m+1)\pi x}{2L} \right]

TODO: complete description of symbols (utilized also later in the script that compares numerical and analytical solutions)

------------------------------------------------------------------
Preparing the input files
------------------------------------------------------------------

All inputs for this case are contained inside a single XML file.
In this tutorial, we focus our attention on the ``Solvers`` tags,
the ``NumericalMethods`` tags, and we will briefly inspect the mesh
and field specification tags.


Solvers: setting up a multiphysics coupling
---------------------------------------------

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

 - the single-physics flow solver, a solver of type ``SinglePhaseFVM`` called here ``SinglePhaseFlow`` (more information on these solvers at :ref:`SinglePhaseFlow`),
 - the small-stress Lagrangian mechanics solver, a solver of type ``SolidMechanicsLagrangianSSLE`` called here ``lagsolve`` for Lagrangian Solver (more information here: :ref:`SolidMechanicsLagrangianFEM`),
 - the coupling solver that will bind the two single-physics solvers above, an object of type ``Poroelastic`` called here ``poroSolve`` (more information at :ref:`PoroelasticSolver`).

Note that the ``name`` attribute of these solvers is
chosen by the user and is not imposed by GEOSX.

The two single-physics solvers are parameterized as explained
in their respective documentation, each with their own tolerances,
verbosity levels, target regions,
and other solver-specific attributes.

Let us focus our attention on the coupling solver.
This solver (``poroSolve``) uses a set of attributes that specifically describe the
coupling for a poroelastic framework.
For instance, we must point this solver to the correct fluid solver (here: ``SinglePhaseFlow``), the correct solid solver (here: ``lagsolve``).
Now that these two solvers are tied together inside the coupling solver,
we have a coupled multiphysics problem defined.
More parameters are required to characterize a coupling.
Here, we specify the coupling type (``FIM``, fully implicit method; a choice among several possible options),
the discretization method (``FE1``, defined further in the input file),
and the target regions (here, we only have one, ``Region2``).


.. literalinclude:: PoroElastic_Terzaghi_FIM.xml
  :language: xml
  :start-after: <!-- SPHINX_POROELASTIC_SOLVER -->
  :end-before: <!-- SPHINX_POROELASTIC_SOLVER_END -->


Multiphysics numerical methods
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

.. literalinclude:: PoroElastic_Terzaghi_FIM.xml
  :language: xml
  :start-after: <!-- SPHINX_POROELASTIC_NUMERICAL_METHODS -->
  :end-before: <!-- SPHINX_POROELASTIC_NUMERICAL_METHODS_END -->


Setting up the mesh
---------------------------------

Last, let us take a closer look at the geometry of this simple problem.
We use the internal mesh generator to create a beam-like mesh,
with one single element along the Y and Z axes, and 21 elements along the X axis.
All the elements are hexahedral elements (C3D8) of the same dimension (2x1x1 meters).

We also define a pair of geometric boxes that will help us
locate and specify our boundary conditions. These boundary conditions are defined under the ``FieldSpecifications`` tag.

.. literalinclude:: PoroElastic_Terzaghi_FIM.xml
  :language: xml
  :start-after: <!-- SPHINX_POROELASTIC_MESH -->
  :end-before: <!-- SPHINX_POROELASTIC_MESH_END -->



To give some physical meaningfulness to our problem,
we specify the material properties of our beam domain in the
``Constitutive`` and ``FieldSpecifications`` sections.
Here, we have a homogeneous shaly material with a 0.3 porosity and an isotropic permeability of about 50mD, completely saturated with water.
The shale has a bulk modulus of 61.9 MPa, and a density of 2,700 kg/m3. Our beam is subject to compression forces exerted at one if its extremities.

As shown in the ``Events`` section, we run this simulation for 2,000 seconds. We use a function specification to change our constraints in time. For more on functions, see :ref:`FunctionManager`.



------------------------------------------------------------------
Running the case and inspecting the results
------------------------------------------------------------------

Running the case
---------------------------------

To run the case, use the following command:

``path/to/geosx -i src/coreComponents/physicsSolvers/multiphysics/integratedTests/PoroElastic_Terzaghi_FIM.xml``

When it is finished, if successful, you should see something like this:

.. code-block:: sh

  Cleaning up events
  Rank 0: Writing out restart file at PoroElastic_Terzaghi_FIM_restart_000000014/rank_0000000.hdf5

  init time = 0.015293s, run time = 0.44605s
  Umpire            HOST high water mark:  540.6 KB


Inspecting the console output
---------------------------------

Depending on the individual level of log verbosity,
coupled solvers may display information about the Newton
convergence of each single-physics solver.


Here, we see for instance the ``RSolid`` and ``RFluid`` at a representative timestep
(residual values for solid and fluid mechanics solvers, respectively)

.. code-block:: sh

   Attempt:  0, NewtonIter:  0
   ( RSolid ) = (2.54e-16) ;     ( Rsolid, Rfluid ) = ( 2.54e-16, 8.44e-06 )
   ( R ) = ( 8.44e-06 ) ;
   Attempt:  0, NewtonIter:  1
   ( RSolid ) = (2.22e-16) ;     ( Rsolid, Rfluid ) = ( 2.22e-16, 2.76e-20 )
   ( R ) = ( 2.22e-16 ) ;
   Last LinSolve(iter,res) = (   1, 1.15e-11 ) ;

As expected, since we are dealing with a linear problem,
the fully implicit solver coverges in a single iteration.



Inspecting results
---------------------------------

TODO: Description to be added.

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


   class terzaghi:

       def __init__(self, hydromechanicalParameters, xMin, xMax, appliedTraction):
           E = hydromechanicalParameters["youngModulus"]
           nu = hydromechanicalParameters["poissonRation"]
           b = hydromechanicalParameters["biotCoefficient"]
           mu = hydromechanicalParameters["fluidViscosity"]
           cf = hydromechanicalParameters["fluidCompressibility"]
           phi = hydromechanicalParameters["porosity"]
           k = hydromechanicalParameters["permeability"]

           K = E / 3.0 / (1.0 - 2.0 * nu) # bulk modulus
           Kv = E * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu )) # uniaxial bulk modulus
           Se = (b - phi) * (1.0 - b) / K + phi * cf          # constrained specific storage

           self.characteristicLength = xMax - xMin
           self.appliedTraction = abs(appliedTraction)
           self.loadingEfficiency = b / (Kv * Se + b**2)
           self.consolidationCoefficient = (k / mu) * Kv / (Se * Kv + b**2)
           self.consolidationTime = self.characteristicLength**2 / self.consolidationCoefficient
           self.initialPressure = self.loadingEfficiency * self.appliedTraction

       def computePressure(self, x, t):
           if  t == 0.0:
               return self.initialPressure
           else:
               cc = self.consolidationCoefficient
               L = self.characteristicLength
               p = nsum(lambda m: 1 / (2 * m + 1)
                 * exp(-((2 * m + 1) ** 2) * (math.pi ** 2) * cc * t / 4 / L / L)
                 * sin((2 * m + 1) * math.pi * x / 2 / L), [0, inf])
               return 4 * self.initialPressure / math.pi * p


   def getHydromechanicalParametersFromXML( xmlFilePath ):
       tree = ElementTree.parse(xmlFilePath)

       param1 = tree.find('Constitutive/PoroLinearElasticIsotropic')
       param2 = tree.find('Constitutive/CompressibleSinglePhaseFluid')
       param3 = tree.findall('FieldSpecifications/FieldSpecification')

       found_porosity = False
       found_permeability = False
       porosity = 0.0
       for elem in param3:
           if elem.get("fieldName") == "permeability" and elem.get("component") == "0":
               permeability = float(elem.get("initialCondition")) * float(elem.get("scale"))
               found_permeability = True

           if elem.get("fieldName") == "referencePorosity":
               porosity = float(elem.get("initialCondition")) * float(elem.get("scale"))
               found_porosity = True

           if found_permeability and found_porosity: break

       hydromechanicalParameters = dict.fromkeys(["youngModulus",
                                                  "poissonRation",
                                                  "biotCoefficient",
                                                  "fluidViscosity",
                                                  "fluidCompressibility",
                                                  "porosity",
                                                  "permeability"])

       hydromechanicalParameters["youngModulus"] = float(param1.get("defaultYoungsModulus"))
       hydromechanicalParameters["poissonRation"] = float(param1.get("defaultPoissonRatio"))
       hydromechanicalParameters["biotCoefficient"] = float(param1.get("BiotCoefficient"))
       hydromechanicalParameters["fluidViscosity"] = float(param2.get("defaultViscosity"))
       hydromechanicalParameters["fluidCompressibility"] = float(param2.get("compressibility"))
       hydromechanicalParameters["porosity"] = porosity
       hydromechanicalParameters["permeability"] = permeability

       return hydromechanicalParameters


   def getAppliedTractionFromXML( xmlFilePath ):
       tree = ElementTree.parse(xmlFilePath)

       param = tree.findall('FieldSpecifications/FieldSpecification')

       found_traction = False
       for elem in param:
           if elem.get("fieldName") == "Traction" and elem.get("component") == "0":
               traction = float(elem.get("scale"))
               found_traction = True

           if found_traction: break

       return traction


   def getDomainMaxMinXCoordFromXML(xmlFilePath):
       tree = ElementTree.parse(xmlFilePath)
       meshElement = tree.find('Mesh/InternalMesh')
       nodeXCoords = meshElement.get("xCoords")
       nodeXCoords = [float(i) for i in nodeXCoords[1:-1].split(",")]
       xMin = nodeXCoords[0]
       xMax = nodeXCoords[-1]
       return xMin, xMax


   def reconstructCellCentroidXCoordsFromXML(xmlFilePath):
       # Reconstruct cell centroid list
       # (temporary solution -- centroids should be available in the hdf5)
       tree = ElementTree.parse(xmlFilePath)
       meshElement = tree.find('Mesh/InternalMesh')
       nodeXCoords = meshElement.get("xCoords")
       nCellX = meshElement.get("nx")
       nCellX = [int(i) for i in nCellX[1:-1].split(",")]
       nodeXCoords = [float(i) for i in nodeXCoords[1:-1].split(",")]

       cellCentroidXCoords = np.empty(sum(nCellX)) # cell centroid in x-direction
       istr = 0
       for i in range(len(nCellX)):
           iend = istr + nCellX[i]
           tmp = np.linspace(nodeXCoords[i], nodeXCoords[i+1], nCellX[i] + 1)
           cellCentroidXCoords[istr:iend] = (tmp[1:] + tmp[:-1]) / 2.
           istr = iend

       return cellCentroidXCoords


   def main():
       # File path
       hdf5FilePath = "pressure_history.hdf5"
       xmlFilePath = "PoroElastic_Terzaghi_FIM.xml"

       # Read HDF5
       hf = h5py.File(hdf5FilePath, 'r')
       time = hf.get('Time')
       time = np.array(time)
       pressure = hf.get('pressure')
       pressure = np.array(pressure)

       # Extract info from XML
       hydromechanicalParameters = getHydromechanicalParametersFromXML(xmlFilePath)
       appliedTraction = getAppliedTractionFromXML(xmlFilePath)

       # Get domain min/max coordinate in the x-direction
       xMin, xMax =getDomainMaxMinXCoordFromXML(xmlFilePath)

       # (temporary solution -- centroids will be read in the hdf5)
       x = reconstructCellCentroidXCoordsFromXML(xmlFilePath)

       # Initialize Terzaghi's analytical solution
       terzaghiAnalyticalSolution = terzaghi(hydromechanicalParameters, xMin, xMax, appliedTraction)

       # Plot analytical (continuous line) and numerical (markers) pressure solution
       x_analytical = np.linspace(xMin, xMax, 51, endpoint=True)
       pressure_analytical = np.empty(len(x_analytical))

       cmap = plt.get_cmap("tab10")
       iplt = -1
       for k in range(0, len(time), 2):
           iplt += 1
           t = time[k,0]
           i = 0
           for xCell in x_analytical:
               xScaled = terzaghiAnalyticalSolution.characteristicLength * (xCell - xMin) / (xMax - xMin)
               pressure_analytical[i] = terzaghiAnalyticalSolution.computePressure(xScaled, t)
               i += 1
           plt.plot(x_analytical, pressure_analytical, color=cmap(iplt), label='t = ' + str(t) + ' s')
           plt.plot(x, pressure[k, :], 'o', color=cmap(iplt))

       plt.grid()
       plt.xlabel('time [s]')
       plt.ylabel('pressure [Pa]')
       plt.legend(bbox_to_anchor=(0.05, 0.5), loc='lower left', borderaxespad=0.)
       plt.show()

   if __name__ == "__main__":
       main()






------------------------------------------------------------------
To go further
------------------------------------------------------------------


**Feedback on this tutorial**

This concludes the poroelastic tutorial.
For any feedback on this tutorial, please submit a `GitHub issue on the project's GitHub page <https://github.com/GEOSX/GEOSX/issues>`_.



**For more details**

  - More on poroelastic multiphysics solvers, please see :ref:`PoroelasticSolver`.
  - More on numerical methods, please see :ref:`NumericalMethodsManager`.
  - More on functions, please see :ref:`FunctionManager`.
