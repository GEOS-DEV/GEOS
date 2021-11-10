.. _TutorialSneddon:


#######################################################
Sneddon's Problem
#######################################################

**Objectives**

At the end of this example you will know:

  - how to define fractures in a porous medium,
  - how to use various solvers (EmbeddedFractures, LagrangianContact and HydroFracture) to solve the mechanics problems with fractures.


**Input file**

This example uses no external input files and everything required is contained within GEOSX input files.

The xml input files for the case with EmbeddedFractures solver are located at:

.. code-block:: console

  inputFiles/efemFractureMechanics/Sneddon_embeddedFrac_base.xml
  inputFiles/efemFractureMechanics/Sneddon_embeddedFrac_benchmark.xml


The xml input files for the case with LagrangianContact solver are located at:

.. code-block:: console

  inputFiles/lagrangianContactMechanics/Sneddon_contactMechanics_base.xml
  inputFiles/lagrangianContactMechanics/Sneddon_contactMechanics_benchmark.xml


The xml input files for the case with HydroFracture solver are located at:

.. code-block:: console

  inputFiles/hydraulicFracturing/Sneddon_hydroFrac_base.xml
  inputFiles/hydraulicFracturing/Sneddon_hydroFrac_benchmark.xml


------------------------------------------------------------------
Description of the case
------------------------------------------------------------------

We compute the displacement field induced by the presence of a pressurized fracture,
of length :math:`L_f`, in a porous medium.

GEOSX will calculate the displacement field in the porous matrix and the displacement
jump at the fracture surface.
We will use the analytical solution for the fracture aperture, :math:`w_n` (normal component of the
jump), to verify the numerical results

.. math::
   w_n (s) = \frac{4(1 - \nu^2)p_f}{E} \, \sqrt{ \frac{L_f^2}{4} - s^2 }

where
- :math:`E` is the Young's modulus
- :math:`\nu` is the Poisson's ratio
- :math:`p_f` is the fracture pressure
- :math:`s` is the local fracture coordinate in :math:`[-\frac{L_f}{2}, \frac{L_f}{2}]`


In this example, we focus our attention on the ``Solvers``, the ``ElementRegions``, and the ``Geometry`` tags.

-----------------------------------------------------------
Mechanics solver
-----------------------------------------------------------

To define a mechanics solver capable of including embedded fractures, we will
define two solvers:

 - a ``SolidMechanicsEmbeddedFractures`` solver, called ``mechSolve``
 - a small-strain Lagrangian mechanics solver, of type ``SolidMechanicsLagrangianSSLE`` called here ``matrixSolver`` (see: :ref:`SolidMechanicsLagrangianFEM`)

Note that the ``name`` attribute of these solvers is chosen by the user and is not imposed by GEOSX. 
It is important to make sure that the ``solidSolverName`` specified in the embedded fractures solver corresponds to the
small-strain Lagrangian solver used in the matrix.

The two single-physics solvers are parameterized as explained in their respective documentation, each with their own tolerances,
verbosity levels, target regions, and other solver-specific attributes.

Additionally, we need to specify another solver of type, ``EmbeddedSurfaceGenerator``, which is used to discretize the fracture planes.

.. literalinclude:: ../../../../../../inputFiles/efemFractureMechanics/Sneddon_embeddedFrac_base.xml
  :language: xml
  :start-after: <!-- SPHINX_SNEDDON_SOLVER -->
  :end-before: <!-- SPHINX_SNEDDON_SOLVER_END -->


To setup a coupling between rock and fracture deformations in LagrangianContact solver, we define three different solvers:

- For solving the frictional contact, we define a Lagrangian contact solver, called here ``lagrangiancontact``. In this solver, we specify ``targetRegions`` that include both the continuum region ``Region`` and the discontinuum region ``Fracture``  where the solver is applied to couple rock and fracture deformations. The contact constitutive law used for the fracture elements is named ``fractureMaterial``,  and is defined later in the ``Constitutive`` section. 

- Rock deformations are handled by a solid mechanics solver ``SolidMechanicsLagrangianSSLE``. The problem runs in ``QuasiStatic`` mode without inertial effects. The computational domain is discretized by ``FE1``, which is defined in the ``NumericalMethods`` section. The solid material is named ``rock`` and its mechanical properties are specified later in the ``Constitutive`` section.

- The solver ``SurfaceGenerator`` defines the fracture region and rock toughness.


.. literalinclude:: ../../../../../../inputFiles/lagrangianContactMechanics/Sneddon_contactMechanics_base.xml
  :language: xml
  :start-after: <!-- SPHINX_SNEDDON_SOLVER -->
  :end-before: <!-- SPHINX_SNEDDON_SOLVER_END -->


Three elementary solvers are combined in the solver ``Hydrofracture`` to model the coupling between fluid flow within the fracture, rock deformation, fracture opening/closure and propagation:

- Rock and fracture deformation are modeled by the solid mechanics solver ``SolidMechanicsLagrangianSSLE``. In this solver, we define ``targetRegions`` that includes both the continuum region and the fracture region. The name of the contact constitutive behavior is also specified in this solver by the ``contactRelationName``, besides the ``solidMaterialNames``.

- The single phase fluid flow inside the fracture is solved by the finite volume method in the solver ``SinglePhaseFVM``.

- The solver ``SurfaceGenerator`` defines the fracture region and rock toughness.


.. literalinclude:: ../../../../../../inputFiles/hydraulicFracturing/Sneddon_hydroFrac_base.xml
  :language: xml
  :start-after: <!-- SPHINX_SNEDDON_SOLVER -->
  :end-before: <!-- SPHINX_SNEDDON_SOLVER_END -->


-------------------------------------------------------------------
Events
-------------------------------------------------------------------
For the case with EmbeddedFractures solver, we add multiple events defining solver applications:

- an event specifying the execution of the ``EmbeddedSurfaceGenerator`` to generate the fracture elements.
- a periodic event specifying the execution of the embedded fractures solver.
- three periodic events specifying the output of simulations results.

.. literalinclude:: ../../../../../../inputFiles/efemFractureMechanics/Sneddon_embeddedFrac_benchmark.xml
  :language: xml
  :start-after: <!-- SPHINX_SNEDDON_EVENTS -->
  :end-before: <!-- SPHINX_SNEDDON_EVENTS_END -->

Similar settings are applied for the other two cases.

--------------------------------------------------------------------
Mesh, material properties, and boundary conditions
--------------------------------------------------------------------

Last, let us take a closer look at the geometry of this simple problem, if using EmbeddedFractures solver.
We use the internal mesh generator to create a large domain
(:math:`40\, m \, \times 40 \,  m \, \times 1 \, m`), with one single element
along the Z axes, 121 elements along the X axis and 921 elements along the Y axis.


.. literalinclude:: ../../../../../../inputFiles/efemFractureMechanics/Sneddon_embeddedFrac_benchmark.xml
   :language: xml
   :start-after: <!-- SPHINX_SNEDDON_MESH -->
   :end-before: <!-- SPHINX_SNEDDON_MESH_END -->


The mesh for the case with LagrangianContact solver was also created using the internal mesh generator, as parametrized in the ``InternalMesh`` XML tag. The mesh discretizes the same compational domain (:math:`40\, m \, \times 40 \,  m \, \times 1 \, m`) with 300 x 300 x 1 eight-node brick elements in the x, y, and z directions respectively. 


.. literalinclude:: ../../../../../../inputFiles/lagrangianContactMechanics/Sneddon_contactMechanics_benchmark.xml
   :language: xml
   :start-after: <!-- SPHINX_SNEDDON_MESH -->
   :end-before: <!-- SPHINX_SNEDDON_MESH_END -->


Similarly, the internal mesh generator was used to discretize the same domain (:math:`40\, m \, \times 40 \,  m \, \times 1 \, m`) and generate the mesh for the case with Hydrofracture solver, which contains 280 x 280 x 1 eight-node brick elements in the x, y, and z directions. 


.. literalinclude:: ../../../../../../inputFiles/hydraulicFracturing/Sneddon_hydroFrac_benchmark.xml
   :language: xml
   :start-after: <!-- SPHINX_SNEDDON_MESH -->
   :end-before: <!-- SPHINX_SNEDDON_MESH_END -->


In all the three cases, eight-node hexahedral elements are defined as ``C3D8`` elementTypes, and their collection forms a mesh
with one group of cell blocks named here ``cb1``. 
Refinement is necessary to conform with the fracture geometry specified in the ``Geometry`` section.


The parameters used in the simulation are summarized in the following table.

  +----------------+-----------------------+------------------+-------------------+
  | Symbol         | Parameter             | Units            | Value             |
  +================+=======================+==================+===================+
  | :math:`K`      | Bulk modulus          | [GPa]            | 16.7              |
  +----------------+-----------------------+------------------+-------------------+
  | :math:`G`      | Shear modulus         | [GPa]            | 10.0              |
  +----------------+-----------------------+------------------+-------------------+
  | :math:`L_f`    | Fracture length       | [m]              | 2.0               |
  +----------------+-----------------------+------------------+-------------------+
  | :math:`p_f`    | Fracture pressure     | [MPa]            | -2.0              |
  +----------------+-----------------------+------------------+-------------------+

Note that the internal fracture pressure has a negative value, due to the negative sign convention for compressive stresses in GEOSX. 

Material properties and boundary conditions are specified in the ``Constitutive`` and ``FieldSpecifications`` sections.

---------------------------------
Adding a fracture
---------------------------------

The static fracture is defined by a nodeset occupying a small region within the computation domain, where the fracture tends to open upon internal pressurization:

- The test case with EmbeddedFractures solver:

.. literalinclude:: ../../../../../../inputFiles/efemFractureMechanics/Sneddon_embeddedFrac_base.xml
  :language: xml
  :start-after: <!-- SPHINX_SNEDDON_GEOMETRY -->
  :end-before: <!-- SPHINX_SNEDDON_GEOMETRY_END -->


- The test case with LagrangianContact solver:

.. literalinclude:: ../../../../../../inputFiles/lagrangianContactMechanics/Sneddon_contactMechanics_base.xml
  :language: xml
  :start-after: <!-- SPHINX_SNEDDON_GEOMETRY -->
  :end-before: <!-- SPHINX_SNEDDON_GEOMETRY_END -->


- The test case with HydroFracture solver:

.. literalinclude:: ../../../../../../inputFiles/hydraulicFracturing/Sneddon_hydroFrac_base.xml
  :language: xml
  :start-after: <!-- SPHINX_SNEDDON_GEOMETRY -->
  :end-before: <!-- SPHINX_SNEDDON_GEOMETRY_END -->

To make these cases identical to the analytical example, fracture propagation is not allowed in this example.

--------------------------------
Time history function
--------------------------------

In the ``Tasks`` section, ``PackCollection`` tasks are defined to collect time history information from fields. 
Either the entire field or specified named sets of indices in the field can be collected. 
In this example, a task is specified to output fracture aperture (normal opening); however, for different solvers, different ``fieldName`` and ``objectPath`` should be called: 

- The test case with EmbeddedFractures solver:

.. literalinclude:: ../../../../../../inputFiles/efemFractureMechanics/Sneddon_embeddedFrac_base.xml
  :language: xml
  :start-after: <!-- SPHINX_SNEDDON_TASK -->
  :end-before: <!-- SPHINX_SNEDDON_TASK_END -->


- The test case with LagrangianContact solver:

.. literalinclude:: ../../../../../../inputFiles/lagrangianContactMechanics/Sneddon_contactMechanics_base.xml
  :language: xml
  :start-after: <!-- SPHINX_SNEDDON_TASK -->
  :end-before: <!-- SPHINX_SNEDDON_TASK_END -->


- The test case with Hydrofracture solver:

.. literalinclude:: ../../../../../../inputFiles/hydraulicFracturing/Sneddon_hydroFrac_base.xml
  :language: xml
  :start-after: <!-- SPHINX_SNEDDON_TASK -->
  :end-before: <!-- SPHINX_SNEDDON_TASK_END -->


These tasks are triggered using the ``Event`` manager with a ``PeriodicEvent`` defined for these recurring tasks. 
GEOSX writes output files named after the string defined in the ``filename`` keyword and formatted as HDF5 files. 
The ``TimeHistory`` file contains the collected time history information from each specified time history collector.
This information includes datasets for the simulation time, element center defined in the local coordinate system, and the time history information. 
A Python script is used to read and plot any specified subset of the time history data for verification and visualization. 


---------------------------------
Running GEOSX
---------------------------------

To run these three cases, use the following commands:

``path/to/geosx -i inputFiles/efemFractureMechanics/Sneddon_embeddedFrac_benchmark.xml``

``path/to/geosx -i inputFiles/lagrangianContactMechanics/Sneddon_contactMechanics_benchmark.xml``

``path/to/geosx -i inputFiles/hydraulicFracturing/Sneddon_hydroFrac_benchmark.xml``

---------------------------------
Inspecting results
---------------------------------

This plot compares the analytical solution (continuous lines) with the numerical solutions (markers) for the normal opening of the pressurized fracture. As shown below, consistently, numerical solutions with different solvers correlate very well with the analytical solution.


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

        param = tree.find('Constitutive/ElasticIsotropic')

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
        origin = [float(i) for i in origin[1:-1].split(",")]

        return length, origin[0]


    def main():
        #-------- EmbeddeFrac File path
        hdf5File1Path = "displacementJump_embeddedFrac.hdf5"
        
        # Read HDF5
        hf = h5py.File(hdf5File1Path, 'r')
        jump = hf.get('displacementJump')
        jump = np.array(jump)
        aperture_EmbeddeFrac = jump[0,:,0]
        x = hf.get('displacementJump elementCenter')
        loc_EmbeddeFrac = x[0,:,1]
        
        #-------- Lagrange Contact File path
        hdf5File2Path = "displacementJump_contactMechanics.hdf5" 
       
        # Read HDF5
        hf = h5py.File(hdf5File2Path, 'r')
        jump = hf.get('localJump')
        jump = np.array(jump) 
        aperture_Contact = jump[0,:,0]
        x = hf.get('localJump elementCenter')
        loc_Contact = x[0,:,1]
         
        #-------- HydroFrac File path
        hdf5File3Path = "displacementJump_hydroFrac.hdf5"
        
        # Read HDF5
        hf = h5py.File(hdf5File3Path, 'r')
        jump = hf.get('elementAperture')
        jump = np.array(jump)
        aperture_HydroFrac = jump[0,:]
        x = hf.get('elementAperture elementCenter')
        loc_HydroFrac = x[0,:,1]
  
        #-------- Extract info from XML
        xmlFilePath = "../../../../../../inputFiles/efemFractureMechanics/Sneddon_embeddedFrac_base.xml"
        
        mechanicalParameters = getMechanicalParametersFromXML(xmlFilePath)
        appliedPressure = getFracturePressureFromXML(xmlFilePath)

        # Get length of the fracture
        length, origin = getFractureLengthFromXML(xmlFilePath)
        
        # Initialize Sneddon's analytical solution
        sneddonAnalyticalSolution = Sneddon(mechanicalParameters, length, appliedPressure)

        # Plot analytical (continuous line) and numerical (markers) aperture solution
        x_analytical = np.linspace(-length, length, 501, endpoint=True)
        aperture_analytical = np.empty(len(x_analytical))
        i=0
        for xCell in x_analytical:
            aperture_analytical[i] = sneddonAnalyticalSolution.computeAperture( xCell )
            i += 1

        fsize = 30
        msize = 15
        lw=6
        fig, ax = plt.subplots(1,figsize=(16, 12))
        cmap = plt.get_cmap("tab10")
        
    
        N1=1
        ax.plot(x_analytical, aperture_analytical*1.0e3, color='k', label='Analytical Solution', lw=lw)
        ax.plot(loc_EmbeddeFrac[0::N1], aperture_EmbeddeFrac[0::N1]*1.0e3, 'o', color=cmap(0), label='Embedded Fracture', markersize=msize*0.6, alpha=0.8)
        ax.plot(loc_Contact[0::N1], aperture_Contact[0::N1]*1.0e3, 's', color=cmap(4), label='Lagrange Contact', markersize=msize*0.6, alpha=0.6)
        ax.plot(loc_HydroFrac[0::N1], aperture_HydroFrac[0::N1]*1.0e3, 'X', color=cmap(1), label='HydroFrac Solver', markersize=msize*0.6, alpha=0.4)    

        ax.grid()
        ax.set_xlabel('Fracture Length [m]', size=fsize, weight="bold")
        ax.set_ylabel('Fracture Aperture [mm]', size=fsize, weight="bold")
        ax.legend(bbox_to_anchor=(0.5, 0.2), loc='center', borderaxespad=0., fontsize=fsize)
        ax.xaxis.set_tick_params(labelsize=fsize)
        ax.yaxis.set_tick_params(labelsize=fsize)
        
        plt.show()

    if __name__ == "__main__":
        main()



------------------------------------------------------------------
To go further
------------------------------------------------------------------

**Feedback on this example**

This concludes the Sneddon example.
For any feedback on this example, please submit a `GitHub issue on the project's GitHub page <https://github.com/GEOSX/GEOSX/issues>`_.
