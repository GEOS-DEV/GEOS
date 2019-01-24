=====================
Data sets examples
=====================

Introduction
===============

GEOSX comes with a suite of simple data sets to get you started. We are going
to explore a suite of simple sets with different geometrical objects. All
meshes will be imported from external files using PAMELA as the mesh importer.


Cube made of hexahedral elements
------------------------------------


This example consists of a simple sugar-cube stack of size 10x10x10.


Problem description
~~~~~~~~~~~~~~~~~~~~

We are going to propagate fluid from one vertical face of a cube to the
opposite side. The displacement is single phase, compressible, subject
to gravity forces.


Looking at the XML file
~~~~~~~~~~~~~~~~~~~~~~~~~

We are going to inspect blocks in the following XML file:
``CoreComponents\physicsSolvers\integratedTests\singlePhaseFlow\pamela_test\3D_10x10x10_compressible_pamela_hex_gravity.xml``

The file contains a number of XML blocks. We will describe the most important.

1. Solver specification

Here, we specify the type of solver (``singlePhaseTPFA``), the ``gravityVector``,
the tolerances of the Newton and Krylov iterations, and the maximum number
of iterations.

.. code-block:: XML

    <Solvers
      gravityVector="0.0,0.0,-9.81">
        <SinglePhaseFlow name="SinglePhaseFlow"
                            verboseLevel="0"
                            gravityFlag="1"
                            fluidName="water"
                            solidName="rock"
                            discretization="singlePhaseTPFA"
                            targetRegions="Domain">
        <SystemSolverParameters name="SystemSolverParameters"
                                krylovTol="1.0e-10"
                                newtonTol="1.0e-6"
                                maxIterNewton="8"/>
      </SinglePhaseFlow>
    </Solvers>


2. Mesh specification

Here, we specify the source of our mesh and some homogeneous properties to apply
one the grid blocks. In this example, the mesh is imported using PAMELA from
an existing file called ``cube_10x10x10_hex.msh``.

.. code-block:: XML

    <Mesh>
        <PAMELAMeshGenerator name="CubeHex"
          file="cube_10x10x10_hex.msh">
          <Property name="permx"
                    path="all"
                    component="0"
                    fieldName="permeability"
                    source="XML"
                    value="2.0e-16"/>

          <Property name="permy"
                    path="all"
                    component="1"
                    fieldName="permeability"
                    source="XML"
                    value="2.0e-16"/>

          <Property name="permz"
                    path="all"
                    component="2"
                    fieldName="permeability"
                    source="XML"
                    value="2.0e-16"/>

          <Property name="poro"
                    path="all"
                    fieldName="referencePorosity"
                    source="PAMELA"
                    nameInPAMELA="poro"/>

          <Property name="initialPressure"
                    path="all"
                    fieldName="pressure"
                    source="TXT"
                    filePath="allPressures.txt"/>
        </PAMELAMeshGenerator>
    </Mesh>



3. Events specification

.. code-block:: XML

    <Events maxTime="100">
      <!-- This event is applied every cycle, and overrides the
      solver time-step request -->
      <PeriodicEvent name="solverApplications"
                     forceDt="1"
                     target="/Solvers/SinglePhaseFlow" />

      <!-- This event is applied every 1.0s.  The targetExactTimestep
      flag allows this event to request a dt modification to match an
      integer multiple of the timeFrequency. -->
      <PeriodicEvent name="outputs"
                     timeFrequency="1"
                     targetExactTimestep="1"
                     target="/Outputs/siloWellPump" />

      <PeriodicEvent name="restarts"
                     timeFrequency="1e99"
                     targetExactTimestep="0"
                     target="/Outputs/sidreRestart"
                     endTime="-1"/>
    </Events>



4. Numerical methods

.. code-block:: XML

  <NumericalMethods>
    <FiniteVolume>
      <TwoPointFluxApproximation name="singlePhaseTPFA"
                                 fieldName="pressure"
                                 boundaryFieldName="facePressure"
                                 coefficientName="permeability"/>
    </FiniteVolume>
  </NumericalMethods>


5. Regions

.. code-block:: XML

    <ElementRegions>
      <ElementRegion name="Domain" cellBlocks="PART00001_POLYHEDRON_POLYHEDRON_GROUP_1_HEX" materialList="water rock"/>
    </ElementRegions>



6. Constitutive equations specifications

.. code-block:: XML

    <Constitutive>
      <CompressibleSinglePhaseFluid name="water"
                                    referencePressure="0.0"
                                    referenceDensity="1000"
                                    compressibility="1e-9"
                                    referenceViscosity="0.001"
                                    viscosibility="0.0"/>
      <PoreVolumeCompressibleSolid name="rock"
                                   referencePressure="0.0"
                                   compressibility="1e-9"/>
    </Constitutive>

7. Boundary pressure conditions

.. code-block:: XML

  <InitialConditions
    </InitialConditions>
  <FieldSpecifications>
    <FieldSpecification name="boundaryPressure"
               objectPath="faceManager"
               fieldName="facePressure"
               scale="1.1e3"
               setNames="left"/>
  </FieldSpecifications>

Done





Running GEOSX
~~~~~~~~~~~~~~~~~~~~~~~~~


Inspecting results
~~~~~~~~~~~~~~~~~~~~~~~~~



Cube made of tetrahedral elements
------------------------------------


Problem description
~~~~~~~~~~~~~~~~~~~~


Looking at the XML file
~~~~~~~~~~~~~~~~~~~~~~~~~


Running GEOSX
~~~~~~~~~~~~~~~~~~~~~~~~~


Inspecting results
~~~~~~~~~~~~~~~~~~~~~~~~~





Cube made of pyramidal elements
------------------------------------


Problem description
~~~~~~~~~~~~~~~~~~~~


Looking at the XML file
~~~~~~~~~~~~~~~~~~~~~~~~~


Running GEOSX
~~~~~~~~~~~~~~~~~~~~~~~~~


Inspecting results
~~~~~~~~~~~~~~~~~~~~~~~~~
