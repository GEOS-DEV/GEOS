.. _Constitutive:

Constitutive Models
============================================

Constitutive models describe relations between various physical quantities.
In a physics simulation they are used to model the response or state of material (solid, fluid, or a mixture) as a function of input variables.

In GEOSX constitutive models are listed in ``<Constitutive>`` block of the input XML file.
Each type of model has its own XML tag and each model must be assigned a unique name via ``name`` attribute.
Names are used to assign models to regions of the physical domain via a ``materialList`` attribute of the ``<ElementRegion>`` node, see :ref:`XML_ElementRegions`.
In some cases, physics solvers must also be assigned specific constitutive models to use, see :ref:`Solvers`.

Typical ``<Constitutive>`` and ``<ElementRegions>`` block will then look like:

.. code-block:: xml

  <Problem>
    ...
    <Constitutive>
      <PoroLinearElasticIsotropic name="shale"
                                  defaultDensity="2700"
                                  defaultBulkModulus="61.9e6"
                                  defaultShearModulus="28.57e6"
                                  BiotCoefficient="1.0"/>

      <CompressibleSinglePhaseFluid name="water"
                                    referencePressure="2.125e6"
                                    referenceDensity="1000"
                                    compressibility="1e-19"
                                    referenceViscosity="0.001"
                                    viscosibility="0.0"/>
    </Constitutive>

    <ElementRegions>
      <ElementRegion name="region1"
                     cellBlocks="cellBlock1"
                     materialList="water shale"/>
    </ElementRegions>
    ...
  </Problem>

There are several types of constitutive models that differ in their input and output variables.
Currently supported constitutive models are:

.. toctree::
   :maxdepth: 2

   SolidMechanicsModels

   FluidModels
 
   RelativePermeabilityModels
 
   CapillaryPressureModels

