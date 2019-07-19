.. _Constitutive:

Constitutive Models
============================================

Constitutive models describe relations between various physical quantities.
In a physics simulation they are used to model the response or state of material (solid, fluid, or a mixture) as a function of input variables.

In GEOSX constitutive models are listed in ``<Constitutive>`` block of the input XML file.
Each type of model has its own XML tag and each model must be assigned a unique name via ``name`` attribute.
Names are used to assign models to regions of the physical domain via a ``materialList`` attribute of the ``<ElementRegion>`` node, see :ref:`XML_ElementRegion`.
In some cases, physics solvers must also be assigned specific constitutive models to use, see :ref:`Solvers`.

Model hierarchy
------------------------------------

There are several types of constitutive models that differ in purpose, input and output variables.
Currently supported constitutive models are:

* Solids

  * Solid mechanics

    * :doc:`/coreComponents/constitutive/docs/LinearElasticIsotropic`
    * :doc:`/coreComponents/constitutive/docs/LinearElasticAnisotropic`

* Fluids

  * Single phase fluids

    * :doc:`/coreComponents/constitutive/docs/CompressibleSinglePhaseFluid`

  * Multiphase fluids

    * :doc:`/coreComponents/constitutive/docs/BlackOilFluid`
    * :doc:`/coreComponents/constitutive/docs/CompositionalMultiphaseFluid`
    * :doc:`/coreComponents/constitutive/docs/CO2-EOS`

* Porous solids

  * Poroelasticity

    * :doc:`/coreComponents/constitutive/docs/PoroLinearElasticIsotropic`
    * :doc:`/coreComponents/constitutive/docs/PoroLinearElasticAnisotropic`

* Solid-fluid interaction

  * Relative permeability

    * :doc:`/coreComponents/constitutive/docs/BrooksCoreyRelativePermeability`
    * :doc:`/coreComponents/constitutive/docs/ThreePhaseRelativePermeability`

  * Capillary pressure

    * :doc:`/coreComponents/constitutive/docs/BrooksCoreyCapillaryPressure`
    * :doc:`/coreComponents/constitutive/docs/VanGenuchtenCapillaryPressure`

Input example
------------------------------------

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
      <ElementRegion name="Region2"
                     cellBlocks="cb1"
                     materialList="water shale"/>
    </ElementRegions>
    ...
  </Problem>
