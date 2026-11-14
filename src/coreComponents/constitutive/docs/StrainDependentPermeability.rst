.. _StrainDependentPermeability:


####################################################
Strain Dependent Permeability Model
####################################################


Overview
======================

In the input file, this model is selected with the ``<StrainDependentPermeability>`` XML tag.

The strain dependent permeability model relates the permeability of the porous rock to its volumetric
strain, following `Chin et al. (2000)`_. Assuming that the volume of the solid grains is conserved
(incompressible grains), the porosity is a function of the volumetric strain:

.. math::
   \phi( \epsilon_v ) = 1 - ( 1 - \phi_0 ) \, e^{ -\epsilon_v }

and each diagonal component of the permeability tensor is updated with a power law of the porosity ratio:

.. math::
   \frac{ k_i }{ k_{0,i} } = \left[ \frac{ 1 - ( 1 - \phi_0 ) \, e^{ -\epsilon_v } }{ \phi_0 } \right]^{ n_i }, \quad i \in \{ x, y, z \}

where:

- :math:`\epsilon_v = \text{tr}( \boldsymbol{\epsilon} )` is the total volumetric strain since the initial state,
- :math:`\phi_0` is the initial porosity, provided by the porosity model,
- :math:`k_{0,i}` is the permeability at the initial state (zero volumetric strain),
- :math:`n_i` is the exponent of the porosity ratio for the permeability component :math:`i`.

The derivative of the permeability with respect to the volumetric strain is also computed:

.. math::
   \frac{ \partial k_i }{ \partial \epsilon_v } = n_i \, k_i \, \frac{ 1 - \phi }{ \phi }


Sign convention
======================

GEOS uses the convention of positive strains in extension. Therefore:

- :math:`\epsilon_v > 0` (dilation) increases the porosity and the permeability,
- :math:`\epsilon_v < 0` (compaction) decreases the porosity and the permeability.

The volumetric strain is measured from the initial state, i.e., after the stress initialization
(the strain of the stress initialization does not modify the permeability). At the initial state,
:math:`\epsilon_v = 0` and the permeability is equal to :math:`k_0`.


Remarks
======================

- This model must be used in a poromechanical simulation, through a ``Porous<Solid>StrainDependentPermeability``
  model (for instance ``PorousElasticIsotropicStrainDependentPermeability``), with either the fully coupled or the
  sequential (fixed-stress) scheme.
- The permeability is cell-centered: it is computed from the volumetric strain averaged over the quadrature points of
  the element.
- If the permeability is not set by a field specification, it is initialized with ``referencePermeabilityComponents``.
  Otherwise, the initial permeability field is used as :math:`k_0`.
- Under strong compaction (:math:`\epsilon_v < \ln( 1 - \phi_0 )`), the computed porosity would become negative:
  it is clamped to a small positive value, and the permeability no longer changes with the volumetric strain.


Parameters
======================

The Strain Dependent Permeability model can be called in the
``<Constitutive>`` block of the input XML file.
This permeability model must be assigned a unique name via the
``name`` attribute.
This name is used to assign the model to regions of the physical
domain via a ``permeabilityModelName`` attribute in the ``<Porous<Solid>StrainDependentPermeability>``
block.

The following attributes are supported:

.. include:: /docs/sphinx/datastructure/StrainDependentPermeability.rst


Example
=======================

.. code-block:: xml

   <Constitutive>
      ...
      <PorousElasticIsotropicStrainDependentPermeability
        name="porousRock"
        solidModelName="rock"
        porosityModelName="rockPorosity"
        permeabilityModelName="rockPerm"/>

      <StrainDependentPermeability
        name="rockPerm"
        referencePermeabilityComponents="{ 1.0e-15, 1.0e-15, 1.0e-15 }"
        permeabilityExponents="{ 3.0, 3.0, 3.0 }"/>
      ...
   </Constitutive>

.. _Chin et al. (2000): https://doi.org/10.2118/58968-PA
