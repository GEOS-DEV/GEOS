.. _PorousSolids:

############################################
Porous Solids
############################################

Overview
========================

CompressibleSolid
========================


.. code-block:: xml

   <Constitutive>
     <CompressibleSolidConstantPermeability name="porousRock"
                                           solidModelName="nullSolid"
                                           porosityModelName="rockPorosity"
                                           permeabilityModelName="rockPermeability"/>
   </Constitutive>


PorousSolid
======================
To run poromechanical problems, the total stress is decomposed into an "effective stress" (driven by mechanical deformations) and a pore fluid
pressure component, following the `Biot theory of poroelasticity <https://doi.org/10.1016/B978-0-08-040615-2.50011-3>`__.
For single-phase flow, or multiphase problems with no capillarity, this decomposition reads

.. math::
   \sigma_{ij} = \sigma\prime_{ij}  - b p \delta_{ij}

where :math:`\sigma_{ij}` is the :math:`ij` component of the total stress tensor,
:math:`\sigma\prime_{ij}` is the :math:`ij` component of the effective (Cauchy) stress tensor,
:math:`b` is Biot's coefficient,
:math:`p` is fluid pressure,
and :math:`\delta` is the Kronecker delta.

.. code-block:: xml

    <Constitutive>
      <PorousElasticIsotropic name="porousRock"
                              porosityModelName="rockPorosity"
                              solidModelName="rockSkeleton"
                              permeabilityModelName="rockPermeability"/>
    </Constitutive>

Note that any of the previously described solid models may then be used to compute the effective stress, leading to either
poro-elastic, poro-plastic, or poro-damage behavior depending on the specific model chose.
