.. _PoroElastic:

############################################
Porous Solids
############################################

Overview
----------------------------------------

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

Note that any of the previously described solid models may then be used to compute the effective stress, leading to either
poro-elastic, poro-plastic, or poro-damage behavior depending on the specific model chose. 

.. note::
   The above effective stress decomposition may also be modified to include a tensorial Biot modulus in the case of anisotropic
   materials, or a saturation-weighted effective stress decomposition in the case of multiphase flow with capillarity.

Parameters
-----------------------

The porous solid models simply append the keyword ``Poro`` in front of the
solid model they derive from: ``PoroElasticIsotropic``, ``PoroDruckerPrager``, and so on.

As an example, the following attributes are supported for a ``PoroElasticIsotropic`` model.  
Note that most of the parameters derive from the underlying solid model,
with a few additional parameters added for the Biot parameters. 

.. include:: /coreComponents/schema/docs/PoroElasticIsotropic.rst

The parameters for a ``PoroDruckerPrager`` model are:

.. include:: /coreComponents/schema/docs/PoroDruckerPrager.rst

Other models will be similar.

Example
-----------------------

.. code-block:: xml

  <Constitutive>
    <PoroElasticIsotropic name="shale"
                          defaultDensity="2700"
                          defaultBulkModulus="61.9e6"
                          defaultShearModulus="28.57e6"
                          BiotCoefficient="1.0"/>
  </Constitutive>


