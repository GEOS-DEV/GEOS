.. _LinearElasticIsotropic:

############################################
Linear elastic isotropic solid model
############################################

Overview
=========================

This model may be used to represents a solid material with a linear elastic isotropic response to loading.
The relationhip between stress and strain is typically represented by Hooke's Law, 
which for the case of linear elasticity, may be expressed as:

.. math::
   \sigma_{ij} + \lambda \epsilon_{kk} + 2 \mu \epsilon_{ij},
   
where :math:`sigma_{ij}` is the :math:`ij` component of the cauchy stress tensor, 
:math:`\epsilon_{ij}` is the :math:`ij` component of the the strain tensor,
:math:`\lambda` is the Lames elastic constant,
and :math:`\mu` is the elastic shear modulus.

Hooke's Law may also be expressed using Voigt notation for the stress and strain tensors as:

.. math::
   \tensor{\sigma} = \tensor{C} \dot \tensor{\epsilon}
   
   
Usage
=========================

The following attributes are supported:

.. include:: /coreComponents/fileIO/schema/docs/LinearElasticIsotropic.rst

Input example
=========================

.. code-block:: xml

  <Constitutive>
    <LinearElasticIsotropic name="shale"
                            defaultDensity="2700"
                            defaultBulkModulus="61.9e6"
                            defaultShearModulus="28.57e6"
  </Constitutive>