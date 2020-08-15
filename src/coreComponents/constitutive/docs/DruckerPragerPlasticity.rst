.. _DruckerPragerPlasticity:

############################################
Drucker-Prager plasticity solid model
############################################

Overview
=========================

This model may be used to represents a solid material with plastic response to loading according to the `Drucker-Prager <https://en.wikipedia.org/wiki/Drucker%E2%80%93Prager_yield_criterion>`__ yield criterion below:


The material behavior initially starts as linear elastic (see :ref:`LinearElasticIsotropic`)
The relationship between stress and strain is typically represented by  
which for the case of linear elasticity, may be expressed as:

.. math::
   \sigma_{ij} = C_{ijkl} \epsilon_{kl},
   
where :math:`\sigma_{ij}` is the :math:`ij` component of the Cauchy stress tensor, 
:math:`\epsilon_{kl}` is the :math:`kl` component of the the strain tensor, and
:math:`C_{ijkl}` is the fourth order  elastic stiffness tensor.

Hooke's Law may also be expressed using `Voigt notation <https://en.wikipedia.org/wiki/Voigt_notation>`__ for the stress and strain tensors as:

.. math::
   \tensor{\sigma} = \tensor{C} \cdot \tensor{\epsilon},
   
or,

.. math::
    \begin{bmatrix}
      \sigma_{11} \\
      \sigma_{22} \\
      \sigma_{33} \\
      \sigma_{23} \\
      \sigma_{13} \\
      \sigma_{12}
    \end{bmatrix}
    = 
    \begin{bmatrix}
      c_{1111} & c_{1122} & c_{1133} & c_{1123} & c_{1113} & c_{1112} \\
      c_{2211} & c_{2222} & c_{2233} & c_{2223} & c_{2213} & c_{2212} \\
      c_{3311} & c_{3322} & c_{3333} & c_{3323} & c_{3313} & c_{3312} \\
      c_{2311} & c_{2322} & c_{2333} & c_{2323} & c_{2313} & c_{2312} \\
      c_{1311} & c_{1322} & c_{1333} & c_{1323} & c_{1313} & c_{1312} \\
      c_{2311} & c_{2322} & c_{2333} & c_{2323} & c_{2313} & c_{2312}
    \end{bmatrix}
    \begin{bmatrix}
      \epsilon_{11} \\
      \epsilon_{22} \\
      \epsilon_{33} \\
      2\epsilon_{23} \\
      2\epsilon_{13} \\
      2\epsilon_{12}
    \end{bmatrix}.

.. _druckerPragerYield:
.. figure:: ../../../coreComponents/constitutive/docs/DPyield.png
   :align: center
   :width: 500
   :figclass: align-center
   
Return Mapping Algorithm 
========================
The application of linear elasticity as presented above is typically restricted to the case of 
`infinitesimal strain <https://en.wikipedia.org/wiki/Infinitesimal_strain_theory>`__.
For the case of `infinitesimal strain <https://en.wikipedia.org/wiki/Infinitesimal_strain_theory>`__, the 
above relations are applied directly. 
For the case of `finite strain <https://en.wikipedia.org/wiki/Finite_strain_theory>`__,
the above relations may be slightly modified to an incremental update and rotation:

.. math::
   \Delta \tensor{\sigma} &= \tensor{C} \cdot \hat{\tensor{D}},\\
   \tensor{\sigma}^{n+1} &= \hat{\tensor{R}}( \tensor{\sigma}^{n} + \Delta \tensor{\sigma} ) \hat{\tensor{R}}^T,
   
where :math:`\hat{\tensor{D}}` is the "incremental rate of deformation tensor" and :math:`\hat{\tensor{R}}` is the incremental rotation tensor, which are 
typically calculated from the `velocity gradient <https://en.wikipedia.org/wiki/Strain-rate_tensor>`__.
This extension into finite strain constitutes a `hypo-elastic update <https://en.wikipedia.org/wiki/Hypoelastic_material>`__,
and the choice of method to calculate :math:`\hat{\tensor{D}}`, and :math:`\hat{\tensor{R}}` determines if the update is objective.
One commonly used method is the `Hughes-Winget <https://onlinelibrary.wiley.com/doi/abs/10.1002/nme.1620151210>`__ algorithm. 


Usage
=========================

The following attributes are supported:

.. include:: /coreComponents/fileIO/schema/docs/DruckerPrager.rst

Input example
=========================

.. code-block:: xml

  <Constitutive>
    <DruckerPrager name="shale"
                              defaultDensity="2700"
                              defaultBulkModulus="1.0"
                              defaultShearModulus="1.0"
                              defaultTanFrictionAngle="0.5"
                              defaultTanDilationAngle="0.5"
                              defaultHardeningRate="0.0"
                              defaultCohesion="1.0" />
  </Constitutive>
