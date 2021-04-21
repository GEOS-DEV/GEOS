.. _PoroLinearElasticIsotropic:

############################################
Linear poroelastic isotropic solid model
############################################

Overview
=========================

This model may be used to represents a porous material with a linear poroelastic isotropic response to coupled deformation-diffusion processes.
The relationship between stress and strain is typically formulated within the framework of the `Biot theory of poroelasticity <https://doi.org/10.1016/B978-0-08-040615-2.50011-3>`__,
which for the case of linear poroelasticity, may be expressed as:

.. math::
   \sigma_{ij} = \sigma\prime_{ij}  - b p \delta_{ij} = \lambda \epsilon_{kk} + 2 \mu \epsilon_{ij} - b p \delta_{ij},

where :math:`\sigma_{ij}` is the :math:`ij` component of the total stress tensor,
:math:`\sigma\prime_{ij}` is the :math:`ij` component of the effective (Cauchy) stress tensor,
:math:`\epsilon_{ij}` is the :math:`ij` component of the the strain tensor,
:math:`\lambda` is the Lames elastic constant,
:math:`\mu` is the elastic shear modulus,
:math:`b` is Biot's coefficient,
:math:`p` is fluid pressure,
and :math:`\delta` is Kronecker delta.

Hooke's Law may also be expressed using `Voigt notation <https://en.wikipedia.org/wiki/Voigt_notation>`__ for the effective stress and strain tensors as:

.. math::
   \tensor{\sigma\prime} = \tensor{C} \cdot \tensor{\epsilon},

or,

.. math::
    \begin{bmatrix}
      \sigma\prime_{11} \\
      \sigma\prime_{22} \\
      \sigma\prime_{33} \\
      \sigma\prime_{23} \\
      \sigma\prime_{13} \\
      \sigma\prime_{12}
    \end{bmatrix}
    =
    \begin{bmatrix}
      2\mu+\lambda  &   \lambda     &   \lambda   & 0   & 0 & 0 \\
          \lambda     &  2\mu+\lambda   &   \lambda   & 0   & 0 & 0 \\
          \lambda     &    \lambda    & 2\mu+\lambda & 0  & 0 & 0 \\
          0         &       0     &       0 &\mu  & 0 & 0 \\
      0         &           0     & 0       & 0   & \mu & 0 \\
      0         &       0     & 0       & 0   & 0 & \mu
    \end{bmatrix}
    \begin{bmatrix}
      \epsilon_{11} \\
      \epsilon_{22} \\
      \epsilon_{33} \\
      2\epsilon_{23} \\
      2\epsilon_{13} \\
      2\epsilon_{12}
    \end{bmatrix}.

Variations
==========
The application of linear poroelasticity as presented above is typically restricted to the case of
`infinitesimal strain <https://en.wikipedia.org/wiki/Infinitesimal_strain_theory>`__.
For the case of `infinitesimal strain <https://en.wikipedia.org/wiki/Infinitesimal_strain_theory>`__, the
above relations are applied directly.
For the case of `finite strain <https://en.wikipedia.org/wiki/Finite_strain_theory>`__,
the above relations may be slightly modified to an incremental update and rotation:

.. math::
   \Delta \tensor{\sigma\prime} &= \tensor{C} \cdot \hat{\tensor{D}},\\
   \tensor{\sigma\prime}^{n+1} &= \hat{\tensor{R}}( \tensor{\sigma\prime}^{n} + \Delta \tensor{\sigma\prime} ) \hat{\tensor{R}}^T,

where :math:`\hat{\tensor{D}}` is the "incremental rate of deformation tensor" and :math:`\hat{\tensor{R}}` is the incremental rotation tensor, which are
typically calculated from the `velocity gradient <https://en.wikipedia.org/wiki/Strain-rate_tensor>`__.
This extension into finite strain constitutes a `hypo-elastic update <https://en.wikipedia.org/wiki/Hypoelastic_material>`__,
and the choice of method to calculate :math:`\hat{\tensor{D}}`, and :math:`\hat{\tensor{R}}` determines if the update is objective.
One commonly used method is the `Hughes-Winget <https://onlinelibrary.wiley.com/doi/abs/10.1002/nme.1620151210>`__ algorithm.


Parameters
=========================

The following attributes are supported:

.. include:: /coreComponents/schema/docs/PoroElasticIsotropic.rst

Example
=========================

.. code-block:: xml

  <Constitutive>
    <PoroLinearElasticIsotropic name="shale"
                                defaultDensity="2700"
                                defaultBulkModulus="61.9e6"
                                defaultShearModulus="28.57e6"
                                BiotCoefficient="1.0"/>
  </Constitutive>
