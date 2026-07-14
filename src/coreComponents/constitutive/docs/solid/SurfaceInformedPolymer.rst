.. _SurfaceInformedPolymer:

############################################
Model: Surface Informed Polymer
############################################

Overview
=========================

``SurfaceInformedPolymer`` is a corotational explicit-MPM polymer model for
large-deformation continuum particles.  It follows the integration philosophy of
``StrainHardeningPolymer``: the MPM solver supplies a corotational strain
increment, the material performs an elastic-plastic stress update, finite
rotation is handled outside the constitutive model, and no implicit tangent is
provided.  The new form is intended to use the same temperature, crystallinity,
softening, hardening, and pressure-asymmetry functions as the finite-thickness
cohesive projection described in :ref:`SurfaceInformedPolymerCohesiveZone`.

The model is designed for polymer response in which the flow stress is governed
by an initial yield strength, a decaying post-yield softening contribution, and a
large-strain chain-stretch hardening contribution.  The yield envelope can also
include a weak pressure-sensitive tension/compression asymmetry.  The pressure
term is non-associated: it changes the scalar yield strength, while the return
mapping remains radial in deviatoric stress and retains the trial mean stress.

Temperature and crystallinity scaling
=====================================

All reference material parameters are interpreted at the transition/reference
temperature ``glassTransitionTemperature``.  The scalar temperature scale is
normalized so that

.. math::

   S_T(T_g)=1.

The implementation evaluates the scale in logarithmic space,

.. math::

   \log S_T =
      a_c \max(T_g-T,0)
      - a_h \max(T-T_g,0)
      - J_T\left[
         \frac{1}{1+\exp[-(T-T_g)/w_T]}-\frac{1}{2}
       \right],

where ``temperatureColdSlope`` is :math:`a_c`, ``temperatureHotSlope`` is
:math:`a_h`, ``temperatureTransitionMagnitude`` is :math:`J_T`, and
``temperatureTransitionWidth`` is :math:`w_T`.  A zero transition magnitude
reduces the expression to a two-sided log-linear temperature scale.

The optional crystallinity correction uses a smooth activation function,

.. math::

   A_X(T)=\frac{1}{1+\exp[-(T-T_g)/w_X]},

where :math:`w_X` is ``crystallinityTransitionWidth``.  The elastic and yield
multipliers are

.. math::

   C_E=1+\beta_E\,(X_c-X_{c,0})A_X(T),
   \qquad
   C_Y=1+\beta_Y\,(X_c-X_{c,0})A_X(T),

with lower bounds imposed internally to keep the scales positive.  Setting the
crystallinity coefficients to zero disables this effect.

The material functions used by the continuum update are

.. math::

   K(T,X_c) = K_g S_T C_E,
   \qquad
   G(T,X_c) = G_g S_T C_E,

.. math::

   Y(T,X_c) = Y_g S_T C_Y,
   \qquad
   S_{soft}(T) = S_{soft,g} S_T,

.. math::

   H(T) = H_g S_T^{p_H},

where :math:`p_H` is ``hardeningScaleExponent``.  This separate hardening
exponent allows post-yield hardening to have a weaker or stronger temperature
sensitivity than yield and softening.

Flow surface
============

The scalar flow surface is

.. math::

   \Phi = q + \eta(T) p_{eff} - \sigma_f^0 \le 0,

where :math:`q` is the von-Mises equivalent stress and :math:`p_{eff}` is the
mean stress used only by the pressure-asymmetry term.  The base flow strength is

.. math::

   \sigma_f^0 =
      Y(T,X_c)
      + S_{soft}(T) \exp\left[-\left(\frac{\kappa}{r_1}\right)^{r_2}\right]
      + H(T) h_\lambda,

with

.. math::

   h_\lambda = \lambda_h^2 - \lambda_h^{-1},
   \qquad
   \lambda_h = \max(\lambda_{chain},1).

Here :math:`\kappa` is the accumulated plastic-strain-like history variable,
``shearSofteningShapeParameter1`` is :math:`r_1`, and
``shearSofteningShapeParameter2`` is :math:`r_2`.

The pressure-asymmetry coefficient is localized around the reference
temperature,

.. math::

   \eta(T) = \eta_0
      \exp\left[-\frac{1}{2}\left(\frac{T-T_g}{w_\eta}\right)^2\right],

where ``pressureAsymmetryAmplitude`` is :math:`\eta_0` and
``pressureAsymmetryWidth`` is :math:`w_\eta`.

Compressive pressure cap
========================

The pressure-asymmetry term is a low-pressure strength correction.  Applying it
linearly at very large confined compression can make the shear strength grow
without bound.  The optional input ``compressivePressureStrengtheningCap`` clips
only the compressive side of the mean stress used in the pressure term,

.. math::

   p_{eff}=\max(p,-p_c),

where :math:`p_c` is the cap magnitude.  A negative cap disables clipping and
recovers the uncapped form.  The cap does not clip or modify the retained
volumetric stress; it only bounds pressure-assisted strengthening in the scalar
yield condition.

Chain stretch and failure
=========================

The chain-stretch measure is based on the right stretch tensor computed from the
particle deformation gradient.  If :math:`\lambda_{max}` is the maximum principal
stretch and :math:`J=\det F`, the implementation uses

.. math::

   \lambda_{chain}=\max\left(\lambda_{max}, J^{-1/3}\lambda_{max}\right).

This preserves ordinary tensile stretch hardening, avoids artificial chain
extension in hydrostatic compression, and still allows constrained compression or
shear under pressure to reach finite chain extensibility.

By default, the material damage flag is set to one when
``maximumStretch`` is exceeded.  An optional temperature-dependent fracture
stretch can be enabled by setting ``fractureStretchLambda0`` to a positive
value:

.. math::

   \lambda_f(T)=\lambda_{min}+\lambda_0
   \exp\left(\frac{T-T_0}{a}\right).

The corresponding input names are ``fractureStretchLambdaMin``,
``fractureStretchLambda0``, ``fractureStretchT0``, and
``fractureStretchTemperatureScale``.  When the optional law is active,
``lambda_f`` replaces ``maximumStretch`` as the damage trigger.  This
maximum-stretch failure is not suppressed by pressure.


Primary XML attributes
======================

The continuum card inherits ``defaultDensity``, ``defaultBulkModulus``, and
``defaultShearModulus`` from ``ElasticIsotropic``.  The additional polymer inputs
are grouped below.

.. list-table:: SurfaceInformedPolymer inputs
   :header-rows: 1

   * - Attribute
     - Purpose
   * - ``defaultYieldStrength``
     - Reference yield strength :math:`Y_g`.
   * - ``shearSofteningMagnitude``
     - Reference magnitude :math:`S_{soft,g}` of the decaying softening term.
   * - ``shearSofteningShapeParameter1``
     - Plastic-strain scale :math:`r_1` for softening.
   * - ``shearSofteningShapeParameter2``
     - Shape exponent :math:`r_2` for softening.
   * - ``strainHardeningSlope``
     - Reference stretch-hardening slope :math:`H_g`.
   * - ``hardeningScaleExponent``
     - Exponent :math:`p_H` in :math:`H=H_gS_T^{p_H}`.
   * - ``maximumStretch``
     - Constant chain stretch at which the damage flag is set to one when the optional exponential fracture-stretch law is disabled.
   * - ``fractureStretchLambdaMin``, ``fractureStretchLambda0``
     - Optional :math:`\lambda_{min}` and :math:`\lambda_0` in :math:`\lambda_f=\lambda_{min}+\lambda_0\exp((T-T_0)/a)`.  A non-positive ``fractureStretchLambda0`` disables the law.
   * - ``fractureStretchT0``, ``fractureStretchTemperatureScale``
     - Optional reference temperature :math:`T_0` and temperature scale :math:`a` for the exponential fracture-stretch law.
   * - ``glassTransitionTemperature``
     - Reference temperature :math:`T_g` for the normalized scale.
   * - ``temperatureColdSlope``, ``temperatureHotSlope``
     - Slopes of :math:`\log S_T` below and above :math:`T_g`.
   * - ``temperatureTransitionMagnitude``, ``temperatureTransitionWidth``
     - Centered smooth transition in :math:`\log S_T`.
   * - ``crystallinity``, ``referenceCrystallinity``
     - Current and reference crystallinity measures.
   * - ``crystallinityTransitionWidth``
     - Temperature width for crystallinity activation.
   * - ``elasticCrystallinityCoeff``
     - Linear crystallinity coefficient for elastic moduli.
   * - ``yieldStrengthCrystallinityCoeff``
     - Linear crystallinity coefficient for yield strength.
   * - ``pressureAsymmetryAmplitude``
     - Amplitude :math:`\eta_0` of the pressure-sensitive strength asymmetry.
   * - ``pressureAsymmetryWidth``
     - Temperature width :math:`w_\eta` for pressure asymmetry.
   * - ``compressivePressureStrengtheningCap``
     - Optional compressive cap :math:`p_c` applied only to the pressure-strengthening term.

Output fields
=============

The model stores plastic strain in engineering Voigt form and also stores a
scalar equivalent-plastic-strain-like history variable.  To make primitive
particle data available in Silo/VisIt, request the relevant GEOS fields through
the MPM solver's ``plottableFields`` input, for example
``particleStress``, ``particlePlasticStrain``, and ``particleDamage``.  The
``particleFileFields`` PFW input controls initial particle-file columns and does
not by itself request Silo output fields.

Example
~~~~~~~

A typical XML block is

.. code-block:: xml

  <Constitutive>
    <SurfaceInformedPolymer
      name="polymer"
      defaultDensity="1.85"
      defaultBulkModulus="260.0"
      defaultShearModulus="5.0"
      defaultYieldStrength="7.0"
      glassTransitionTemperature="300.0"
      temperatureColdSlope="0.02"
      temperatureHotSlope="0.04"
      temperatureTransitionMagnitude="1.0"
      temperatureTransitionWidth="5.0"
      referenceCrystallinity="0.0"
      crystallinity="0.0"
      elasticCrystallinityCoeff="0.0"
      yieldStrengthCrystallinityCoeff="0.0"
      shearSofteningMagnitude="2.0"
      shearSofteningShapeParameter1="0.2"
      shearSofteningShapeParameter2="2.0"
      strainHardeningSlope="2.0"
      hardeningScaleExponent="0.85"
      pressureAsymmetryAmplitude="0.0"
      pressureAsymmetryWidth="4.0"
      compressivePressureStrengtheningCap="26.0"
      maximumStretch="3.0"
      fractureStretchLambdaMin="1.0"
      fractureStretchLambda0="0.0"
      fractureStretchT0="300.0"
      fractureStretchTemperatureScale="1.0" />
  </Constitutive>

Notes and limitations
=====================

``SurfaceInformedPolymer`` is an explicit-MPM constitutive model.  It does not
provide an implicit tangent and should not be used as a drop-in implicit solid
model without additional development.  The compressive pressure cap is a
bounded-extrapolation control for the pressure-asymmetry term; predictive
response at very high pressure also requires an equation of state and pressure
hardening calibration appropriate to the material and loading path.
