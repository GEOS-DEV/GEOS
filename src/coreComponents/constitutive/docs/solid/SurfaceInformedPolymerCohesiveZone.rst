.. _SurfaceInformedPolymerCohesiveZone:

############################################################
Model: Surface Informed Polymer Cohesive Zone
############################################################

Overview
=========================

``SurfaceInformedPolymerCohesiveZone`` is a finite-thickness cohesive-zone
projection of the :ref:`SurfaceInformedPolymer` continuum model.  It is intended
for cohesive interfaces that represent a thin polymer layer rather than a
zero-thickness brittle interface.  The cohesive law uses the same reference
thermal scale, crystallinity multipliers, softening law, stretch-hardening law,
pressure-asymmetry function, compressive pressure cap, and maximum-stretch
failure criterion as the continuum model.

The cohesive model operates on normal and tangential displacement jumps, but
interprets those jumps as nominal strains in a film of thickness ``thickness``.
This allows the cohesive-zone parameters to be derived from the corresponding
continuum polymer parameters and a layer thickness.

Film kinematics
===============

For normal jump :math:`\delta_n`, tangential jump magnitude
:math:`\delta_t`, and film thickness :math:`h_0`, the nominal film strains are

.. math::

   \epsilon_n = \frac{\delta_n}{h_0},
   \qquad
   \gamma = \frac{\delta_t}{h_0}.

For the chain-stretch calculation, the implementation constructs a reduced film
deformation gradient,

.. math::

   F_{film} =
   \begin{bmatrix}
      1 & \gamma & 0 \\
      0 & 1+\epsilon_n & 0 \\
      0 & 0 & 1
   \end{bmatrix},

and evaluates the same chain-stretch measure used by the continuum model,

.. math::

   \lambda_{chain}=\max\left(\lambda_{max}, J^{-1/3}\lambda_{max}\right).

If ``maximumStretch`` is exceeded, the cohesive damage flag is set to one and
zero traction is returned.  This finite-extensibility failure is pressure
independent.

Elastic split and flow surface
==============================

The normal response is split into a retained volumetric part and a deviatoric
normal part.  With temperature- and crystallinity-scaled moduli :math:`K` and
:math:`G`,

.. math::

   p = K\epsilon_n,
   \qquad
   s_n^{tr} = \frac{4G}{3}(\epsilon_n-\epsilon_n^p),
   \qquad
   \tau^{tr}=G(\gamma-\gamma^p).

The trial equivalent deviatoric stress is

.. math::

   q^{tr}=\sqrt{\left(\frac{3s_n^{tr}}{2}\right)^2+3(\tau^{tr})^2}.

The cohesive return is made to the same scalar surface as the continuum model,

.. math::

   q + \eta(T)p_{eff} \le \sigma_f^0,

where

.. math::

   \sigma_f^0 =
      Y(T,X_c)
      + S_{soft}(T)\exp\left[-\left(\frac{\kappa}{r_1}\right)^{r_2}\right]
      + H(T)h_\lambda.

Only the deviatoric normal/shear components are returned.  The volumetric film
stress :math:`p` is retained so that a nearly incompressible finite-thickness
polymer layer remains stiff in constrained normal loading.

Compressive pressure cap
========================

The cohesive model accepts the same ``compressivePressureStrengtheningCap`` as
the continuum model.  The cap applies only to the pressure-asymmetry term,

.. math::

   p_{eff}=\max(p,-p_c).

It does not cap the normal traction carried by the film and does not suppress
maximum-stretch failure.  For consistency, use the same cap in the continuum and
cohesive cards when the cohesive-zone law represents a thin layer of the same
polymer.


Primary XML attributes
======================

The cohesive card uses the same polymer parameters as the continuum card, but
its elastic inputs are named ``bulkModulus`` and ``shearModulus`` because the
cohesive law is not an ``ElasticIsotropic`` solid model.  Use the same reference
units and temperature convention as the corresponding continuum card.

.. list-table:: SurfaceInformedPolymerCohesiveZone inputs
   :header-rows: 1

   * - Attribute
     - Purpose
   * - ``thickness``
     - Physical film thickness :math:`h_0` used to convert jumps to strains.
   * - ``bulkModulus``, ``shearModulus``
     - Reference film moduli at ``glassTransitionTemperature``.
   * - ``defaultYieldStrength``
     - Reference yield strength :math:`Y_g`.
   * - ``shearSofteningMagnitude``
     - Reference magnitude of the decaying softening term.
   * - ``shearSofteningShapeParameter1``
     - Plastic-strain-like scale :math:`r_1` for softening.
   * - ``shearSofteningShapeParameter2``
     - Shape exponent :math:`r_2` for softening.
   * - ``strainHardeningSlope``
     - Reference stretch-hardening slope :math:`H_g`.
   * - ``hardeningScaleExponent``
     - Exponent :math:`p_H` in the hardening temperature scale.
   * - ``maximumStretch``
     - Constant chain stretch at which the cohesive law fails when the optional exponential fracture-stretch law is disabled.
   * - ``fractureStretchLambdaMin``, ``fractureStretchLambda0``
     - Optional :math:`\lambda_{min}` and :math:`\lambda_0` in :math:`\lambda_f=\lambda_{min}+\lambda_0\exp((T-T_0)/a)`.
   * - ``fractureStretchT0``, ``fractureStretchTemperatureScale``
     - Optional reference temperature :math:`T_0` and temperature scale :math:`a` for the exponential fracture-stretch law.
   * - ``glassTransitionTemperature``
     - Reference temperature :math:`T_g`.
   * - ``temperatureColdSlope``, ``temperatureHotSlope``
     - Slopes of :math:`\log S_T` below and above :math:`T_g`.
   * - ``temperatureTransitionMagnitude``, ``temperatureTransitionWidth``
     - Centered smooth transition in :math:`\log S_T`.
   * - ``crystallinity``, ``referenceCrystallinity``
     - Current and reference crystallinity measures.
   * - ``elasticCrystallinityCoeff``
     - Linear crystallinity coefficient for the film moduli.
   * - ``yieldStrengthCrystallinityCoeff``
     - Linear crystallinity coefficient for yield strength.
   * - ``pressureAsymmetryAmplitude``, ``pressureAsymmetryWidth``
     - Magnitude and temperature width of pressure-sensitive asymmetry.
   * - ``compressivePressureStrengtheningCap``
     - Optional compressive cap applied only to pressure strengthening.

Sign convention
===============

The internal film stress is tensile-positive.  The GEOS cohesive normal-traction
convention is opposite for opening, so a tensile film normal stress is returned
as a negative cohesive normal stress.  The shear traction is returned with the
corresponding cohesive-zone sign convention used by the MPM cohesive-force
assembler.

Example
~~~~~~~

A typical XML block is

.. code-block:: xml

  <Constitutive>
    <SurfaceInformedPolymerCohesiveZone
      name="polymerCZ"
      thickness="0.1"
      bulkModulus="260.0"
      shearModulus="5.0"
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

The cohesive law is attached to an MPM interface with a ``CohesiveZoneRegion``
and a ``ReferenceCohesiveZones`` event.  The particle input must provide the
surface flag, surface normal, surface position, and cohesive-zone tag fields
needed by the reference-configuration cohesive-zone initialization.

Notes and limitations
=====================

This model is a reduced thin-film approximation, not a general three-dimensional
continuum replacement.  Its pressure and volumetric response are based on the
normal film strain and are intended to match a constrained layer response.  Use a
continuum polymer region when the polymer thickness is large enough to resolve
with particles and a cohesive-zone region when the polymer is best represented as
a finite-thickness interface.
