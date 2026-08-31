.. _PhaseFieldFractureSolver:

##############################
Phase Field Fracture Solver
##############################

Introduction
===========================================

This section describes the phase-field fracture model implemented in GEOS. In a phase-field
description of fracture, a sharp crack surface is replaced by a smooth damage field, so that
crack nucleation, propagation, branching and merging are obtained from the solution of a
boundary value problem, without any explicit tracking of the crack geometry and without a
propagation criterion.

The model is solved as a coupled problem, and therefore involves three solver blocks:

* a coupling solver of type ``PhaseFieldFracture``, which drives the two single-physics solvers
  in a staggered fashion,
* a mechanics solver of type ``SolidMechanicsLagrangianFEM`` (see
  :ref:`SolidMechanicsLagrangianFEM`), which solves the quasi-static momentum balance using the
  degraded stiffness,
* a damage solver of type ``PhaseFieldDamageFEM``, which solves the phase-field evolution
  equation.

The behavior of the model is controlled by the constitutive model attached to the target
regions, which selects the type of fracture, the local dissipation function and the
tension-compression asymmetry. A poromechanical variant, ``PhaseFieldPoromechanics``, couples
the same damage solver with a single-phase poromechanics solver.


Theory
=========================

Regularized variational formulation
--------------------------------------

The model derives from the variational approach to fracture, in which the total energy of a
cracked body is the sum of a stored elastic energy and of the energy required to create the
crack surface. The sharp surface :math:`\Gamma` is regularized over a band of characteristic
width :math:`L`, described by a damage field :math:`d(x,t)`, with :math:`d = 0` in the intact
material and :math:`d = 1` in the fully broken material. The regularized energy reads

.. math::
   E(u, d) = \int_{\Omega} \left[ g(d) \, \psi^{+}(\varepsilon) + \psi^{-}(\varepsilon) \right] d\Omega
           + \frac{G_c}{c_0} \int_{\Omega} \left[ \frac{\omega(d)}{L} + L \, | \nabla d |^2 \right] d\Omega

where :math:`G_c` is the critical fracture energy, :math:`g(d)` the degradation function,
:math:`\omega(d)` the local part of the crack density function, and

.. math::
   c_0 = 4 \int_{0}^{1} \sqrt{\omega(s)} \, ds

is a normalization constant ensuring that the second integral converges to :math:`G_c |\Gamma|`
as :math:`L \rightarrow 0`. The strain energy density is split into a crack-driving part
:math:`\psi^{+}` and a part :math:`\psi^{-}` that is not degraded, which is what prevents cracks
from opening in compression.


Governing equations
--------------------------

Stationarity of the energy with respect to the displacement gives the quasi-static momentum
balance

.. math::
   \nabla \cdot \sigma + b = 0,
   \qquad
   \sigma = g(d) \, \sigma^{+} + \sigma^{-},
   \qquad
   \sigma^{\pm} = \frac{\partial \psi^{\pm}(\varepsilon)}{\partial \varepsilon}

and stationarity with respect to the damage gives the phase-field equation

.. math::
   g^{\prime}(d) \, \mathcal{H}^{+}
   + \frac{G_c}{c_0 L} \left( \omega^{\prime}(d) - 2 L^2 \nabla^2 d \right)
   + \eta \, \dot{d} = 0

The first term is the release of stored strain energy, the second is the generation of fracture
dissipation, and the third is an optional viscous regularization of the damage evolution, with
:math:`\eta` set by ``viscousRegularizationCoeff`` and discretized with a backward Euler scheme.

The crack driving force :math:`\mathcal{H}^{+}` is a history field, which makes the damage
monotonic in time and prevents an unloading path from healing the material:

.. math::
   \mathcal{H}^{+} = \max_{\theta \in [0,t]} \left( \psi^{+}(\varepsilon(\theta)), \, \psi_c \right)

where :math:`\psi_c` is an energy threshold below which no damage develops. Taking the maximum
over the loading history is what makes crack growth irreversible: an unloading path cannot heal
the material. In addition, setting ``irreversibilityFlag="1"`` activates a constraint that
freezes the damage of every node that has reached ``damageUpperBound``, by turning the
corresponding degrees of freedom into Dirichlet conditions, so that fully broken zones stop
evolving.


Energy and stress split
--------------------------

The decomposition of the strain energy density determines how the material responds in
compression, and is selected through the constitutive model:

``DamageElasticIsotropic``
  No split: :math:`\psi^{+} = \psi` and :math:`\psi^{-} = 0`, so the full elastic response is
  degraded. Cracks then carry no stress in compression either.

``DamageSpectralElasticIsotropic``
  Split based on the principal strains :math:`\varepsilon_i` (Miehe et al., 2010):

  .. math::
     \psi^{+} = \frac{\lambda}{2} \langle \text{tr} \, \varepsilon \rangle_{+}^2
              + \mu \sum_{i=1}^{3} \langle \varepsilon_i \rangle_{+}^2

  with :math:`\langle \cdot \rangle_{\pm}` the positive and negative parts. Only the tensile
  part of the spectrum drives the crack and is degraded.

``DamageVolDevElasticIsotropic``
  Split based on the volumetric and deviatoric parts of the strain (Amor et al., 2009):

  .. math::
     \psi^{+} = \frac{K}{2} \langle \text{tr} \, \varepsilon \rangle_{+}^2
              + \mu \, \varepsilon^{dev} : \varepsilon^{dev}

  The deviatoric response is always degraded, whereas the volumetric response is degraded only
  in extension.


Local dissipation function
-----------------------------

Two crack density functions are available, selected by ``localDissipationOption``:

.. list-table::
   :widths: 24 20 20 36
   :header-rows: 1

   * - Option
     - :math:`\omega(d)`
     - :math:`c_0`
     - Comment
   * - ``Linear`` (AT1)
     - :math:`d`
     - :math:`8/3`
     - Has a purely elastic phase up to :math:`\psi_c`, and a damage profile with
       finite support, of half-width :math:`2L`.
   * - ``Quadratic`` (AT2)
     - :math:`d^2`
     - :math:`2`
     - Damage grows as soon as the material is loaded (:math:`\psi_c = 0`), and the
       damage profile has an infinite support.

With the linear model, the crack driving force is floored by the threshold, which is what
produces the elastic phase; with the quadratic model the threshold is not used.


Fracture models
--------------------------

The ``fractureModelType`` attribute of the constitutive model selects the degradation function
and the energy threshold:

``Brittle``
  Quadratic degradation function

  .. math::
     g(d) = (1 - \epsilon) (1 - d)^2 + \epsilon

  where :math:`\epsilon` is ``degradationLowerLimit``, a small residual stiffness that keeps the
  mechanical problem well posed once the crack is formed. The threshold is computed internally:
  :math:`\psi_c = 3 G_c / 16 L` for the linear dissipation, and :math:`\psi_c = 0` for the
  quadratic one. This is the standard brittle model, for which the strength of the material is
  not an input but a consequence of :math:`G_c` and :math:`L`.

``Cohesive``
  Quasi-quadratic (Lorentz-type) degradation function

  .. math::
     g(d) = \frac{(1-d)^2}{(1-d)^2 + m \, d \, (1 + p \, d)},
     \qquad
     m = \frac{3 G_c}{8 L \psi_c},
     \qquad
     p = 1

  with a user-defined threshold :math:`\psi_c` given by ``criticalStrainEnergy``, which can be
  calibrated on the tensile strength of the material. This model gives a cohesive response whose
  peak stress is, to a large extent, insensitive to :math:`L`. It requires the linear
  dissipation and a strictly positive ``criticalStrainEnergy``.

``Nucleation``
  Brittle degradation function, augmented by an external driving force added to the phase-field
  equation, so that the nucleation of cracks follows a prescribed three-dimensional strength
  envelope rather than the strength implied by :math:`G_c` and :math:`L`. The envelope is
  calibrated with ``defaultTensileStrength``, ``defaultCompressiveStrength`` and
  ``defaultDeltaCoefficient``, which must all be provided. The threshold becomes
  :math:`\psi_c = 3 G_c / 16 L + \tfrac{1}{2} c_e`, where :math:`c_e` is the external driving
  force.


.. note::
   The default value of ``fractureModelType`` is not the same for all three constitutive
   models: it is ``Brittle`` for ``DamageElasticIsotropic`` and
   ``DamageVolDevElasticIsotropic``, but ``Cohesive`` for
   ``DamageSpectralElasticIsotropic``. Selecting the spectral split without setting the
   attribute explicitly therefore requests the cohesive model, and the run stops with
   *criticalStrainEnergy must be positive when the Cohesive crack model is used* unless a
   positive threshold is provided. It is good practice to always set ``fractureModelType`` and
   ``localDissipationOption`` explicitly (the latter defaults to ``Linear`` for the three
   models).


Valid combinations
--------------------------

The three categories are not independent. The following combinations are supported, and the
others are rejected at initialization:

.. list-table::
   :widths: 22 22 18 18 20
   :header-rows: 1

   * - Fracture model
     - Local dissipation
     - No split
     - Spectral
     - Volumetric-deviatoric
   * - ``Brittle``
     - ``Linear`` (AT1)
     - yes
     - yes
     - yes
   * - ``Brittle``
     - ``Quadratic`` (AT2)
     - yes
     - yes
     - yes
   * - ``Cohesive``
     - ``Linear`` (AT1)
     - yes
     - yes
     - yes
   * - ``Cohesive``
     - ``Quadratic`` (AT2)
     - no
     - no
     - no
   * - ``Nucleation``
     - ``Linear`` (AT1)
     - yes
     - no
     - no
   * - ``Nucleation``
     - ``Quadratic`` (AT2)
     - no
     - no
     - no


Solution strategy
=========================

The coupled problem is solved with a staggered (sequential) scheme: at each time step, the
momentum balance is solved at fixed damage, then the phase-field equation is solved at fixed
displacement. The number of staggered passes is controlled by the ``NonlinearSolverParameters``
of the coupling solver: ``subcycling="0"`` performs a single pass per step, which is the
operator-split algorithm of Miehe et al. (2010), while enabling subcycling iterates until the
two fields are converged.

Because the damage localizes in a band whose width is set by :math:`L`, the mesh must resolve
that band, and a mesh coarser than :math:`L` will simply not develop a crack. In the benchmarks
distributed with GEOS, the element size in the region swept by the crack is :math:`L/7.5` for
the tension case and :math:`L/3.75` for the shear case of
:ref:`ValidationStudiesPhaseField`. The regularization length is therefore a numerical
parameter as much as a material one, and must always be reported together with the mesh
resolution.

Two attributes help with the robustness of the mechanical solve once cracks are formed.
``degradationLowerLimit`` leaves a small residual stiffness in the broken material, so that the
mechanical problem stays well posed. ``damageUpperBound``, combined with
``irreversibilityFlag="1"``, freezes the damage of the nodes that reach it; its default value of
1.5 is never reached, so it has to be lowered to be active, and the benchmarks use 0.95. A
non-zero ``viscousRegularizationCoeff`` can also be used to stabilize unstable propagation, at
the price of a rate-dependent response, the rate-independent limit being recovered for
:math:`\eta = 0`.

Since the loading is applied quasi-statically, the pseudo-time of the simulation is the loading
parameter. Unstable propagation usually requires a much smaller time step than the initial
loading phase, which is conveniently handled with two ``PeriodicEvent`` blocks with different
``maxEventDt``.


Parameters
=========================

The following attributes are supported in the ``PhaseFieldFracture`` coupling solver:

.. include:: /docs/sphinx/datastructure/PhaseFieldFracture.rst

The damage solver ``PhaseFieldDamageFEM`` supports:

.. include:: /docs/sphinx/datastructure/PhaseFieldDamageFEM.rst

and the constitutive model, here in its spectral-split version, supports:

.. include:: /docs/sphinx/datastructure/DamageSpectralElasticIsotropic.rst


Example
=========================

An example of a valid set of solver blocks is given here, taken from the single-edge notched
benchmark:

.. literalinclude:: ../../../../../inputFiles/phaseField/benchmark/singleEdgeNotch/PhaseFieldFracture_SingleEdgeNotch_base_iterative.xml
  :language: xml
  :start-after: <!-- SPHINX_PHASEFIELD_SOLVER -->
  :end-before: <!-- SPHINX_PHASEFIELD_SOLVER_END -->

.. literalinclude:: ../../../../../inputFiles/phaseField/benchmark/singleEdgeNotch/PhaseFieldFracture_SingleEdgeNotch_base_iterative.xml
  :language: xml
  :start-after: <!-- SPHINX_DAMAGE_SOLVER -->
  :end-before: <!-- SPHINX_DAMAGE_SOLVER_END -->

with the associated constitutive model:

.. literalinclude:: ../../../../../inputFiles/phaseField/benchmark/singleEdgeNotch/PhaseFieldFracture_SingleEdgeNotch_base_iterative.xml
  :language: xml
  :start-after: <!-- SPHINX_MATERIAL -->
  :end-before: <!-- SPHINX_MATERIAL_END -->


Validation and verification
============================

The following examples exercise the model against reference solutions:

* :ref:`ExampleSingleEdgeNotchTension`, mode-I propagation in a single-edge notched block,
* :ref:`ExampleSingleEdgeNotchShear`, mixed-mode propagation on the same specimen,
* :ref:`ExampleThreePointsBendingWithHoles`, crack path driven by the interaction with
  geometric features.


References
=========================

* C. Miehe, M. Hofacker and F. Welschinger, "A phase field model for rate-independent crack
  propagation: Robust algorithmic implementation based on operator splits",
  *Computer Methods in Applied Mechanics and Engineering*, 199:2765-2778, 2010.
* H. Amor, J.-J. Marigo and C. Maurini, "Regularized formulation of the variational brittle
  fracture with unilateral contact: Numerical experiments",
  *Journal of the Mechanics and Physics of Solids*, 57:1209-1229, 2009.
* B. Bourdin, G. A. Francfort and J.-J. Marigo, "The variational approach to fracture",
  *Journal of Elasticity*, 91:5-148, 2008.
