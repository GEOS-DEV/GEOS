.. _CeramicDamageModel:

############################################
Ceramic Damage Model
############################################

Overview, Summary
========================
This damage model is intended for use with damage-field partitioning (DFG) within the
MPM solver, but can also be used without DFG by any solver. It is only appropriate for
schemes implementing explicit time integration. The model is really a hybrid plasticity/
damage model in the sense that we assume damaged material behaves like granular material
and hence follows a modified Mohr-Coulomb law. The modifications are that at low pressures,
the shape of the yield surface is modified to resemble a maximum principal stress criterion,
and at high pressures the shape converges on the von Mises yield surface. The damage
parameter results in softening of the deviatoric response i.e. causes the yield surface to
contract. Furthermore, damage is used to scale back tensile pressure: :math:`p = (1 - d) * p_{Trial}`.
:math:`p_{Trial}` is calculatd as :math:`p_{Trial} = -k * log(J)`, where the Jacobian J of the material motion is
integrated and tracked by this model.




Theory
========================
The following description is sourced from a manuscript in prep (:cite:t:`HomelEtAlInPrep`). 

Strength
------------------------  

The isotropic damage model is designed to capture key features of brittle
material strength while requiring few input parameters. The pressure-dependent
yield surface is defined by the tensile strength :math:`Y_t`, the unconfined
compressive strength :math:`Y_c`, and a high-pressure shear strength
:math:`Y_{\text{max}}`. The tensile and compressive strengths are perturbed by
Weibull scaling of the input parameters :math:`\overline{Y_t}` and
:math:`\overline{Y_c}`, respectively. The parameter :math:`Y_{\text{max}}`
represents the strength in the ductile limit and is not scaled by the Weibull
function.

.. math::

   Y_t = \overline{Y_t}
   \left(
     \frac{\overline{V}}{V}
     \frac{\log R}{\log(1/2)}
   \right)^{1/m},
   \qquad R \in (0,1]

.. math::

   Y_c = \overline{Y_c}
   \left(
     \frac{\overline{V}}{V}
     \frac{\log R}{\log(1/2)}
   \right)^{1/m},
   \qquad R \in (0,1]


If the scaled value of :math:`Y_c` exceeds :math:`Y_{\text{max}}`, which may occur
for small values of the random variable :math:`R` or when the element volume
:math:`V` is much smaller than the reference volume :math:`\overline{V}`, the
compressive strength is capped and rescaled. In this case, we set

.. math::
   Y_c \leftarrow 0.999\,Y_{\text{max}},
   \qquad
   Y_t \leftarrow \frac{\overline{Y_t}}{\overline{Y_c}}\,Y_{\text{max}}.

Perturbations of the yield surface resulting from these control points are shown
in :numref:`fig:damagemodel`b.

The isotropic damage model employs a **hyperelastic volumetric response** and a
**hyperelastic deviatoric response**, assuming linear elasticity and a
*radial return* plasticity algorithm :cite:`Zhou2006ExpandingRing`, with a yield
function that evolves according to a scalar damage variable. A simple piecewise
pressure-dependent strength function is defined, as illustrated in Fig. 1. 


.. _fig:damageEnvelop:

.. figure:: ceramic_damage_piecewise.png
   :align: center
   :width: 75%

Figure 1: Yield surface for the isotropic damage model. Pressure-dependent strength without third-invariant dependence
(:math:`\Gamma` = 1) defined by input values :math:`Y_{c}` and :math:`Y_{t}` , evolving with scalar damage 0 ≤ D ≤ 1 



.. _fig:strengthScaling:

.. figure:: ceramic_damage_strength_scaling.png
   :align: center
   :width: 75%

Figure 2: Weibull scaling of strength through variation of Yt and Yc, subject to the limitation :math:`Y_{c}` < :math:`Y_{max}`.
pressure.


Damage-Dependent Strength
------------------------

As damage accumulates, the material progressively
loses cohesive strength and transitions to a frictional response characterized
by a slope :math:`\mu` at low pressure. This slope also defines the
brittle–ductile transition pressure,

.. math::
   p_2 = \frac{Y_{\text{max}}}{\mu}.


.. math::

   Y = \frac{1}{\Gamma} Y_0
   =
   \frac{1}{\Gamma}
   \begin{cases}
     Y_{1,0}, & \text{if } p_0 < p \le p_1, \\
     Y_{2,0}, & \text{if } p_1 < p \le p_2, \\
     Y_{3,0}, & \text{if } p_2 < p.
   \end{cases}


where :math:`Y_1` linearly interpolates between the unconfined compressive and
tensile strength points, :math:`Y_2` is a blending function, and :math:`Y_3` is
the high-pressure, constant shear strength limit. The pressure-dependent shear
strength is necessary to account for the effects of friction on microcrack
surfaces in brittle materials; however, it is not sufficient to describe tensile
failure mechanisms at low confining stress, where a maximum-principal-stress
failure criterion is more appropriate.

Following :cite:`brannon2009kayenta`, to account for both low-pressure tensile
failure and high-pressure shear failure, a *third-invariant* dependence is
introduced that reduces the strength in triaxial compression (TXC) relative to
triaxial extension (TXE). The resulting yield strength is written as

.. math::
   Y = \frac{1}{\Gamma(\theta)}\,Y_0,

where :math:`\Gamma` is a function of the Lode angle :math:`\theta`.


.. math::

   \theta = \frac{1}{3}\arcsin\!\left[
     \frac{-J_3}{2}\left(\frac{3}{J_2}\right)^{3/2}
   \right].

where the invariants of the stress deviator are

.. math::
   J_2 = \frac{1}{2}\,\mathrm{Tr}\!\left[\boldsymbol{S}^2\right],
   \qquad
   J_3 = \frac{1}{3}\,\mathrm{Tr}\!\left[\boldsymbol{S}^3\right],

with

.. math::
   \boldsymbol{S}
   = \boldsymbol{\sigma}
   - \frac{1}{3}\,\mathrm{Tr}\!\left[\boldsymbol{\sigma}\right]\boldsymbol{I}.

.. math::

   \Gamma =
   \frac{
     4\beta^2(1-\psi^2) + (2\psi - 1)^2
   }{
     (2\psi - 1)\sqrt{4\beta^2(1-\psi^2) + 5\psi^2 - 4\psi}
     + 2\beta(1-\psi^2)
   }.

where :math:`\beta = \cos\!\left(\theta + \frac{\pi}{6}\right)`. The TXC/TXE
strength ratio is scaled based on the slope of the pressure-dependent yield
surface.

.. math::

   \psi = \frac{1}{1+\frac{1}{3}\,\overline{\partial Y_0/\partial p}}.


For the model to produce a tensile strength consistent with the input parameter
:math:`Y_t`, it is necessary to define the interpolation functions
:math:`Y_1` and :math:`Y_2` in terms of a *scaled tensile strength* that accounts
for the :math:`\psi` ratio applied by the :math:`\Gamma` function in triaxial
extension. In addition, a constraint is imposed to ensure that the slope defined
by :math:`Y_t` and :math:`Y_c` is greater than the slope of the fully damaged
yield surface, :math:`\mu`.

.. math::

   Y_{t,0}
   =
   \min\!\left(
     \frac{1}{2}
     \left(
       -Y_c + \sqrt{Y_c}\sqrt{Y_c + 8Y_t}
     \right),
     \frac{3Y_c - \mu Y_c}{3 + \mu}
   \right),
   \qquad
   Y_{t,0} \in \left[\tfrac{1}{2}Y_t,\,2Y_t\right].

With this definition, the intercept of the linear undamaged yield surface is

.. math::

   p_0
   =
   -\frac{2 Y_c Y_{t,0}}
          {3\left(Y_c - Y_{t,0}\right)}.

The pressure-dependent branch of the undamaged yield surface is then defined as

.. math::

   Y_{1,0}
   =
   m_1\left(
     p - (1-\mathcal{D})\,p_0
   \right).


where the slope :math:`m_1` is defined as

.. math::

   m_1
   =
   (1-\mathcal{D})
   \frac{Y_c - Y_{t,0}}{Y_c/3 + Y_{t,0}/3}
   + \mathcal{D}\,\mu.

To smoothly interpolate between the points :math:`(p_1,y_1)` and
:math:`(p_2,y_2)`, we define the function :math:`h(p)` such that
:math:`h'(p_1)=m_1` and :math:`h'(p_2)=0`:

.. math::

   h\!\left(p; p_1, y_1, m_1, p_2, y_2\right)
   =
   (y_1 - y_2)
   \left(
     \frac{p - p_2}{p_1 - p_2}
   \right)^{\,m_1\frac{p_1 - p_2}{y_1 - y_2}}
   + y_2.

The interpolation branch of the undamaged yield surface is then given by

.. math::

   Y_{2,0}
   =
   h\!\left(p; p_1, y_1, m_1, p_2, y_2\right).

Here, :math:`y_1 = Y_c`, :math:`y_2 = Y_{\text{max}}`,
:math:`p_1 = Y_c/3`, and :math:`p_2 = Y_{\text{max}}/\mu`. This function
has slope :math:`m_1` at :math:`(p_1,y_1)` and slope
:math:`m_2 = 0` at :math:`(p_2,y_2)`.

However, when using :math:`\partial Y_{2,0}/\partial p` to evaluate
:math:`\psi`, the yield surface becomes discontinuous
at :math:`p=p_2` for :math:`\mathcal{D}=1`, where the slope is
non-smooth. To produce a continuous yield surface that transitions smoothly
from a maximum-principal-stress criterion at low pressure
(:math:`\mathcal{D}=0`) to a frictional shear surface at high pressure
(:math:`\mathcal{D}=1`), we therefore define



.. math::

   \overline{\partial Y_0/\partial p}
   =
   \begin{cases}
     m_1,
     & \text{if } p_0 < p \le p_1, \\[6pt]
     m_1\left(1 - 3\chi^2 + 2\chi^3\right),
     & \text{if } p_1 < p \le p_2, \\[6pt]
     0,
     & \text{if } p_2 < p.
   \end{cases}

where

.. math::
   \chi = \frac{p - p_1}{p_2 - p_1}.


This formulation is sufficient for the present solution, which employs a simple
pressure cutoff for the vertex treatment and a radial return algorithm. In
models with a different flow rule, it may be necessary to modify this approach
to ensure convexity of the yield surface for all values of :math:`p` and
:math:`\mathcal{D}`.


.. _fig:damagemodel:

.. figure:: ceramic_damage_plots.png
   :align: center
   :width: 75%

Figure 3: Meridional plane of the isotropic damage model yield surface showing triaxial
compression (TXC, top) and triaxial extension (TXE, bottom) strength.
**Left:** yield surface constructed with :math:`\psi` evaluated using the actual
derivative :math:`\partial Y_0/\partial p`.
**Right:** yield surface constructed with :math:`\psi` evaluated using the
smoothed derivative :math:`\overline{\partial Y_0/\partial p}`.



Implementation
========================

``smallStrainUpdateHelper`` Function Overview
------------------------

Updates the Cauchy stress (Voigt form) for a small-strain isotropic elastic material with a
pressure-dependent brittle yield surface and a scalar damage variable :math:`d \in [0,1]`.
The update uses an elastic predictor followed by a yield check and (if needed) a plastic return.
Damage contracts the yield surface (softening) and scales tensile pressure response through
a Jacobian-based pressure computation.

Specifically, this function handles:

#. Elastic trial stress update (via ElasticIsotropicUpdates)
#. Jacobian tracking and Jacobian-based pressure computation :math:`p=-K\ln(J)`
#. Pressure-dependent shear strength with optional third-invariant (Lode-angle) scaling
#. Yield detection and plastic return to the pressure-dependent yield surface
#. Damage evolution via either (a) time-to-failure or (b) energy regularization
#. Optional crack-tip stress concentration modifier
#. Diagnostic plastic strain update (hypoelastic-style)

.. list-table:: A 10-Step Process
   :widths: 5 95
   :header-rows: 1

   * - Step
     - Description
   * - 1
     - **Elastic predictor (trial stress):** Call ElasticIsotropicUpdates::smallStrainUpdate to overwrite the old stress with an elastic trial stress for the full strain increment.
   * - 2
     - **Update Jacobian:** Update the tracked Jacobian using volumetric strain increment,
       :math:`J \leftarrow J\,\exp(\Delta\varepsilon_{xx}+\Delta\varepsilon_{yy}+\Delta\varepsilon_{zz})`.
   * - 3
     - **Compute effective bulk in tension:** If :math:`J>1` (tension), use :math:`K_\mathrm{eff}=(1-d)K` so that unloading from a damaged tensile vertex smoothly approaches :math:`p\to0` as :math:`J\to1`.
   * - 4
     - **Compute trial pressure:** Compute pressure from the Jacobian,
       :math:`p_\mathrm{trial}=-K_\mathrm{eff}\ln(J)`.
   * - 5
     - **Scale strengths (porosity and strengthScale):** Compute
       :math:`Y_t = s\,Y_t^\mathrm{input}(1-\phi)` and :math:`Y_c = s\,Y_c^\mathrm{input}(1-\phi)`,
       enforce limits relative to :math:`Y_\max` and the implied TXC/TXE ratio.
   * - 6
     - **Set tensile vertex pressure:** Compute intact vertex pressure :math:`p_{\min,0}` from :math:`Y_c` and :math:`Y_{t0}` and scale by damage:
       :math:`p_{\min}=(1-d)p_{\min,0}`. If :math:`p_\mathrm{trial}<p_{\min}`, the state is treated as yielding at the vertex.
   * - 7
     - **Compute deviatoric invariants:** Decompose trial stress using twoInvariant::stressDecomposition to obtain deviatoric direction, von Mises magnitude, and compute :math:`J_2` and :math:`J_3`.
   * - 8
     - **Optional crack-tip stress concentration:** If enabled and distanceToCrackTip>0, compute a crack-tip stress concentration factor (stored in crackTipStressConcentration) using a fracture-process-zone radius based on toughness and nominal intact strength.
   * - 9
     - **Yield check and stress update:**
       If not yielding: recompose stress using the Jacobian-based pressure and keep the deviatoric trial state.
       If yielding: call plasticReturn to reconstruct stress on the yield surface (or at the vertex).
   * - 10
     - **Damage + diagnostics:**
       If energy criterion is enabled: iterate (fixed-point / bisection) on damage so dissipation matches :math:`G_c/\ell`, possibly activating surfaceFlag.
       Else: use time-to-failure update :math:`d \leftarrow \min(d + \Delta t / (\ell/c_\mathrm{crack}),\,1)`, applied below a brittle–ductile transition pressure.
       Update diagnostic plastic strain via computePlasticStrainIncrement.

.. raw:: html

   <div style="margin-top: 6em;"></div>

``getStrength`` (Strength Evalutation) Function Overview
------------------------
This function computes the pressure-dependent shear strength :math:`Y(p,d)` used to decide yielding
and to perform the plastic correction. The strength is piecewise in pressure and may be scaled
by a third-invariant (Lode-angle) factor when ``thirdInvariantDependence==1``.
A crack-tip stress concentration factor can further reduce the effective continuum strength.

.. list-table:: A 6-Step Process
   :widths: 5 95
   :header-rows: 1

   * - Step
     - Description
   * - 1
     - **Define pressure breakpoints:** Compute

       .. math::

          p_1 = Y_c/3,
          \qquad
          p_2 = Y_\max/\mu
   * - 2
     - **Low-pressure branch** (:math:`p \le p_1`) **:**
       Evaluate a modified Mohr–Coulomb-like strength via
       ``ceramicY10(p,d,mu,Yt0,Yc)``.
   * - 3
     - **Compute slope proxy for scaling:** Compute :math:`df/dp` using
       ``ceramicdY10dp`` (low branch) or a smoothed variant from
       ``ceramicdY20dp`` (transition).
   * - 4
     - **Optional third-invariant scaling:** If enabled, compute
       :math:`1/\Gamma` via ``thirdInvariantStrengthScaling(J2,J3,dfdp)``
       and multiply the strength by :math:`1/\Gamma`.
   * - 5
     - **Intermediate pressures** (:math:`p_1 < p < p_2`) **:**
       Blend strength smoothly from the low-pressure value at :math:`p_1`
       to :math:`Y_\max` at :math:`p_2` using an exponential form
       (with slope control).
   * - 6
     - **High-pressure cap** (:math:`p \ge p_2`) **:**
       Return the capped strength

       .. math::

          Y = Y_\max

       (optionally scaled by :math:`1/\Gamma`).

       Finally apply the crack-tip modifier

       .. math::

          Y \leftarrow Y/\text{stressConcentration}.





.. raw:: html

   <div style="margin-top: 6em;"></div>



``plasticReturn`` Function Overview
------------------------
This function reconstructs the end-of-step stress so it lies on the yield surface at the
current pressure and damage level. It enforces a tensile vertex and otherwise scales the
deviatoric magnitude back to the admissible strength while preserving deviatoric direction.

.. list-table:: A 5-Step Process
   :widths: 5 95
   :header-rows: 1

   * - Step
     - Description
   * - 1
     - **Check tensile vertex:** If

       .. math::

          p_\mathrm{trial} \le (1-d)\,p_{\min,0},

       set

       .. math::

          p=(1-d)\,p_{\min,0}

       and return an isotropic stress state (zero deviator).
   * - 2
     - **Compute strength at pressure:** Otherwise compute :math:`Y(p,d)` using
       ``getStrength(d,...)`` (including optional third-invariant and crack-tip scaling).
   * - 3
     - **Scale deviatoric magnitude:** Set

       .. math::

          q_\mathrm{new}=\min(q_\mathrm{trial},\,Y).
   * - 4
     - **Recompose stress:** Call ``twoInvariant::stressRecomposition`` with pressure
       and the deviatoric direction to construct the updated stress.
   * - 5
     - **Compute elastic strain energy:** Compute end-of-step elastic strain energy

       .. math::

          \Psi_e = \tfrac12 p^2/K + q_\mathrm{new}^2/(6G),

       used by the energy failure criterion.


.. raw:: html

   <div style="margin-top: 6em;"></div>




``computePlasticStrainIncrement`` Function Overview
------------------------
This function computes a diagnostic plastic strain increment using a hypoelastic decomposition.
It is intended primarily as a plotting variable: the constitutive update enforces yielding by
scaling stresses, and this routine infers the corresponding plastic part of the strain.

.. list-table:: The 6-Step Process
   :widths: 5 95
   :header-rows: 1

   * - Step
     - Description
   * - 1
     - **Compute stress increment:** Form

       .. math::

          \Delta\sigma = \sigma_\mathrm{trial} - \sigma_\mathrm{new}

       (the stress removed by plastic correction).
   * - 2
     - **Isotropic–deviatoric decomposition:** Decompose :math:`\Delta\sigma`
       into isostatic and deviatoric parts using
       ``twoInvariant::stressDecomposition``.
   * - 3
     - **Compute elastic strain increment:** Estimate
       :math:`\Delta\varepsilon_e` using bulk and shear moduli
       (with Voigt shear factors handled explicitly).
   * - 4
     - **Compute plastic strain increment:**

       .. math::

          \Delta\varepsilon_p
          =
          \Delta\varepsilon - \Delta\varepsilon_e.
   * - 5
     - **Rotation handling:** Unrotate old plastic strain, add increment,
       then re-rotate to the end configuration
       (with Voigt shear scaling corrections).
   * - 6
     - **Store updated plastic strain:** Write the updated value to the
       ``plasticStrain`` state variable.


.. raw:: html

   <div style="margin-top: 6em;"></div>




Parameters, User Inputs
========================

.. include:: /coreComponents/schema/docs/CeramicDamage.rst

.. raw:: html

   <div style="margin-top: 6em;"></div>



References
========================

.. bibliography::
   :style: plain
   :filter: docname in docnames
