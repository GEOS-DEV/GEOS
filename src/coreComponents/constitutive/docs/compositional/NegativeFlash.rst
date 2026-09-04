Negative two-phase flash
------------------------

The negative flash is a robust, rigorous approach to solving the phase equilibrium problem using an Equation of State (EoS). 

Stability analysis
~~~~~~~~~~~~~~~~~~

Before executing a flash calculation, the model must ascertain if the homogeneous single-phase mixture is thermodynamically unstable at the specified pressure :math:`P` and temperature :math:`T`. This is performed via a tangent plane distance (TPD) analysis developed by Michelsen (Michelsen, 1982a; Michelsen, 1982b). The TPD assesses whether the Gibbs free energy of the mixture can be minimised by splitting into multiple phases. A phase with feed composition :math:`z` is stable if and only if:

.. math::
    g(\mathbf{y}) = \sum_{i=1}^{N_c} y_i \left( \ln y_i + \ln \phi_i(\mathbf{y}) - \ln z_i - \ln \phi_i(\mathbf{z}) \right) \ge 0

for any permissible trial composition :math:`\mathbf{y}`, where :math:`\phi_i` denotes the fugacity coefficient of component :math:`i`. 

To determine stability, testing is initiated from basic starting points using Wilson's K-values (Wilson, 1969; Whitson and Torp, 1981) to generate a lighter trial mixture (:math:`y_i = z_i / K_i`) and a heavier trial mixture (:math:`y_i = z_i K_i`). Wilson's K-values are defined as:

.. math::
    K_i = \frac{P_{ci}}{P} \exp \left( 5.37 (1 + \omega_i) \left( 1 - \frac{T_{ci}}{T} \right) \right)

where :math:`P_{ci}`, :math:`T_{ci}`, and :math:`\omega_i` are the critical pressure, critical temperature, and acentric factor of component :math:`i`, respectively.

The stationarity condition, :math:`\ln Y_i + \ln \phi_i(\mathbf{y}) - h_i = 0` (where :math:`h_i = \ln z_i + \ln \phi_i(\mathbf{z})` and :math:`Y_i` are unnormalized trial phase moles), is solved using successive substitution. If both initial states converge to a solution where :math:`g(\mathbf{y}) \ge 0`, the mixture is stable; otherwise, it is unstable.

Phase labeling
~~~~~~~~~~~~~~

If the stability test confirms a single stable phase, the two-phase flash calculation is bypassed. To appropriately label the single phase as liquid or vapour for relative permeability evaluations, the model utilises a "Li-temperature" correlation. This pseudo-critical temperature is weighted by the critical volumes (:math:`V_{ci}`):

.. math::
    T_{cp} = \frac{\sum_{i=1}^{N_c} T_{ci} V_{ci} z_i}{\sum_{i=1}^{N_c} V_{ci} z_i}

If :math:`T_{cp} < T`, the mixture is logically categorised as a vapour; otherwise, it is labelled as a liquid.

Flash calculation
~~~~~~~~~~~~~~~~~

If the mixture is deemed unstable, the system must determine the phase split by ensuring thermodynamic equilibrium, meaning the fugacities in the liquid and vapour phases must be equal (:math:`\phi_{iL} = \phi_{iV}`). 

Based on the material balance and the equilibrium constants (:math:`K_i = y_i/x_i`), the mole fractions in the liquid (:math:`x_i`) and vapour (:math:`y_i`) phases are given by:

.. math::
    x_i = \frac{z_i}{1 + (K_i - 1)V}, \quad y_i = \frac{K_i z_i}{1 + (K_i - 1)V}

The vapour molar fraction, :math:`V`, is determined by solving the Rachford-Rice equation alongside the EoS:

.. math::
    F(V) = \sum_{i=1}^{N_c} (x_i - y_i) = \sum_{i=1}^{N_c} \frac{z_i(1 - K_i)}{1 + (K_i - 1)V} = 0

The flash algorithm:

1. Initialize: An initial set of K-values is chosen using Wilson's formula.
2. Solve for V: The Rachford-Rice equation is solved for the vapour fraction :math:`V` (using successive substitution, followed by Newton iterations).
3. Compute compositions: The corresponding liquid (:math:`x_i`) and vapour (:math:`y_i`) mole fractions are computed.
4. Evaluate EoS: These compositions are used to calculate the component fugacities :math:`\phi_{iL}` and :math:`\phi_{iV}` via the Equation of State.
5. Check convergence: Convergence is reached when :math:`\sum_{i=1}^{N_c} (\ln \phi_{iL} - \ln \phi_{iV})^2 < \epsilon`.
6. Update K-values: If not converged, the algorithm employs successive substitution iterations, constantly updating the K-values using fugacity coefficients derived from the EoS:

.. math::
    K_i^{(new)} = K_i^{(old)} \frac{\phi_{iL}}{\phi_{iV}} = K_i^{(old)} \exp \left( \ln \phi_i^L - \ln \phi_i^V \right)

The term "negative" flash arises because the algorithm temporarily permits the vapour fraction, :math:`V`, to converge to values slightly outside the physically meaningful bounds of :math:`[0, 1]`. This mathematical relaxation prevents the solver from getting artificially trapped at phase boundaries, vastly improving convergence robustness near the critical point or saturation envelopes. Upon convergence, if :math:`V` is negative or greater than unity, the system truncates to a single-phase solution.

Recommendations
~~~~~~~~~~~~~~~

The negative two-phase flash is the recommended standard for rigorous compositional modelling of miscible gas injection, primary depletion of volatile oils, and scenarios where fluid compositions change drastically over time. An extended three-phase variant is available for immiscible water systems.

Parameters
~~~~~~~~~~~~~~~~

* ``stabilityThreshold``: Tangent plane distance below which a mixture is unstable (default: -1.0e-8). Unit: [dimensionless].
* ``stabilityTolerance``: Tolerance for stationarity in the stability test. Unit: [dimensionless].
* ``stabilityMaxIterations``: Maximum successive substitution steps for stability analysis. Unit: [dimensionless].
* ``flashTolerance``: Convergence tolerance for the fugacity ratio error. Unit: [dimensionless].
* ``flashMaxIterations``: Maximum successive substitution steps for the flash solve. Unit: [dimensionless].