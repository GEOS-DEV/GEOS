.. _NegativeTwoPhaseFlash:

Negative two-phase flash
------------------------
When a cell is identified as having an unstable mixture, it is necessary to determine the amounts in the liquid
and vapor phases through phase splitting. This phase split is calculated by ensuring that the two phases are
in thermodynamic equilibrium. For a system to be in thermodynamic equilibrium, the fugacities of each component
in both the liquid and vapor phases must be equal:

.. math::
    \phi_{iL} = \phi_{iV} \hspace{1cm} i=1,2,3,\ldots,N_c

where :math:`\phi_{iL}` is the fugacity of component :math:`i` in the liquid phase and :math:`\phi_{iL}` is the
fugacity of component :math:`i` in the vapor phase.

Fugacities are functions of temperature, pressure, and composition:

.. math::
    \phi_{iL} = \phi_{iL}(T, p, x_i) \hspace{1cm} i=1,2,3,\ldots,N_c

and 

.. math::
    \phi_{iV} = \phi_{iV}(T, p, y_i) \hspace{1cm} i=1,2,3,\ldots,N_c

and are calculated directly from an equation of state.

Equilibrium constants, also known as K-values, are defined for each component as:

.. math::
    K_i = \frac{y_i}{x_i}

where :math:`x_i` is the mole fraction of component :math:`i` in the liquid phase and :math:`y_i` is the
mole fraction of component :math:`i` in the vapor phase. If we denote :math:`V` as the mole fraction
of the vapor phase, the material balance indicates that the mole fractions of each component in the liquid
and vapor phases are given by:

.. math::
    x_i = \frac{z_i}{1 + (K_i - 1)V}

and

.. math::
    y_i = \frac{K_i z_i}{1 + (K_i - 1)V}

The value of :math:`V` corresponding to a given set of K-values is determined by solving the
so called Rachford and-Rice equation:

.. math::
    F(V) = \sum_{i=1}^{N_c} \left(x_i - y_i\right) = \sum_{i=1}^{N_c} \frac{z_i(1 - K_i)}{1 + (K_i - 1)V} = 0

The flash calculation process is as follows:

#. Once the mixture is confirmed to be stable, an initial set of K-values is chosen, typically using Wilson's formula.

#. Given :math:`z_i` and :math:`K_i`, the Rachford-Rice equation is solved to determine the molar fraction of vapor,  :math:`V`. This is initially solved using successive substitution, followed by Newton iterations once the residual is sufficiently reduced.

#. After  :math:`V` is calculated, the corresponding liquid and vapor mole fractions, :math:`x_i` and :math:`y_i`, are computed.

#. These phase compositions are then used to calculate the component fugacities :math:`\phi_{iL}` and :math:`\phi_{iV}` in the liquid and vapor phases using the equation of state.

#. Convergence is reached when the fugacities are equal for all components. The convergence criterion is defined as:

   .. math::
       \sum_{i=1}^{N_c} \left( \phi_{iL} - \phi_{iV} \right)^2 < \varepsilon
   
   where :math:`\varepsilon` is the convergence tolerance.

#. If convergence is not achieved, successive substitution is used to update the set of K-values for the next iteration. The new K-values at iteration  :math:`t+1` are given by:

   .. math::
       K_i^{(t+1)} = K_i^{(t)} \frac{\phi_{iL}}{\phi_{iV}}
