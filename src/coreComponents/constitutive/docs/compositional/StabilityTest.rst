.. _StabilityTest:

Stability test
-------------------------------
The first step is to determine if the provided mixture with total molar fractions :math:`z_c` is stable
as a single phase at the current pressure :math:`p` and temperature :math:`T`. However, this can only
be confirmed through stability testing.

The stability of a mixture is traditionally assessed using the Tangent Plane Distance (TPD) criterion
developed by Michelsen (1982a). This criterion states that a phase with composition :math:`z` is stable
at a specified pressure :math:`p` and temperature :math:`T` if and only if 

.. math::
  g(y) = \sum_{i=1}^{N_c}y_i\left(\ln y_i + \ln \phi_i(y) - \ln z_i - \ln \phi_i(z) \right) \geq 0

for any permissible trial composition :math:`y`, where :math:`\phi_i` denotes the fugacity
coefficient of component :math:`i`. 

To determine stability of the mixture this testing in initiated from a basic starting point, based on
Wilson K-values, to get both a lighter and a heavier trial mixture. The two trial mixtures are
calculated as :math:`y_i = z_i/K_i` and :math:`y_i = z_iK_i` where :math:`K_i` are defined by

.. math::
  K_i = \frac{P_{ci}}{p}\exp\left( 5.37( 1 + \omega_i ) \left(1-\frac{T_{ci}}{T}\right)\right)
  
where :math:`P_{ci}` and :math:`T_{ci}` are respectively, the critical pressure and temperature of
component :math:`i` and :math:`\omega_i` is the accentric factor of component :math:`i`.

The stability problem is solved by observing that a necessary condition is that :math:`g(y)` must
be non-negative at all its stationary points. The stationarity criterion can be expressed as

.. math::
  \ln y_i + \ln \phi_i(y) - h_i = k \hspace{1cm} i=1,2,3,\ldots,N_c

where :math:`h_i = \ln z_i + \ln \phi_i(z)` is a constant parameter dependent on the feed composition
:math:`z` and :math:`k` is an undetermined constant. This constant can be further incorporated into
the equation by defining the unnormalized trial phase moles :math:`Y_i` as

.. math::
  Y_i = \exp(-k)y_i

which reduces the stationarity criterion to

.. math::
  \ln Y_i + \ln \phi_i(y) - h_i = 0

with the mole fractions :math:`y_i` related to the trial phase moles :math:`Y_i` by

.. math::
  y_i = Y_i/\sum_jY_j

With the two starting mixtures, the stationarity condition is solved using successive substitution to
determine the stationary points. If both initial states converge to a solution which has :math:`g(y)\geq 0`
then the mixture is deemed to be stable, otherwise it is deemed unstable.

Phase labeling
----------------------------------------
Once it is confirmed that the fluid with composition :math:`z` is stable as a single phase at the current
pressure and temperature, it must be labeled as either 'liquid' or 'vapor'. This is necessary only to apply
the correct relative permeability function for calculating the phase's flow properties. The properties of the
fluid (density, viscosity) are unchanged by the assignment of the label.

Determining the mixture's true critical point is the most rigorous method for labeling. It is however expensive
and may not always be necessary. As such, a simple correlation for pseudo-critical temperature is used and this
is expected to be sufficiently accurate for correct phase labeling, except under some specific conditions.

The Li-correlation is a weighted average of the component critical temperatures and is used to determine the label
applied to the mixture. The Li pseudo-critical temperature is calcaulated as

.. math::
  T_{cp} = \frac{\sum_{i=1}^{N_c}T_{ci}V_{ci}z_{i}}{\sum_{i=1}^{N_c}V_{ci}z_{i}}

where :math:`V_{ci}` and :math:`T_{ci}` are respectively the critical volume and temperature of component
:math:`i`. This is compared to the current temperature :math:`T` such that if :math:`T_{cp}<T` then the mixture
is labeled as vapor and as liquid otherwise.
