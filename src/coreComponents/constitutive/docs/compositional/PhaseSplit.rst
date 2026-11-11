Phase split
===========

The primary purpose of the phase split (or flash) calculation is to determine the number of fluid phases present at equilibrium, the fraction of the total mixture occupying each phase, and the detailed molar composition of these phases. This is a fundamental requirement for compositional modelling, as the physical properties (density, viscosity) of a fluid depend heavily on its phase state.

Phase split equations
---------------------

The phase split relies on two fundamental physical principles: the conservation of mass (material balance) and thermodynamic equilibrium.

**Material Balance**

For a two-phase system (liquid and vapour), the overall mole fraction of a component, :math:`z_i`, must equal the sum of its molar contributions in the liquid phase, :math:`x_i`, and the vapour phase, :math:`y_i`. Introducing the vapour phase fraction, :math:`V`, and the equilibrium ratio (K-value), :math:`K_i = y_i / x_i`, the material balance is expressed by the Rachford-Rice equation (Rachford and Rice, 1952):

.. math::
    \sum_{i=1}^{N_c} \frac{z_i (K_i - 1)}{1 + V (K_i - 1)} = 0

The root of this equation provides the vapour fraction :math:`V`. The phase compositions are subsequently derived via:

.. math::
    x_i = \frac{z_i}{1 + V(K_i - 1)}
    
.. math::
    y_i = K_i x_i

**Thermodynamic Equilibrium**

At thermodynamic equilibrium, the chemical potential, or equivalently the fugacity, of each component must be identical across all existing phases:

.. math::
    f_i^L(P, T, x) = f_i^V(P, T, y)

Using the fugacity coefficient, :math:`\phi_i`, where :math:`f_i = x_i \phi_i P`, the equilibrium condition dictates the K-values:

.. math::
    K_i = \frac{\phi_i^L}{\phi_i^V}

The phase split models in this framework either rigorously solve for these fugacity coefficients using an Equation of State or use simplified, pre-tabulated K-values.

