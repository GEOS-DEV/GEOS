Miscellaneous
=============

Unit and Variable Conversions
-----------------------------

Thermodynamic models, including Equations of State and Flash calculations, are fundamentally derived from molar balances and molar properties. Consequently, the internal computations within the compositional fluid framework are strictly executed using molar variables (e.g., molar fractions, molar densities).

Mass Formulation Support
-----------------------------

While the core thermodynamics rely on molar formulations, many reservoir simulators frame their governing conservation equations in terms of mass to ensure strict mass conservation across numerical schemes. 

To bridge this gap seamlessly, the fluid model supports an automated conversion layer:

1. Entry Conversion: Upon entering the fluid property update routine, if the overall fluid state is provided in mass fractions, the framework dynamically converts these inputs into molar fractions using the designated component molecular weights.
2. Internal Computation: All stability, flash, density, and viscosity calculations proceed purely in molar space.
3. Exit Transformation: Before the update routine concludes, the resulting phase compositions, phase fractions, and molar densities are transformed back into mass variables. Additionally, rigorous chain-rule transformations are applied to all derivatives (with respect to pressure, temperature, and composition) to ensure the Jacobian matrices populated by the simulator remain mathematically consistent with the mass-based primary variables.

Note on Component Limits
--------------------------

When compositional fluid models are evaluated in standalone mode (for example, using the ``PVTDriver`` task (:ref:`PVTDriver`)), the maximum number of allowed components is 9. However, when these fluid models are actively coupled with a flow solver during a full simulation, the maximum number of components is restricted to 5.

