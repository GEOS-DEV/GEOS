K-value flash
-------------

The K-value flash model offers a computationally simplified alternative to the rigorous equation of state flash calculation. Rather than iteratively evaluating complex equations of state to find fugacity coefficients, this model determines phase splits using pre-tabulated K-values. 

The model supports both two-phase and three-phase variants. The two-phase flash calculates equilibrium partitioning between a liquid and a vapour phase. The three-phase flash introduces an additional aqueous phase, allowing for the predefined partitioning of components into water alongside the hydrocarbon phases.

Methodology
~~~~~~~~~~~

In this formulation, the equilibrium ratio for each component in the secondary phases relative to the primary liquid phase, :math:`K_i = y_i / x_i`, is treated strictly as a predefined function of pressure and temperature.

Users input these functions as either tables or symbolic functions. To optimize performance during the simulation, the model discretizes these functions during initialization to construct a multi-dimensional hypercube of K-values spanning the phase, component, pressure, and temperature dimensions.

The pressure and temperature coordinates for this hypercube grid can be explicitly chosen and specified by the user. If they are omitted, the model automatically inspects the provided table functions, extracts all unique coordinate points from the underlying data, and constructs the grid dynamically. If only a single point is found, an artificial offset is added to permit interpolation.

Equilibrium calculation
~~~~~~~~~~~~~~~~~~~~~~~

During the simulation step, the model performs bilinear interpolation on the pre-built hypercube to evaluate the K-values and their partial derivatives at the current pressure and temperature. 

Because the K-values are explicitly known, the iterative fugacity updates are bypassed. Instead, the vapour fraction, :math:`V`, is found by directly solving the Rachford-Rice equation (Rachford and Rice, 1952) using successive substitution and Newton iterations:

.. math::
    \sum_{i=1}^{N_c} \frac{z_i (K_i - 1)}{1 + V(K_i - 1)} = 0

Once the vapour fraction is determined, the phase compositions are directly evaluated:

.. math::
    x_i = \frac{z_i}{1 + V(K_i - 1)}, \quad y_i = K_i x_i

This direct solution approach intrinsically supports a "negative" flash. The mathematical solver allows the vapour fraction, :math:`V`, to temporarily converge to values outside the physically meaningful bounds of :math:`[0, 1]`. If the resulting :math:`V \le 0` or :math:`V \ge 1`, the system cleanly truncates the compositions to represent a single-phase liquid or single-phase vapour solution, respectively.

Recommendations and pitfalls
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The K-value flash is highly computationally efficient and is recommended for thermal recovery processes (like steam injection) or relatively stable black-oil systems where fluid compositions remain within a narrow, predictable envelope.

Pitfalls: In reality, thermodynamic K-values are heavily dependent on the overall mixture composition, not just pressure and temperature. Relying on tabulated K-values implies assuming a fixed compositional path. Care should be taken when using this model for highly variable compositional processes such as miscible gas injection (e.g., CO2 or hydrocarbon gas flooding), or strong depletion of gas condensates, where composition dependency is critical to accurate phase behaviour.

Model parameters
~~~~~~~~~~~~~~~~

The K-value flash fluid models are assigned using specific XML blocks that pair the flash model with viscosity and density models.

Here is an example demonstrating how to specify the parameters for a two-phase system using the K-value flash paired with the Phillips brine model:

.. code-block:: xml

        <CompositionalTwoPhaseKValueFluidPhillipsBrine
            name="FLUID"
            phaseNames="{ gas, water }"
            equationsOfState="{ PengRobinson, SoreideWhitson }"
            componentNames="{ CH4, CO2, H2S, H2O }"
            componentMolarWeight="{ 1.60425e-02, 4.40095e-02, 3.40809e-02, 1.80153e-02 }"
            componentCriticalPressure="{ 4.59920e+06, 7.37730e+06, 9.00000e+06, 2.20640e+07 }"
            componentCriticalTemperature="{ 1.90564e+02, 3.04128e+02, 3.73100e+02, 6.47096e+02 }"
            componentCriticalVolume="{ 9.86278e-05, 9.41185e-05, 9.81354e-05, 5.59480e-05 }"
            componentAcentricFactor="{ 1.14200e-02, 2.23940e-01, 1.00500e-01, 3.44300e-01 }"
            componentBinaryCoeff="{
                { 0.0000, 0.0000, 0.0000, 0.4850 },
                { 0.0000, 0.0000, 0.0000, 0.1896 },
                { 0.0000, 0.0000, 0.0000, 0.1353 },
                { 0.4850, 0.1896, 0.1353, 0.0000 }
            }"
            kValueTables="{ KV_CH4, KV_CO2, KV_H2S, KV_H2O }"
            waterCompressibility="4.1483E-10"
            salinity="2.3" /> 

Parameters
~~~~~~~~~~~~~~~~

* ``kValueTables``: List of function names linking to the pre-tabulated K-values. For an :math:`N`-phase system, this requires :math:`(N-1) \times N_c` function names.
* ``pressureCoordinates``: List of pressure values for interpolation grid generation (optional).
* ``temperatureCoordinates``: List of temperature values for interpolation grid generation (optional).