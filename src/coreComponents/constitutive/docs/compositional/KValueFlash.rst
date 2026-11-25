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

* ``kValueTables``: List of function names linking to the pre-tabulated K-values. For an :math:`N`-phase system, this requires :math:`(N-1) \times N_c` function names.
* ``pressureCoordinates``: List of pressure values for interpolation grid generation (optional). Unit: [Pa].
* ``temperatureCoordinates``: List of temperature values for interpolation grid generation (optional). Unit: [K].