.. _CO2-EOS:

##################################################################################
CO2-brine model 
##################################################################################


Summary
=======

The CO2-brine model implemented in GEOSX includes two components (CO2 and H2O) that are transported by one or two fluid phases (the brine phase and the CO2 phase).
We refer to the brine phase with the subscript :math:`\ell` and to the CO2 phase with the subscript :math:`g` (although the CO2 phase can be in supercritical, liquid, or gas state).
The water component is only present in the brine phase, while the CO2 component can be present in the CO2 phase as well as in the brine phase.
Thus, considering the phase component fractions, :math:`y_{c,p}` (i.e., the fraction of the mass of phase :math:`p` represented by component :math:`c`) the following partition matrix determines the component distribution within the two phases:

.. math::
    \begin{bmatrix}
    y_{CO2,g} & y_{CO2,\ell} \\
         0 & 1            \\
    \end{bmatrix}

The update of the fluid properties is done in two steps:

1) The phase fractions (:math:`\nu_p`) and phase component fractions (:math:`y_{c,p}`) are computed as a function of pressure (:math:`p`), temperature (:math:`T`), component fractions (:math:`z_c`), and a constant salinity.

2) The phase densities (:math:`\rho_p`) and phase viscosities (:math:`\mu_p`) are computed as a function of pressure, temperature, the updated phase component fractions, and a constant salinity.

Once the phase fractions, phase component fractions, phase densities, phase viscosities--and their derivatives with respect to pressure, temperature, and component fractions--have been computed, the :ref:`CompositionalMultiphaseFlow` proceeds to the assembly of the accumulation and flux terms.
Note that the current implementation of the flow solver is isothermal and that the derivatives with respect to temperature are therefore discarded.

The models that are used in steps 1) and 2) are reviewed in more details below.

Step 1: Computation of the phase fractions and phase component fractions (flash)
================================================================================

At initialization, GEOSX performs a preprocessing step to construct a two-dimensional table storing the values of CO2 solubility in brine as a function of pressure, temperature, and a constant salinity.
The user can parameterize the construction of the table by specifying the salinity and by defining the pressure (:math:`p`) and temperature (:math:`T`) axis of the table in the form:

+------------+---------------+-----------------+-----------------+------------------+-----------------+-----------------+------------------+----------+
| FlashModel | CO2Solubility | :math:`p_{min}` | :math:`p_{max}` | :math:`\Delta p` | :math:`T_{min}` | :math:`T_{max}` | :math:`\Delta T` | Salinity | 
+------------+---------------+-----------------+-----------------+------------------+-----------------+-----------------+------------------+----------+

Note that the pressures are in Pascal, and the temperatures are in degree Celsius.
The temperature must be between 10 and 350 degrees Celsius.
The table is populated using the model of Duan and Sun (2003).
Specifically, we solve the following nonlinear CO2 equation of state (equation (A1) in Duan and Sun, 2003) for each pair :math:`(p,T)` to obtain the reduced volume, :math:`V_r`.

.. math::
   \frac{p_r V_r}{T_r} &= 1 + \frac{a_1 + a_2/T^2_r + a_3/T^3_r}{V_r} 
   + \frac{a_4 + a_5/T^2_r + a_6/T^3_r}{V^2_r} + \frac{a_7 + a_8/T^2_r + a_9/T^3_r}{V^4_r} \\
   &+ \frac{a_{10} + a_{11}/T^2_r + a_{12}/T^3_r}{V^5_r} 
   + \frac{a_{13}}{T^3_r V^2_r} \big( a_{14} + \frac{a_{15}}{V^2_r} \big) \exp( - \frac{a_{15}}{V^2_r} )

where :math:`p_r = p / p_{crit}` and :math:`T_r = T / T_{crit}` are respectively the reduced pressure and the reduced temperature.
We refer the reader to Table (A1) in Duan and Sun (2003) for the definition of the coefficients :math:`a_i` involved in the previous equation. 
Using the reduced volume, :math:`V_r`, we compute the fugacity coefficient of CO2, :math:`\ln_{\phi}(p,T)`, using equation (A6) of Duan and Sun (2003).
To conclude this preprocessing step, we use the fugacity coefficient of CO2 to compute and store the solubility of CO2 in brine, :math:`s_{CO2}`, using equation (6) of Duan and Sun (2003):

.. math::
   \ln \frac{ x_{CO2} P }{ s_{CO2} } = \frac{\Phi_{CO2}}{RT} - \ln_{\phi}(p,T) + \sum_c 2 \lambda_c m + \sum_a 2 \lambda_a m + \sum_{a,c} \zeta_{a,c} m^2

where :math:`\Phi_{CO2}` is the chemical potential of the CO2 component, :math:`R` is the gas constant, and :math:`m` is the salinity.
The mole fraction of CO2 in the vapor phase, :math:`x_{CO2}`, is computed with equation (4) of Duan and Sun (2003).
Note that the first, third, fourth, and fifth terms in the equation written above are approximated using equation (7) of Duan and Sun (2003) as recommended by the authors.

During the simulation, Step 1 starts with a look-up in the precomputed table to get the CO2 solubility, :math:`s_{CO2}`, as a function of pressure and temperature.
Then, we compute the phase fractions as:

.. math::
   \nu_{\ell} &= \frac{1 + s_{CO2}}{1 + z_{CO2} / ( 1 - z_{CO2} ) } \\
   \nu_{g} &= 1 - \nu_{\ell}

We conclude Step 1 by computing the phase component fractions as:

.. math::
   y_{CO2,\ell} &= \frac{ s_{CO2} }{ 1 + s_{CO2} } \\
   y_{H2O,\ell} &= 1 - y_{CO2,\ell} \\
   y_{CO2,g} &= 1 \\
   y_{H2O,g} &= 0 
    
   
Step 2: Computation of the phase densities and phase viscosities
================================================================

CO2 phase density and viscosity
-------------------------------

In GEOSX, the computation of the CO2 phase density and viscosity  is entirely based on look-up in precomputed tables.
The user defines the pressure (in Pascal) and temperature (in degrees Celsius) axis of the density table in the form:

+------------+----------------------+-----------------+-----------------+------------------+-----------------+-----------------+------------------+
| DensityFun | SpanWagnerCO2Density | :math:`p_{min}` | :math:`p_{max}` | :math:`\Delta p` | :math:`T_{min}` | :math:`T_{max}` | :math:`\Delta T` |
+------------+----------------------+-----------------+-----------------+------------------+-----------------+-----------------+------------------+

This correlation is valid for pressures less than :math:`8 \times 10^8` Pascal and temperatures less than 800 degrees Celsius.  
Using these parameters, GEOSX internally constructs a two-dimensional table storing the values of density as a function of pressure and temperature.
This table is populated as explained in the work of Span and Wagner (1996) by solving the following nonlinear Helmholtz energy equation for each pair :math:`(p,T)` to obtain the value of density, :math:`\rho_{g}`:

.. math::
   \frac{p}{RT\rho_{g}} = 1 + \delta \phi^r_{\delta}( \delta, \tau )

where :math:`R` is the gas constant, :math:`\delta := \rho_{g} / \rho_{crit}` is the reduced CO2 phase density, and :math:`\tau := T_{crit} / T` is the inverse of the reduced temperature.
The definition of the residual part of the energy equation, denoted by :math:`\phi^r_{\delta}`, can be found in equation (6.5), page 1544 of Span and Wagner (1996).
The coefficients involved in the computation of :math:`\phi^r_{\delta}` are listed in Table (31), page 1544 of Span and Wagner (1996).   
These calculations are done in a preprocessing step.

The pressure and temperature axis of the viscosity table can be parameterized in a similar fashion using the format:

+--------------+----------------------+-----------------+-----------------+------------------+-----------------+-----------------+------------------+
| ViscosityFun | FenghourCO2Viscosity | :math:`p_{min}` | :math:`p_{max}` | :math:`\Delta p` | :math:`T_{min}` | :math:`T_{max}` | :math:`\Delta T` |
+--------------+----------------------+-----------------+-----------------+------------------+-----------------+-----------------+------------------+

This correlation is valid for pressures less than :math:`3 \times 10^8` Pascal and temperatures less than 1220 degrees Celsius.  
This table is populated as explained in the work of Fenghour and Wakeham (1998) by computing the CO2 phase viscosity, :math:`\mu_g`, as follows:

.. math::
   \mu_{g} = \mu_{0}(T) + \mu_{excess}( \rho_{g}, T ) + \mu_{crit}( \rho_{g}, T )  
   
The "zero-density limit" viscosity, :math:`\mu_{0}(T)`, is computed as a function of temperature using equations (3), (4), and (5), as well as Table (1) of Fenghour and Wakeham (1998).
The excess viscosity, :math:`\mu_{excess}( \rho_{g}, T )`, is computed as a function of temperature and CO2 phase density (computed as explained above) using equation (8) and Table (3) of Fenghour and Wakeham (1998).
We currently neglect the critical viscosity, :math:`\mu_{crit}`.
These calculations are done in a preprocessing step.

During the simulation, the update of CO2 phase density and viscosity is simply done with a look-up in the precomputed tables. 

Brine density and viscosity 
---------------------------

The computation of the brine density involves a tabulated correlation presented in Phillips et al. (1981). 
The user specifies the (constant) salinity and defines the pressure and temperature axis of the brine density table in the form:

+------------+----------------------+-----------------+-----------------+------------------+-----------------+-----------------+------------------+----------+
| DensityFun | BrineCO2Density      | :math:`p_{min}` | :math:`p_{max}` | :math:`\Delta p` | :math:`T_{min}` | :math:`T_{max}` | :math:`\Delta T` | Salinity | 
+------------+----------------------+-----------------+-----------------+------------------+-----------------+-----------------+------------------+----------+

The pressure must be in Pascal and must be less than :math:`5 \times 10^7` Pascal.
The temperature must be in degree Celsius and must be between 10 and 350 degrees Celsius.
Using these parameters, GEOSX performs a preprocessing step to construct a two-dimensional table storing the brine density, :math:`\rho_{\ell,table}` for the specified salinity as a function of pressure and temperature using the expression:

.. math::
 
   \rho_{\ell,table} &= A + B x + C x^2 + D x^3 \\
   x &= c_1 \exp( a_1 m ) + c_2 \exp( a_2 T ) + c_3 \exp( a_3 P )

We refer the reader to Phillips et al. (1981), equations (4) and (5), pages 14 and 15 for the definition of the coefficients involved in the previous equation.
This concludes the preprocessing step.

Then, during the simulation, the brine density update proceeds in two steps.
First, a table look-up is performed to retrieve the value of density, :math:`\rho_{\ell,table}`.
Then, in a second step, the density is modified using the method of Garcia (2001) to account for the presence of CO2 dissolved in brine as follows:

.. math::

   \rho_{\ell} = \rho_{\ell,table} + M_{CO2} c_{CO2} - c_{CO2} \rho_{\ell,table} V_{\phi}

where :math:`M_{CO2}` is the molecular weight of CO2, :math:`c_{CO2}` is the concentration of CO2 in brine, and :math:`V_{\phi}` is the apparent molar volume of dissolved CO2.
The CO2 concentration in brine is obtained as:

.. math::

   c_{CO2} = \frac{y_{CO2,\ell} \rho_{\ell,table}}{M_{H2O}(1-y_{CO2,\ell})} 

where :math:`M_{H2O}` is the molecular weight of water. 
The apparent molar volume of dissolved CO2 is computed as a function of temperature using the expression:

.. math::

   V_{\phi} = 37.51 - 9.585 \times 10^{-2} T + 8.740 \times 10^{-4} T^2 - 5.044 \times 10^{-7} T^3

The brine viscosity is controlled by a salinity parameter provided by the user in the form:

+--------------+----------------+----------+
| ViscosityFun | BrineViscosity | Salinity |
+--------------+----------------+----------+

During the simulation, the brine viscosity is updated as a function of temperature using the analytical relationship of Phillips et al. (1981):

.. math::
   \mu_{\ell} = a T + b

where the coefficients :math:`a` and :math:`b` are defined as:

.. math::
   a &= 0.00089 \times 0.000629 (1.0 - \exp( -0.7 m ) ) \\
   b &= 0.00089 (1.0 + 0.0816 m + 0.0122 m^2 + 0.000128 m^3) 
   
where :math:`m` is the user-defined salinity.
   
Parameters
=========================

The model is represented by ``<CO2BrineFluid>`` node in the input.

The following attributes are supported:

.. include:: ../../../coreComponents/schema/docs/CO2BrineFluid.rst

Supported phase names are:

======== ===========
Value     Comment
======== ===========
gas      CO2 phase
water    Water phase
======== ===========

Supported component names are:

============= ===============
Value         Component
============= ===============
co2,CO2       CO2 component
water,liquid  Water component
============= ===============

Example
=======

.. code-block:: xml

    <Constitutive>
        <CO2BrineFluid
          name="fluid"
          phaseNames="{ gas, water }"
          componentNames="{ co2, water }"
          componentMolarWeight="{ 44e-3, 18e-3 }"
          phasePVTParaFiles="{ pvtgas.txt, pvtliquid.txt }"
          flashModelParaFile="co2flash.txt"/>
    </Constitutive>

In the XML code listed above, "co2flash.txt" parameterizes the CO2 solubility table constructed in Step 1.
The file "pvtgas.txt" parameterizes the CO2 phase density and viscosity tables constructed in Step 2, while
the file "pvtliquid.txt" parameterizes the brine density and viscosity tables.
    
References
==========

- Z. Duan and R. Sun, `An improved model calculating CO2 solubility in pure
  water and aqueous NaCl solutions from 273 to 533 K and from 0 to 2000 bar.
  <https://doi.org/10.1016/S0009-2541(02)00263-2>`__, Chemical Geology,
  vol. 193.3-4, pp. 257-271, 2003.

- R. Span and W. Wagner, `A new equation of state for carbon dioxide covering
  the fluid region from the triple-point temperature to 1100 K at pressure up
  to 800 MPa <https://aip.scitation.org/doi/abs/10.1063/1.555991>`__, J. Phys.
  Chem. Ref. Data, vol. 25, pp. 1509-1596, 1996.

- A. Fenghour and W. A. Wakeham, `The viscosity of carbon dioxide
  <https://aip.scitation.org/doi/abs/10.1063/1.556013>`__, J. Phys. Chem. Ref.
  Data, vol. 27, pp. 31-44, 1998.

- S. L. Phillips et al., `A technical databook for geothermal energy
  utilization <https://escholarship.org/content/qt5wg167jq/qt5wg167jq.pdf>`__,
  Lawrence Berkeley Laboratory report, 1981.

- J. E. Garcia, Density of aqueous solutions of CO2. No. LBNL-49023.
  Lawrence Berkeley National Laboratory, Berkeley, CA, 2001.
