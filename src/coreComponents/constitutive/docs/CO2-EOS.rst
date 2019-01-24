##################################################################################
Calculation of thermodynamic properties of CO2 and brine 
##################################################################################


Summary
=================================

The equation-of-state and viscosity of CO2 under both sub- and super-critical conditions are computed as a function of both pressure and temperature, using the empirical equations developed by Span & Wagner (1996) and Fenghour & Wakeman (1998), respectively.

The brine density, which depends on pressure, temperature, and salinity, is calculated based on the correlation developed by Phillips et al.(1981) 


Function Description
=================================

**Function of CO2 density calculation**


double PVTPackage::CO2Model::ComputeMassDensity(double P, double T)

``input parameters:`` 

* P - pressure (Pa),  P < 800 MPa
* T - temperature (K) 200K <= T <= 1100K

``return parameter:`` 

* density (kg/m^3)


**Function of CO2 viscosity calculation**

double PVTPackage::CO2Model::ComputeVisc(double P, double T)

``input parameters:``

* P - pressure (Pa),  P < 300 MPa
* T - temperature (K) 200K <= T <= 1500K 

``return parameter:`` 

* viscosity (Pa.s)

**Function of brine density calculation**

double PVTPackage::BrineModel::ComputeMassDensity(double P, double T, double salinity)

``input parameters:``

* P - pressure (Pa),  P <= 50 MPa
* T - temperature (C), 10C <= T <= 350C 
* salinity - brine salinity, (molal) salinity <= 5 molal 
  
``return parameter:`` 

* density (kg/m^3)

References
=================================

- A. Fenghour & W. A. Wakeman, “The viscosity of carbon dioxide”, J. Phys. Chem. Ref. Data, vol. 27, pp. 31-44, 1998.

- S. L. Phillips et al., "A technical databook for geothermal energy utilization", Lawrence Berkeley Laboratory report, 1981.

- R. Span & W. Wagner., “A new equation of state for carbon dioxide covering the fluid region from the triple-point temperature to 1100 K at pressure up to 800 MPa”, J. Phys. Chem. Ref. Data, vol. 25, pp. 1509-1596, 1996. 









