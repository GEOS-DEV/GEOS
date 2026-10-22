.. _EquilibriumInitialCondition:

####################################################
Hydrostatic Equilibrium Initial Condition
####################################################

Overview
======================

The user can request an initialization procedure enforcing a hydrostatic equilibrium for flow simulations and for coupled flow and mechanics simulations.
The hydrostatic initialization is done by placing one or more **HydrostaticEquilibrium** tag(s) in the **FieldSpecifications** block of the XML input file.
This initialization procedure is described below in the context of single-phase and compositional multiphase flow.
At the end of this document, we compare the hydrostatic equilibrium method to another initialization method, based on the input of x-y-z tables.

Hydrostatic Equilibration
==============================
The objective of the hydrostatic equilibration is to achieve a fluid configuration in which all phases are static. We set the Darcy velocity of each phase, :math:`\mathbf{u}_j`, to zero:

.. math::
   \mathbf{u}_j = -\lambda_j \mathbf{K} \left ( \nabla p_j - \rho_j \mathbf{g} \right ) = 0

Here, :math:`\lambda_j`, :math:`p_j` and :math:`\rho_j` denote the mobility, pressure and density of phase :math:`j`. :math:`\mathbf{K}` is the second order absolute permeability tensor and the gravity vector is represented by :math:`\mathbf{g}`.
Alternatively, equating the pressure gradient with the gravitational gradient results in a no-flow phases' configuration and assuming the gravity vector is aligned with the z-axis:

.. math::
   \frac{\partial p_j}{\partial z} = \rho_j g

where $g = -9.81 \frac{m}{s^2}$.
A mixed discretization is utilized where the pressure derivative is approximated using backward difference while density is evaluated as an average between the reference elevation, :math:`z_{ref}` and current elevation, :math:`z`.
The average density and the discretized pressure equation are given:

.. math::
   {\rho_j}_{avg} = \frac{{\rho_j}_{ref} + \rho_j \left( p(z), T(z), z_i(z) \right )}{2}

.. math::
   p_j(z) = p_{ref} - {\rho_j}_{avg} g \left (z_{ref} - z \right )

The pressure equation is nonlinear because :math:`{\rho_j}_{avg}(z)` is a function of :math:`p_j(z)`. We use the fixed-point method to iteratively solve the pressure equation.

For single phase cases, pressure is computed at the elevation points and then mapped into cell centers. 
However, multiphase systems require phase contact enforcement and composition correction.

Phase Contact Enforcement
------------------------------

Field data usually provide the depths of phase contacts with a high degree of certainty. Therefore, users want to enforce phase contacts precisely. 
Here, the method by which phase contacts are enforced is discussed. Firstly, we define phase contact as the depth at which the capillary pressure, :math:`p_c` is zero:

.. math::
   p_c(z_{contact}) = p_{nw} - p_{w} = 0

$z_{contact}$, $p_{nw}$ and $p_w$ refer to phase contact depth, non-wetting phase pressure and wetting phase pressure respectively.
A nested fixed-point iterative method is used to enforce phase contacts. Given a datum, where pressure of the "primary phase" is known, a phase contact and elevation points between the datum and the phase contact, the algorithm begins by marching from the datum to the phase contact. 
At the phase contact, pressure of "secondary phase" is equated to "primary phase" pressure and the "secondary phase" pressure is then corrected. The process is repeated until the pressure of the two phases at the phase contact are equal upto a tolerace.

Composition Correction
==============================

Once phase pressures are computed during the hydrostatic equilibration step, capillary pressures are then determined. A capillary pressure model relates capillary pressure to the wetting phase saturation. 
Therefore, phase saturations can be back-calculated according to the capillary pressure model. To achieve this, a Newton-based solver is implemented that inverts capillary pressure to obtain the corresponding phase saturations. 
The capillary pressure inversion solver can handle three-phase systems and supports inversion for any generic capillary pressure model.

On the other hand, hydrostatic equilibration requires phase densities. Given a pressure, temperature and overall compositions, the fluid is flashed to obtain phase densities during each iteration of hydrostatic pressure computation. 
The other outputs of the flash are phase compositions, :math:`x_{ij}` and phase fractions, :math:`\gamma_j`. Phase saturations can then be computed from the phase densities and phase fractions. 
Lets denote these saturations as :math:`S_j^{flash}`` and the saturations from the capillary pressure inversion as :math:`S_j^{p_c}`. To account for capillary effects, :math:`S_j^{flash}` has to be equal to :math:`S_j^{p_c}`. 
This is achieved by recombining phases at the corrected overall compositions, without affecting the hydrostatic or thermodynamic equilibrium---the post-correction flash must yield the same density and phase compositions as before the correction. 
The choice of modifying the overall compositions stems from the fact that they are the least certain among the inputs to the flash. Phase fractions are modified using :math:`S_j^{p_c}` and phase densities from pre-correction flash:

.. math::
   \gamma_j^{corr} = \frac{S_j^{p_c} \rho_j}{\sum_j S_j^{p_c} \rho_j}

The phases are then recombined using the corrected phase fractions and the phase compositions from pre-correction flash:

.. math::
   z_i^{corr} = \sum_j \gamma_j^{corr} x_{ij}

Note that volume change upon mixing is ignored.

Single-phase flow parameters
==============================

For single-phase flow, the **HydrostaticEquilibrium** initialization procedure requires the following user input parameters:

* ``datumElevation``: the elevation (in meters) at which the datum pressure is enforced. The user must ensure that the datum elevation is within the elevation range defined by the input mesh. GEOS issues a warning if this is not the case.

* ``datumPressure``: the pressure value (in Pascal) enforced by GEOS at the datum elevation. 

* ``objectPath``: the path defining the groups on which the hydrostatic equilibrium is computed. We recommend using ``ElementRegions`` to apply the hydrostatic equilibrium to all the cells in the mesh. Alternatively, the format ``ElementRegions/NameOfRegion/NameOfCellBlock`` can be used to select only a cell block on which the hydrostatic equilibrium is computed.

.. note::
   In GEOS, the z-axis is positive going upward, this is why the attributes listed in this page are expressed as a function of elevation, not depth. 

Using these parameters and the pressure-density constitutive relationship, GEOS uses a fixed-point iteration scheme to populate a table of hydrostatic pressures as a function of elevation. The fixed-point iteration scheme uses two optional attributes: ``equilibriumTolerance``, the absolute tolerance to declare that the algorithm has converged, and ``maxNumberOfEquilibrationTolerance``, the maximum number of iterations for a given elevation in the fixed point iteration scheme.

In addition, the elevation spacing of the hydrostatic pressure table is set with the optional ``elevationIncrementInHydrostaticPressureTable`` parameter (in meters), whose default value is 0.6096 meters. 
Then, once the table is fully constructed, the hydrostatic pressure in each cell is obtained by interpolating in the hydrostatic pressure table using the elevation at the center of the cell.

.. note::
   The initialization algorithm assumes that the ``gravityVector`` (defined in the **Solvers** XML tag) is aligned with the z-axis. If this is not the case, GEOS terminates the simulation when the **HydrostaticEquilibrium** tag is detected in the XML file. 

Compositional multiphase flow parameters
==========================================

For compositional multiphase flow, the **HydrostaticEquilibrium** initialization procedure follows the same logic but requires more input parameters.
In addition to the required ``datumElevation``, ``datumPressure``, and ``objectPath`` parameters listed above, the user must specify:

* ``componentNames``: the names of the components present in the fluid model. This field is used to make sure that the components provided to **HydrostaticEquilibrium** are consistent with the components listed in the fluid model of the **Constitutive** block. 

* ``componentFractionVsElevationTableNames``: the names of :math:`n_c` tables (where :math:`n_c` is the number of components) specifying the component fractions as a function of elevation. There must be one table name per component, and the table names must be listed in the same order as the components in ``componentNames``. 

* ``temperatureVsElevationTableName``: the names of the table specifying the temperature (in Kelvin) as a function of elevation.

* ``initialPhaseName``: the name of the phase initially saturating the domain. The other phases are assumed to be at residual saturation at the beginning of the simulation. 

* ``phaseContacts``: the elevation of the phase contacts. There must be :math:`n_p - 1` phase contacts (where :math:`n_p` is the number of phases). The phase contacts must be in descending order.

These parameters are used with the fluid density model (depending for compositional flow on pressure, component fractions, and in some cases, temperature) to populate the hydrostatic pressure table, and later initialize the pressure in each cell.

.. note::
   The current initialization algorithm has an important limitation and does not support initial phase contacts (e.g., water-oil, gas-oil, or water-gas contacts). The implementation assumes only one mobile phase in the initial system, identified by the ``initialPhaseName`` attribute. The other phases are assumed at residual saturation. As a result, the system may not be at equilibrium if there is initially more than one mobile phase in the system (for instance if the domain is saturated with gas at the top, and water at the bottom, for instance).

.. note::
   As in the single-phase flow case, GEOS terminates the simulation if **HydrostaticEquilibrium** tag is present in an XML file defining a ``gravityVector`` not aligned with the z-axis.

The full list of parameters is provided below:

.. include:: /docs/sphinx/datastructure/HydrostaticEquilibrium.rst


Examples
=======================

For single-phase flow, a typical hydrostatic equilibrium input looks like:

.. code-block:: xml
	     
   <FieldSpecifications>
   
      <HydrostaticEquilibrium
        name="equil"
        objectPath="ElementRegions"      
        datumElevation="5"
        datumPressure="1e6"/>
      
   </FieldSpecifications>

For compositional multiphase flow, using for instance the CO2-brine flow model, a typical hydrostatic equilibrium input looks like:

.. code-block:: xml

   <FieldSpecifications>		
	     
      <HydrostaticEquilibrium
        name="equil"
        objectPath="ElementRegions"      
        datumElevation="28.5"
        datumPressure="1.1e7"
        initialPhaseName="water"
        componentNames="{ co2, water }"
        phaseContacts="{ 50 }"
        componentFractionVsElevationTableNames="{ initCO2CompFracTable,
                                                  initWaterCompFracTable }"
        temperatureVsElevationTableName="initTempTable"/>

   </FieldSpecifications>

In this case, a possible way to provide the three required tables is:

.. code-block:: xml

   <Functions>

     <TableFunction
       name="initCO2CompFracTable"
       coordinates="{ 0.0, 10.0, 20.0, 30.0 }"
       values="{ 0.04, 0.045, 0.05, 0.055 }"/>

     <TableFunction
       name="initWaterCompFracTable"
       coordinates="{ 0.0, 10.0, 20.0, 30.0 }"
       values="{ 0.96, 0.955, 0.95, 0.945 }"/>

     <TableFunction
       name="initTempTable"
       coordinates="{ 0.0, 15.0, 30.0 }"
       values="{ 358.15, 339.3, 333.03 }"/>
     
   </Functions>

Note that the spacing of the two component fraction tables must be the same, but the spacing of the temperature table can be different.

Expected behavior and comparison with another initialization method
=====================================================================

As illustrated in :ref:`TutorialFieldCase`, users can also use multiple **FieldSpecification** tags to impose initial fields, such as the pressure, component fractions, and temperature fields.
To help users select the initialization method that best meets their needs, we summarize and compare below the two possible ways to initialize complex, non-uniform initial fields for compositional multiphase simulations in GEOS.

Initialization using **HydrostaticEquilibrium**
-----------------------------------------------

This is the initialization procedure that we have described in the first sections of this page.
In **HydrostaticEquilibrium**, the initial component fractions and temperatures are provided as a function of elevation only, and the hydrostatic pressure is computed internally before the simulation starts.
The typical input was illustrated for a CO2-brine fluid model in the previous paragraph.

Expected behavior:

* If **FieldSpecification** tags specifying initial pressure, component fractions, and/or temperature are included in an XML input file that also contains the **HydrostaticEquilibrium** tag, the **FieldSpecification** tags are ignored by GEOS. In other words, only the pressure, component fractions, and temperature fields defined with the **HydrostaticEquilibrium** tag as a function of elevation are taken into account.

* In the absence of source/sink terms and wells, the initial flow residual should be smaller than :math:`10^-6`. Similarly, in coupled simulations, the residual of the mechanical problem should be close to zero.

Initialization using **FieldSpecification** tags
------------------------------------------------

This is the initialization method illustrated in :ref:`TutorialFieldCase`.
The user can impose initial pressure, component fractions, and temperature fields using **FieldSpecification** tags, such as, for a two-component CO2-brine case:

.. code-block:: xml

   <FieldSpecifications>		
	     
     <FieldSpecification
       name="initialPressure"
       initialCondition="1"
       setNames="{ all }"
       objectPath="ElementRegions"
       fieldName="pressure"
       scale="1"
       functionName="initialPressureTableXYZ"/>

     <FieldSpecification
       name="initialCO2CompFraction"
       initialCondition="1"
       setNames="{ all }"
       objectPath="ElementRegions"
       fieldName="globalCompFraction"
       component="0"
       scale="1"
       functionName="initialCO2CompFracTableXYZ"/>

     <FieldSpecification
       name="initialWaterCompFrac"
       initialCondition="1"
       setNames="{ all }"
       objectPath="ElementRegions"
       fieldName="globalCompFraction"
       component="1"
       scale="1"
       functionName="initialWaterCompFracTableXYZ"/>

     <FieldSpecification
       name="initialTemperature"
       initialCondition="1"
       setNames="{ all }"
       objectPath="ElementRegions"
       fieldName="temperature"
       scale="1"
       functionName="initialTemperatureTableXYZ"/>
       
   </FieldSpecifications>		       

In this input method, ``initialPressureTableXYZ``, ``initialCO2CompFracTableXYZ``, ``initialWaterCompFracTableXYZ``, and ``initialTemperatureTableXYZ`` are tables describing these initial fields as a function of the x, y, and z spatial coordinates.
Then, the cell-wise values are determined by interpolating in these tables using the coordinates of the center of the cells.

Expected behavior:

* In this approach, it is the responsibility of the user to make sure that these initial fields satisfy a hydrostatic equilibrium. If not, the model will equilibrate itself during the first time steps of the simulation, possibly causing large initial changes in pressure, component fractions, and temperature.

* If the initial state imposed by the **FieldSpecification** tags is not at equilibrium, the displacements produced by coupled flow and mechanics simulations should be interpreted with caution, as these displacements are computed with respect to a non-equilibrium initial state.

* This method is suited to impose initial fields in complex cases currently not supported by the **HydrostaticEquilibrium** tag (e.g., in the presence of phase contacts, capillary pressure, etc). Specifically, the user can equilibrate the model using other means (such as using another simulator, or running a few steps of GEOS), retrieve the equilibrated values, convert them into x-y-z tables, and impose them in the new GEOS simulations using **FieldSpecification** tags.  
