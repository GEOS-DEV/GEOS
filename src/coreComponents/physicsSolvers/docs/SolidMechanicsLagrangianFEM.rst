.. _SolidMechanicsLagrangianFEM:

#####################################
Lagrangian solid mechanics FEM solver
#####################################

Definition of Terms
===================

.. math::
   i,j,k &\equiv \text {indices over spatial dimensions} \notag \\
   a,b,c &\equiv \text {indices over nodes} \notag \\
   l,m &\equiv \text {indices over volumetric elements} \notag \\
   q,r,s   &\equiv \text {indices over faces} \notag \\
   n &\equiv \text {indices over time} \notag \\
   \Omega &\equiv \text {Volume of continuum body} \notag \\
   \Omega_{crack} &\equiv \text {Volume of open crack} \notag \\
   \Gamma &\equiv \text {External surface of } \Omega \notag \\
   \Gamma_t &\equiv \text {External surface where tractions are applied} \notag \\
   \Gamma_u &\equiv \text {External surface where kinematics are specified} \notag \\
   \Gamma_{crack} &\equiv \text {entire surface of crack} \notag \\
   \Gamma_{cohesive} &\equiv \text {surface of crack subject to cohesive tractions} \notag \\
   \eta_0 &\equiv \text {set of all nodes} \notag \\
   \eta_f & \equiv \text {set of all nodes on flow mesh}  \notag \\
   m &\equiv \text{mass} \notag \\
   \kappa_k &\equiv \text{all elements connected to element k} \notag \\
   \phi & \equiv \text { porosity} \notag \\
   p_f & \equiv \text { fluid pressure} \notag \\
   \mathbf{u} & \equiv \text { displacement} \notag \\
   \mathbf{q} & \equiv \text { volumetric flow rate} \notag \\
   \mathbf{T} & \equiv \text { Cauchy stress} \notag \\
   \rho & \equiv \text { density in the current configuration} \notag \\
   \mathbf{x}& \equiv \text { current position} \notag \\
   \mathbf{w}&\equiv \text { aperture, or gap vector} \notag 

Overview
=========================
The `SolidMechanics_LagrangianFEM` solves the equations of motion as given by:

.. math::
   T_{ij,j} + \rho(b_{i}-\ddot{x}_{i}) = 0. 

The equations of motion are discritized using the Finite Element Method, 
which leads to a discrete set of residual equations:

.. math::   
   (R_{solid})_{ai}=\int\limits_{\Gamma_t} \Phi_a t_i   dA + \int\limits_{\Gamma_c} \left( c_i -p_f n_i  \right)  dA  - \int\limits_\Omega \Phi_{a,j} T_{ij}   dV +\int\limits_\Omega \Phi_a \rho(b_{i}-\Phi_b\ddot{x}_{ib})  dV = 0


Usage
=========================

The following attributes are supported in the input block for `SolidMechanics_LagrangianFEM`:

.. include:: /coreComponents/fileIO/schema/docs/SolidMechanics_LagrangianFEM.rst

Input example
=========================

.. code-block:: xml

  <Solvers>
    <SolidMechanics_LagrangianFEM name="lagsolve"
                                  timeIntegrationOption="QuasiStatic"
                                  discretization="FE1"
                                  targetRegions="Region2">
      <SystemSolverParameters name="solverParams0"
                              useMLPrecond="1"
                              scalingOption="0"
                              krylovTol="1.0e-8"
                              newtonTol="1.0e-4"
                              maxIterNewton="8"
                              verbosityFlag="0"/>
     </SolidMechanics_LagrangianFEM>
  </Solvers>