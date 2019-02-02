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

Newmark Method
---------------
The Newmark update consists of the following equations:

.. math::
   M a^{n+1} + C v^{n+1} + K u^{n+1} &= F_{n+1},  \\
   u^{n+1} &= u^n + v^{n+1/2},  \\
   u^{n+1} &= u^n + \left( v^{n} + \inv{2} \left[ (1-2\beta) a^n + 2\beta a^{n+1} \right] \Delta t \right) \Delta t, \\
   v^{n+1} &= v^n + \left[(1-\gamma) a^n + \gamma a^{n+1} \right] \Delta t.

As intermediate quantities we can form an estimate for the end of step displacement and midstep velocity by
assuming zero end-of-step acceleration.

.. math::
   \tilde{u}^{n+1} &= u^n + \left( v^{n} + \inv{2}  (1-2\beta) a^n  \Delta t \right) \Delta t = u^n + \hat{\tilde{u}}\\
   \tilde{v}^{n+1} &= v^n + (1-\gamma) a^n  \Delta t =  v^n + \hat{\tilde{v}} 

This gives

.. math::
   u^{n+1} &= \tilde{u}^{n+1} + \beta a^{n+1} \Delta t^2 \\
   v^{n+1} &= \tilde{v}^{n+1} + \gamma a^{n+1} \Delta t 

using \eqns{eqn:d}{eqn:v}, we can get acceleration and velocity in terms of displacement

.. math::
   a^{n+1} &= \frac{1}{\beta \Delta t^2} \left(u^{n+1} - \tilde{u}^{n+1} \right)  = \frac{1}{\beta \Delta t^2} \left( \hat{u} - \hat{\tilde{u}} \right) \\
   v^{n+1} &= \tilde{v}^{n+1} + \frac{\gamma}{\beta \Delta t} \left(u^{n+1} - \tilde{u}^{n+1} \right) = \tilde{v}^{n+1} + \frac{\gamma}{\beta \Delta t} \left(\hat{u} - \hat{\tilde{u}} \right)

plugging these into \eqn{eqn:newmarkEOM}

.. math::
   M \left(\frac{1}{\beta \Delta t^2} \left(\hat{u} - \hat{\tilde{u}} \right)\right) + C \left( \tilde{v}^{n+1} + \frac{\gamma}{\beta \Delta t} \left(\hat{u} - \hat{\tilde{u}} \right) \right) + K u^{n+1} &= F_{n+1}  \\

If we use simple Rayliegh damping

.. math::
   C = a_{mass} M + a_{stiff} K

Repose this in context of a nonlinear residual problem

.. math::
    (R_{ss}^e)_{ai} &= 
        \int\limits_{\Gamma_t^e} \Phi_a t_i   dA  
        + \int\limits_{\Gamma_c} \Phi_a c_i   dA \notag \\
        &- \int\limits_{\Omega^e} \Phi_{a,j} \left(T_{ij}^{n+1}+  a_{stiff} \left(\pderiv{T_{ij}^{n+1}}{\hat{u}_{bk}} \right)_{elastic} \left( \tilde{v}_{bk}^{n+1} + \frac{\gamma}{\beta \Delta t} \left(\hat{u}_{bk} - \hat{\tilde{u}}_{bk} \right) \right) \right)  dV \notag \\
        &+\int\limits_{\Omega^e} \Phi_a \rho \left(b_{i}- \Phi_b  \left( a_{mass} \left( \tilde{v}_{bi}^{n+1} + \frac{\gamma}{\beta \Delta t} \left(\hat{u}_{bi} - \hat{\tilde{u}}_{bi} \right) \right) + \frac{1}{\beta \Delta t^2}  \left( \hat{u}_{bi} - \hat{\tilde{u}}_{bi} \right) \right) \right)  dV ,\notag \\
    \pderiv{(R_{ss}^e)_{ai}}{\hat{u}_{bj}} &= 
        \int\limits_{\Gamma_c^e} \Phi_a \pderiv{c_i}{u_{bj}}   dA \notag \\
        &- \int\limits_{\Omega^e} \Phi_{a,k} \left(\pderiv{T_{ik}^{n+1}}{\hat{u}_{bj}}+  a_{stiff} \frac{\gamma}{\beta \Delta t} \left(\pderiv{T_{ik}^{n+1}}{\hat{u}_{bj}} \right)_{elastic} \right)   dV \notag \\
        &- \left( \frac{\gamma a_m}{\beta \Delta t} + \frac{1}{\beta \Delta t^2}  \right) \int\limits_{\Omega^e} \rho \Phi_a \Phi_c    \pderiv{ \hat{u}_{ci} }{\hat{u}_{bj}}dV 

The choice of initial guess of the update of the displacement/velocity should be somewhere between zero end of step acceleration, and setting the end of step acceleration equal to the beginning of step acceleration.



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