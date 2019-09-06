.. _ProppantTransport:

#####################################
Proppant Transport Solver
#####################################

Introduction
=========================

This solver is developed for simulating proppant transport in hydraulic fractures. The behavior of proppant transport is described by a continuum formulation. Here we briefly outline the usage, governing equations and numerical implementation of the proppant transport model in GEOSX.

Usage
=========================




Theory
=========================

The following mass balance and constitutive equations are solved inside fractures,

proppant-fluid mixture:
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. math::
   \frac{\partial}{\partial t}(\rho_m) + \boldsymbol{\nabla} \cdot (\rho_m \boldsymbol{u_m}) = 0,

where the proppant-fluid mixture velocity :math:`\boldsymbol{u_m}` is approximated by the Darcy's law as, 

.. math::
   \boldsymbol{u}_m = -\frac{K_f}{\mu_m}(\nabla p - \rho_m \boldsymbol{g}),

and :math:`p` is pressure, :math:`\rho_m` and :math:`\mu_m` are density and viscosity of the mixed fluid , respectively,  and :math:`\boldsymbol{g}` is the gravity vector. The fracture permeability :math:`K_f` is determined based on fracture aperture :math:`a` as

.. math::
   K_f =  \frac{a^2}{12 \mu_m}

   
proppant:
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. math::
   \frac{\partial}{\partial t}(c) + \boldsymbol{\nabla} \cdot (c \boldsymbol{u}_p) = 0,

in which :math:`c` and :math:`\boldsymbol{u}_p` represent the volume fraction and velocity of the proppant particles. 

The density and velocity of the mixed fluid are further expressed as,

.. math::
   \rho_m = (1 - c) \rho_f + c \rho_p,

and

.. math::
   \rho_m \boldsymbol{u}_m = (1 - c) \rho_f \boldsymbol{u}_f + c \rho_p \boldsymbol{u}_p,

in which :math:`\rho_f` and :math:`\boldsymbol{u}_f` are the density and velocity of the carrying fluid, and :math:`\rho_p` is the density of the proppant particles.   


proppant slip velocity:
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The proppant particle and carrying fluid velocities are related by the slip velocity :math:`\boldsymbol{u}_{slip}`,

 .. math::
    \boldsymbol{u}_{slip} = \boldsymbol{u}_p - \boldsymbol{u}_f.

The slip velocity between the proppant and carrying fluid includes gravitational and collisional components, which take account of particle settling and collision effects, respectively.

The gravitational component of the slip velocity is written as a form as

.. math::
    \boldsymbol{u}_{slip} = F(c) \boldsymbol{u}_{settling},


where :math:`\boldsymbol{u}_{settling}` is the settling velocity for a single particle and calculated based on the Stokes drag law,

.. math::
    \boldsymbol{u}_{settling} = ( \rho_p - \rho_f)  \frac{d{_p}^{2}}{18 \mu_f}\boldsymbol{g},

:math:`d_p` is the particle diameter, and :math:`F(c)` is the correction term to account for hindered settling effects as a result of particle-particle interactions,
         
.. math::
    F(c) = e^{-\lambda_s c},

with :math:`\lambda_s` as an empirical constant set to 5.9 by default. 
    
The collisional component of the slip velocity is modeled by defining :math:`\lambda`, the ratio of the particle velocity to the volume averaged mixture velocity as a function of the proppant concentration. From this the particle slip velocity is related to the mixed fluid velocity by,

.. math::
    \boldsymbol{u}_{slip} =  \frac{1-\lambda c}{1 - c} \boldsymbol{u}_{m}

We use a simple expression proposed by Barree & Conway (1995) to define :math:`\lambda` as,

.. math::
    \lambda=  \left[\alpha - |c - c_{slip} |^{\beta} \right]\, 

where :math:`\alpha` and :math:`\beta` are empirical constants, :math:`c_{slip}` is the volume fraction exhibiting the greatest particle slip. By default the model parameters are set to the values given in (Barree & Conway, 1995): :math:`\alpha= 1.27`, :math:`c_{slip} =0.1` and :math:`\beta =  1.5`. This model can be extended to account for the transition to the particle pack as the proppant concentration approaches the jamming transition.

other features of proppant transport:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The GEOSX also has the option to model

a. proppant pack-formation when proppant concentration is above the maximum packing fraction. Under this condition fluid can continue to flow through the settled proppant pack. The pack permeability `K` is determined based on the Kozeny-Carmen relationship:  

.. math::
   K = \frac{(sd_p)^2}{180}\frac{(1-c)^{3}}{c^{2}},

where :math:`c` is the solid fraction for proppant packing, :math:`s` is the sphericity and :math:`d_p` is the particle diameter.

b. proppant bridging when proppant particle size is close to or larger than fracture aperture. The aperture at which bridging occurs, :math:`h_{b}`, is determined from

.. math::
   h_{b} = \begin{cases} (1.0 + \frac{\lambda_{b} - 1.0}{0.17}c) d_p & \quad c < 0.17 \\ \lambda_{b} d_p & \quad c \geq 0.17 \\ \end{cases},   

in which :math:`\lambda_{b}` is the bridging factor and set to 3 by default, and :math:`d_p` is the particle diameter. 
   
c. the effect of the proppant concentration on the viscosity of the bulk fluid. The effective slurry viscosity, :math:`\mu_m`, is calculated based on the Stokes-Einstein model (Einstein, 1906),

.. math::
     \mu_{m} =  \mu_{f}\left (1 + \frac{5}{2}c \right).


Note that continued model development and improvement are underway and will incorporate additional empirical correlations or functions into GEOSX to support the above calculations.       
   
   
Spatial Discretization
=======================

The above governing equations are discretized using a cell-centered two-point flux approximation (TPFA) finite volume method. We use an upwind scheme to approximate proppant transport across cell interfaces.

Temporal Discretization
=======================

An implicit time integration scheme (backward Euler) is employed to solve slurry flow and proppant transport equations. 

Solution Strategy
=======================

The discretized non-linear equations at each time step are solved by the Newton-Raphson method. Each nonlinear iteration step requires the solution of a set of linear algebraic equations.


References
=================================

- R. D. Barree & M. W. Conway. "Experimental and numerical modeling of convective proppant transport", JPT. Journal of petroleum technology, 47(3):216-222, 1995.

- A. Einstein. Eine neue bestimmung der molekuldimensionen. Annalen der Physik, 324(2):289-306, 1906.  
