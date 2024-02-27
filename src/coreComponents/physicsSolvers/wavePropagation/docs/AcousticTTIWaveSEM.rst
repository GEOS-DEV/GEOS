
#####################################
Acoustic TTI Wave Propagation Solvers
#####################################

Available Solvers
=================


Governing equations
-------------------

Rotation Matrix
+++++++++++++++

Tilted transverse anisotropy (TTI) is a generalization of Vertical transverse anisotropy (VTI). It is a rotation of the transverse axis. 
The coordinate system are rotated to recover the VTI equations. 
The tilted medium is parametrized by its tilt angle :math:`\phi` and its azimuth angle :math:`\theta`. 
Coordinates in the rotated system are denoted by :math:`\check{\vtixx}=(\check{x},\check{y},\check{z})`. 
Following \cite{ZhaZhaZha2011}, the rotation matrices are then given by

.. math::
  :label: eq-tti-R

  \begin{pmatrix}
    \check{x}\\
    \check{y}\\
    \check{z}
  \end{pmatrix}=
  \underbrace{\begin{pmatrix}
    \cos\phi\cos\theta & \sin\phi\cos\theta & -\sin\theta\\
    -\sin\phi & \cos\phi & 0\\
    \cos\phi\sin\theta & \sin\phi\sin\theta & \cos\theta\\
  \end{pmatrix}}_{=\ttiR(\theta,\phi)}
  \begin{pmatrix}
    x\\
    y\\
    z
  \end{pmatrix}

In practice however, the angles :math:`\phi` and :math:`\theta` are not the input of GEOS. They are computed using the inline :math:`D_x` and crossline :math:`D_y` dips - which are the inputs in GEOS - and the below formula:

.. math::
  \left\{\begin{aligned}
  \theta &= \atan(\sqrt{(D_x^2 + D_y^2}) \\
  \phi   &= \atan2(D_y, D_x)
  \end{aligned}\right.

In addition, the angle :math:`\phi` satisfies:

* if :math:`(\phi \leq 0)` then :math:`\phi = \phi + 2\pi`
* if :math:`\theta < 0.001 \frac{\pi}{180}` then :math:`\phi = 0`

  * else if :math:`(0 <\phi < \pi)` then :math:`\phi = \phi + \pi`

  *  else if :math:`(\pi \leq\phi \leq 2\pi)` then :math:`\phi = \phi - \pi`


.. math::
  

``AcousticTTIFletcherWavePropagationSEM``
+++++++++++++++++++++++++++++++++++++++++

This solver is the TTI version of :ref:`AcousticVTIFletcherWavePropagationSEM <sec-vti-fletcher>` described in :footcite:`FlecherDuFowler2009`. As for the VTI case this solver has a tunable parameter :math:`\vtif` (or equivalently :math:`\sigma`, see :ref:`paragraph on the quantities <sec-tti-quantities>`) to manage the spurious artefact s-wave. The only change with the VTI equation :eq:`eqVTIFletcher` is the rotation matrix.

.. math::
  :label: eqTTIFletcher

  \left\{
    \begin{aligned}
      \frac{1}{\rho v_p^2} \frac{\partial^2 p}{\partial t^2} &= (1+2\varepsilon) \vtiAxy\ttiR\nabla\cdot\left(\frac{1}{\rho}\vtiAxy\ttiR\nabla p\right) +\vtiAz\ttiR\nabla\frac{1}{\rho}\vtiAz\ttiR\nabla q - (f_{\mathrm{vti}}-1)\vtiAz\ttiR\nabla\cdot\left(\frac{1}{\rho}\vtiAz\ttiR\nabla (p-q)\right) + f,\\
      \frac{1}{\rho v_p^2} \frac{\partial^2 q}{\partial t^2} &= (1+2\delta) \vtiAxy\ttiR\nabla\cdot\left(\frac{1}{\rho}\vtiAxy\ttiR\nabla p\right)      +\vtiAz\ttiR\nabla\frac{1}{\rho}\vtiAz\ttiR\nabla q + (f_{\mathrm{vti}}-1)\vtiAxy\ttiR\nabla\cdot\left(\frac{1}{\rho}\vtiAxy\ttiR\nabla  (p-q)\right) + f.
   \end{aligned}
  \right.


``AcousticTTIZhangWavePropagationSEM``
++++++++++++++++++++++++++++++++++++++

This second one is based on Zhang et. al set of equations of :footcite:`ZhangZhangZhang2011` (see paragraph :ref:`AcousticVTIZhangWavePropagationSEM <sec-vti-zhang>`):

.. math::
  :label: eqTTIZhang

  \left\{
    \begin{aligned}
       \frac{1}{\rho v_p^2} \frac{\partial^2 p}{\partial t^2} &= (1+2\varepsilon) \vtiAxy\ttiR\nabla\cdot\left(\frac{1}{\rho}\vtiAxy\ttiR\nabla p\right) + \sqrt{1+2\delta}\vtiAz\ttiR\nabla\cdot\left(\frac{1}{\rho}\vtiAz\ttiR\nabla q\right) + f,\\
       \frac{1}{\rho v_p^2} \frac{\partial^2 q}{\partial t^2} &= \sqrt{1+2\delta} \vtiAxy\ttiR\nabla\cdot\left(\frac{1}{\rho}\vtiAxy\ttiR\nabla p\right)      +\vtiAz\ttiR\nabla\left(\frac{1}{\rho}\vtiAz\ttiR\nabla q\right) + f.
    \end{aligned}
  \right.


.. _sec-tti-quantities:
Operators and quantities
------------------------

* :math:`p`: unknown pressure wave
* :math:`q`: auxiliary wave
* :math:`v_p [\textrm{m}.\textrm{s}^{-1}]`: pressure wave celerity in the medium.
* :math:`\rho [\textrm{kg}.\textrm{m}^{-3}]`: density of the medium.
* :math:`\varepsilon` and :math:`\delta`: anisotropy Thomsen parameters. Default value is 0 for both.
* :math:`\sigma_{\mathrm{vti}}` and :math:`f_{\mathrm{vti}} = 1 - \frac{\varepsilon - \delta}{\sigma_{\mathrm{vti}}}`: tunable parameters in Fletcher's equation to manage artefact spurious share wave. Even though :math:`\sigma_{\mathrm{vti}}` is the input in GEOS code, the equations are here written in term of :math:`f_{\mathrm{vti}}`. :math:`\sigma_{\mathrm{vti}}` is set by default to 0.75.
* :math:`f`: source term (Rickers)
* :math:`\nabla = [\partial x, \partial y, \partial z]^T`: nabla operator
* :math:`X\cdot Y = X^TY`: scalar product 
* :math:`\vtiAxy` and :math:`\vtiAz`: truncation matrices such that :math:`\vtiAxy\nabla = [\partial_x,\partial_y,0]^T` and :math:`\vtiAz\nabla=[0,0,\partial_z]^T`:

.. math::
  :label: tti-AxyAz

 \vtiAxy = \begin{pmatrix}
  1&0&0\\
  0&1&0\\
  0&0&0
  \end{pmatrix}\quad\text{and}\quad
 \vtiAz = \begin{pmatrix}
  0&0&0\\
  0&0&0\\
  0&0&1
  \end{pmatrix}.

* :math:`\ttiR`: Rotation matrix defined by :eq:`eq-tti-R`
* :math:`\theta`: Tilt angle
* :math:`\phi`: Azimuth angle


Restrictions and Remarks
------------------------

- The computational domain must be a rectangular cuboid. Its boundary :math:`\Gamma` is divided into the lateral surfacess :math:`\Gammaxy` and the top/bottom surfaces :math:`\Gammaz`. The reason is to simply the boundary terms that arise in the weak formulation, as explained :ref:`for the VTI case <sec-vti-rectangularcuboid>`.
- The anisotropic parameters are (currently) constant per element.
- If :math:`\theta` (resp. :math:`\phi`) is very small, it is then set to 0, for numerical stability.
- For numerical stability reasons, the ``AcousticVTIZhangWavePropagationSEM`` solver needs the following conditions to be satisfied:

  1. The :math:`\delta` parameter to be smooth in the domain (no sharp contrast). This must be achieved in the model by the user.
  2. :math:`\varepsilon \geq \delta` everywhere. If GEOS encounters the relation :math:`\delta > \varepsilon`, it will automatically set :math:`\delta = \varepsilon`.
   



Damping methods
---------------

Currently, the only damping method available for TTI is the VTI - ABC described in section :ref:`sec-vti-abc`, the angles :math:`\phi` and :math:`\theta` are assumed to vanish on the boundary :math:`\Gamma`. Hence the VTI ABC is plugged.


Additional ``Fields``
---------------------

The solvers use the exact same parameters as the VTI (see :ref:`sec-vti-fields`) and isotropic ones. The tilt and azimuth angles are not direct input parameters: they are computed from the inline and crossline dips, which are the two new inputs compared to VTI. These inputs are implemented in GEOS as ``Fields``.



.. list-table:: TTI Fields
   :header-rows: 1

   * - Name
     - Manager
     - Default
     - Description
   * - ``acousticDipx``
     - Cell
     - 0
     - Inline Dips
   * - ``acousticDipy``
     - Cell
     - 0
     - Crossine Dips




Mathematical analysis
=====================

The procedure is the same as for the VTI case in paragraph :ref:`sec-vti-math` and is here briefly described for the ``AcousticTTIZhangWavePropagationSEM`` solver. 

Assumptions
-----------

* The anisotropic parameters are assumed to be constant per element. They will hence not be differentiated after the integration by part.
* The domain :math:`\Omega` is assumed to be a rectangular cuboid with boundary :math:`\Gamma = \Gamma_{xy} \bigcup \Gamma_{z}` where :math:`\Gamma_{xy}` is the lateral surface and :math:`\Gamma_{z}` represent the top and bottom surfaces.
* The dips (or tilt and azimuth angles)  are assumed to be constant per element.

.. _sec-tti-rectangularcuboid:
Rectangular Cuboid Domain (default)
-----------------------------------

As for the VTI case, the domain is assumed to be rectangular cuboid even though we do not recover the Neumann boundary condition. The weak formulation reads as

.. math::
  \left\{
      \begin{aligned}
        &\text{Find } p,q\in C^2([0, +\infty])\times H^1(\Omega)  \text{ such that, }\forall p',q'\in H^1(\Omega)\times C^2([0, +\infty]),\\
        &\begin{multlined}[t]
           \int_{\Omega} \frac{1}{\rho v_p} \frac{\partial^2 p}{\partial t^2} p' \diff \mathbf{x} =
          - \int_{\Omega} \frac{(1+2\varepsilon)}{\rho}\vtiAxy\ttiR\nabla p \cdot\vtiAxy\ttiR\nabla p'\diff \mathbf{x}
          - \int_{\Omega} \frac{\sqrt{1+2\delta}}{\rho}\vtiAz\ttiR\nabla q\cdot\vtiAz\ttiR\nabla p'\diff \mathbf{x}\\
          + \int_{\Gamma} \frac{(1+2\varepsilon)}{\rho} (\vtiAxy\ttiR\nabla p)\cdot(\vtiAxy\ttiR\mathbf{n}) p'\diff s 
          + \int_{\Gamma} \frac{1}{\rho} (\vtiAz\ttiR\nabla q)\cdot(\vtiAz\ttiR\mathbf{n})p'\diff s
          + \int_{\Omega} f p'\diff \mathbf{x},
        \end{multlined}\\
      &\begin{multlined}[t]
        \int_{\Omega} \frac{1}{\rho v_p} \frac{\partial^2 q}{\partial t^2} q' \diff \mathbf{x} = 
        - \int_{\Omega}\frac{\sqrt{1+2\delta}}{\rho}\vtiAxy\ttiR\nabla p\cdot\vtiAxy\ttiR\nabla q' \diff \mathbf{x} 
        - \int_{\Omega}\frac{1}{\rho}\vtiAz\ttiR\nabla q\cdot\vtiAz\ttiR\nabla q' \diff \mathbf{x}  \\
        + \int_{\Gamma}\frac{\sqrt{1+2\delta}}{\rho} (\vtiAxy\ttiR\nabla p)\cdot(\vtiAxy\ttiR\mathbf{n})  q' \diff s 
        + \int_{\Gamma}\frac{1}{\rho}(\vtiAz\ttiR\nabla q)\cdot(\vtiAz\ttiR\mathbf{n}) q' \diff s
        + \int_{\Omega}f  q'\diff \mathbf{x},
      \end{multlined}
    \end{aligned}
    \right.

If the tilt angles are set to zero on the boundaries (*i.e.* :math:`\ttiR|_{\Gamma} = \mathbf{I}`) and if the domain is a rectangle cuboid, then we recover Neumann condition on the boundaries, for example

.. math::

  \int_{\Gamma}(\vtiAxy\ttiR\nabla p)\cdot(\vtiAxy\ttiR\mathbf{n})  q' \diff s
  =\int_{\Gammaxy} \dn(p)  q' \diff s


Absorbing Boundary Conditions (ABC)
-----------------------------------

On the boundary :math:`\Gamma`, the TTI angles are assumed to vanish and thus the rotation matrix is the identity :math:`\ttiR|_{\Gamma}=\mathbf{I}`. The VTI ABC :eq:`eq-vti-alpha` is then be applied by default.


Perfectly Matched Layer (PML)
-----------------------------

See VTI section :ref:`sec-vti-plm`.

Currently, this option is not supported by GEOS.


Initial condition
-----------------
See VTI section :ref:`sec-vti-initial`.


Space discretization
--------------------

The unknown :math:`p` and :math:`q` are discretized using spectral element method or order :math:`r` leading to respectively the unknown vectors :math:`\vtipb` and :math:`\vtiqb` of :math:`\vtiVhr`. 
The following matrices are introduced where :math:`\Phi_I` and :math:`\Phi_J` refer to the basis functions associated to the :math:`I^{\textrm{th}}` and :math:`J^{\textrm{th}}` degree of freedom respectively. First, the mass and damping (or mass on the boundary) matrices

.. math::

  \left\{
  \begin{aligned}
    \vtiMass(\beta) &= \left(\vtiMass_{I,J}(\beta)\right)_{I,J},& \vtiMass_{I,J}(\beta) & = \vtiint{\Omega}{\beta(\vtixx)\Phi_J(\vtixx)\Phi_I(\vtixx)}{\vtixx},\\
    \vtiDamp(\beta) &= \left(\vtiDamp_{I,J}(\beta)\right)_{I,J},& \vtiDamp_{I,J}(\beta) & = \vtiint{\Gamma}{\beta(s(\vtixx))\Phi_J(s(\vtixx))\Phi_I(s(\vtixx))}{s},\\
    \vtiDampxy(\beta) &= \left(\vtiDamp^{xy}_{I,J}(\beta)\right)_{I,J},& \vtiDamp^{xy}_{I,J}(\beta) & = \vtiint{\Gammaxy}{\beta(s(\vtixx))\Phi_J(s(\vtixx))\Phi_I(s(\vtixx))}{s},\\
    \vtiDampz(\beta) &= \left(\vtiDamp^{z}_{I,J}(\beta)\right)_{I,J},& \vtiDamp^{z}_{I,J}(\beta) & = \vtiint{\Gammaz}{\beta(s(\vtixx))\Phi_J(s(\vtixx))\Phi_I(s(\vtixx))}{s}.
  \end{aligned}
  \right.

Second, the stiffness and rotated stiffness matrices are defined by

.. math::
  :label: eq-tti-stiff

  \left\{
  \begin{aligned}
    \vtiStiff(\beta) &=  \left(\ttiStiff_{I,J}(\beta)\right)_{I,J},& \vtiStiff_{I,J}(\beta) & = \vtiint{\Omega}{\beta(\vtixx)\nabla \Phi_J(\vtixx)\cdot\nabla\Phi_I(\vtixx)}{\vtixx},\\
    \ttiStiffxy(\beta) &= \left(\ttiStiff^{xy}_{I,J}(\beta)\right)_{I,J},& \ttiStiff^{xy}_{I,J}(\beta) & = \vtiint{\Omega}{\beta(\vtixx)\vtiAxy\ttiR\Phi_J(\vtixx)\cdot\vtiAxy\ttiR\nabla\Phi_I(\vtixx)}{\vtixx},\\
    \ttiStiffz(\beta) &= \left(\ttiStiff^z_{I,J}(\beta)\right)_{I,J},& \ttiStiff^z_{I,J}(\beta) & = \vtiint{\Omega}{\beta(\vtixx)\vtiAz\ttiR\nabla\Phi_J(\vtixx)\cdot\vtiAz\ttiR\nabla \Phi_I(\vtixx)}{\vtixx}.
  \end{aligned}
  \right.



The discretized weak formulation is then given by



.. math::

  \left\{
    \begin{aligned}
      &\begin{multlined}
        \vtiMass\left(\frac{1}{\rho v_p^2}\right)\frac{\partial^2 \vtipb}{\partial t^2}
        + \ttiStiffxy\left(\frac{1+2\varepsilon}{\rho}\right) \vtipb
        + \ttiStiffz\left(\frac{\sqrt{1+2\delta}}{\rho}\right) \vtiqb \\
        + \vtiDampxy\left(\frac{\alpha(1+2\varepsilon)}{\rho}\right)\frac{\partial \vtipb}{\partial t}
        + \vtiDampz\left( \frac{\alpha\sqrt{1+2\delta}}{\rho}\right) \frac{\partial \vtiqb}{\partial t}
        = \vtifbp,
      \end{multlined}\\
      &\begin{multlined}
        \vtiMass\left(\frac{1}{\rho v_p^2}\right) \frac{\partial^2 \vtiqb }{\partial t^2}
        + \ttiStiffxy\left(\frac{\sqrt{1+2\delta}}{\rho}\right) \vtipb
        + \ttiStiffz \left(\frac{1}{\rho}\right) \vtiqb\\
        + \vtiDampxy\left(\frac{\alpha\sqrt{1+2\delta}}{\rho}\right)\frac{\partial \vtipb}{\partial t}
        + \vtiDampz \left(\frac{\alpha                }{\rho}\right)\frac{\partial \vtiqb}{\partial t}
        = \vtifbq.
      \end{multlined}
    \end{aligned}
  \right.


Time discretization
-------------------

The leapfrog scheme with :math:`\vtidt` as a time step leads to the following approximation to compute time step :math:`n+1` from :math:`n` and :math:`n-1`:

.. math::
  \left\{
    \begin{aligned}
      \frac{\partial^2\vtipb}{\partial t^2} &\approx  \frac{\vtipb^{n+1} - 2\vtipb^{n} +\vtipb^{n-1}}{\vtidt^2},\\
      \frac{\partial \vtipb}{\partial t} &\approx \frac{\vtipb^{n+1}- \vtipb^{n-1}}{2 \vtidt},
    \end{aligned}
  \right.
  \quad\text{ and }\quad
  \left\{
    \begin{aligned}
      \frac{\partial^2\vtiqb}{\partial t^2} &\approx  \frac{\vtiqb^{n+1} - 2\vtiqb^{n} +\vtiqb^{n-1}}{\vtidt^2},\\
      \frac{\partial\vtiqb}{\partial t} &\approx \frac{\vtiqb^{n+1}- \vtiqb^{n-1}}{2 \vtidt}.
    \end{aligned}
  \right.

.. _sec-tti-wf-final:
Weak Formulations (Final Form)
------------------------------


The weak formulation for both solvers finally read as followm where the ABC parameter :math:`\alpha` given by equation :eq:`vti-alpha`.


Weak formulation for ``AcousticTTIFletcherWavePropagationSEM``
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

.. math::

  \left\{
    \begin{aligned}
      &\text{Find } \vtipb,\vtiqb\in \mathbb{R}^{d_r}  \text{ such that,}\\
      &\begin{multlined}[t]
        \left[\frac{1}{\vtidt^2}\vtiMass\left(\frac{1}{\rho v_p^2}\right)
        +\frac{1}{2\vtidt}\vtiDampxy\left(\frac{\alpha(1+2\varepsilon)}{\rho}\right)
        - \frac{1}{2\vtidt}\vtiDampz\left(\frac{\alpha(\vtif-1)}{\rho}\right)
        \right]\vtipb^{n+1} 
        + \frac{1}{2\vtidt} \vtiDampz\left(\frac{\alpha \vtif}{\rho}\right) \vtiqb^{n+1}=  
         \frac{2}{\vtidt^2}\vtiMass\left(\frac{1}{\rho v_p^2}\right)  \vtipb^{n} 
        -\frac{1}{ \vtidt^2}   \vtiMass\left(\frac{1}{\rho v_p^2}\right)\vtipb^{n-1} \\
        -\ttiStiffxy(1+2\varepsilon) \vtipb^n 
        + \ttiStiffz(\vtif-1)\vtipb^n
        -\ttiStiffz(\vtif)\vtiqb^n 
        + \left[\frac{1}{2 \vtidt}\vtiDampxy\left(\frac{\alpha(1+2\varepsilon)}{\rho}\right)
        - \frac{1}{2 \vtidt}\vtiDampz\left(\frac{\alpha(\vtif-1))}{\rho}\right)\right]\vtipb^{n-1}
        +\frac{1}{2 \vtidt}\vtiDampz\left(\frac{\alpha \vtif}{\rho}\right)\vtiqb^{n-1}
        + \vtiMass\vtifbp^n,
      \end{multlined}\\
      &\begin{multlined}[t]
        \left[\frac{1}{\vtidt^2}\vtiMass\left(\frac{1}{\rho v_p^2}\right) 
        + \frac{1}{2\vtidt}\vtiDampz \left(\frac{\alpha}{\rho}\right)
        - \frac{1}{2\vtidt}\vtiDampxy\left(\frac{\alpha(\vtif-1)}{\rho}\right)\right]\vtiqb^{n+1}
        + \frac{1}{2\vtidt}\vtiDampxy\left(\frac{\alpha (\vtif + 2\delta)}{\rho}\right) \vtipb^{n+1}
        =  \frac{1}{2\vtidt^2}\vtiMass\left(\frac{1}{\rho v_p^2 }\right) \vtiqb^{n} -\frac{1}{\vtidt^2}\vtiMass\left(\frac{1}{\rho v_p^2 }\right)\vtiqb^{n-1} \\
        - \ttiStiffxy(2\delta+\vtif) \vtipb^n 
        + \ttiStiffxy(\vtif-1) \vtiqb^n
        - \ttiStiffz \vtiqb^n
        +\frac{1}{2 \vtidt}\vtiDampxy\left(\frac{\alpha(2\delta + \vtif)}{\rho}\right) \vtipb^{n-1}
        +\frac{1}{2 \vtidt}\left[- \vtiDampxy\left(\frac{\alpha(\vtif-1)}{\rho}\right)
        +\frac{1}{2 \vtidt}\vtiDampz\left(\frac{\alpha}{\rho}\right) \right]\vtiqb^{n-1}
        + \vtiMass\vtifbq^n.
      \end{multlined}
    \end{aligned}
  \right.





Weak formulation for ``AcousticVTIZhangWavePropagationSEM``
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


.. math::
  
  \left\{\begin{aligned}
    &\begin{multlined}
    \left[\frac{1}{\vtidt^2}\vtiMass\left(\frac{1}{\rho v_p^2}\right) 
    + \frac{1}{2\vtidt}\vtiDampxy\left(\frac{\alpha(1+2\varepsilon)}{\rho}\right)\right]\vtipb^{n+1} 
    + \frac{1}{2\vtidt}\vtiDampz\left(\frac{\alpha\sqrt{1+2\delta}}{\rho}\right) \vtiqb^{n+1} =
     \left[\frac{2}{\vtidt^2}\vtiMass\left(\frac{1}{\rho v_p^2}\right) 
    - \ttiStiffxy\left(\frac{1+2\varepsilon  }{\rho}\right)\right] \vtipb^n\\
    - \ttiStiffz \left(\frac{\sqrt{1+2\delta}}{\rho}\right) \vtiqb^n 
    + \left[-\frac{1}{\vtidt^2}\vtiMass\left(\frac{1}{\rho v_p^2}\right)
    + \frac{1}{2\vtidt}\vtiDampxy\left(\frac{\alpha(1+2\varepsilon)}{\rho}\right)\right] \vtipb^{n-1}
    + \frac{1}{2\vtidt}\vtiDampz \left(\frac{\alpha\sqrt{1+2\delta}}{\rho}\right) \vtiqb^{n-1}
    + \vtifbp,
    \end{multlined}\\
    &\begin{multlined}
  \left[\frac{1}{\vtidt^2}\vtiMass\left(\frac{1}{\rho v_p^2}\right) 
    + \frac{1}{2\vtidt}\vtiDampz \left(\frac{\alpha                }{\rho}\right)\right]\vtiqb^{n+1}
    + \frac{1}{2\vtidt}\vtiDampxy\left(\frac{\alpha\sqrt{1+2\delta}}{\rho}\right) \vtipb^{n+1} =
    \left[\frac{2}{\vtidt^2}\vtiMass\left(\frac{1}{\rho v_p^2}\right)
    - \ttiStiffz \left(\frac{1               }{\rho}\right)\right]\vtiqb^n \\
    - \ttiStiffxy\left(\frac{\sqrt{1+2\delta}}{\rho}\right)\vtipb^n
    + \left[-\frac{1}{\vtidt^2}\vtiMass\left(\frac{1}{\rho v_p^2}\right)
    +\frac{1}{2\vtidt}\vtiDampz  \left(\frac{\alpha                }{\rho}\right) \right] \vtiqb^{n-1}
    + \frac{1}{2\vtidt}\vtiDampxy\left(\frac{\alpha\sqrt{1+2\delta}}{\rho}\right)\vtipb^{n-1}
    + \vtifbq.
    \end{multlined}
    \end{aligned}\right.





Implementation
==============

Tilted Stiffness Matrices
-------------------------

Following section :ref:`sec-vti-gen-stiff`, the tilted (or Rotated) stiffness matrices :math:`\ttiStiffxy` and :math:`\ttiStiffz` can be computed as for the VTI case:

.. math::
  \left\{\begin{aligned}
    \ttiStiffpxy &= \int_{\hat{K}} \left(\ttiBpxy(\hat{\vtixx})\nabla\hat{\Phi}_{\vtijb}(\hat{\vtixx})\right) \cdot\nabla\hat{\Phi}_{\vtiib}(\hat{\vtixx})\diff\hat{\vtixx},\\
    \ttiStiffpz &= \int_{\hat{K}} \left(\ttiBpz(\hat{\vtixx})\nabla\hat{\Phi}_{\vtijb}(\hat{\vtixx})\right) \cdot\nabla\hat{\Phi}_{\vtiib}(\hat{\vtixx})\diff\hat{\vtixx},
  \end{aligned}\right.

where the matrices  :math:`\ttiBpxy` and :math:`\ttiBpz` are given by

.. math::
  \left\{\begin{aligned}
  \ttiBpxy(\hat{\vtixx}) &= \abs{\vtijacp(\hat{\vtixx})}\vtiJacp^{-1}(\hat{\vtixx})\ttiR(\vtixx)^T\vtiAxy\ttiR(\vtixx)\vtiJacp^{-T}(\hat{\vtixx}),\\
   \ttiBpz(\hat{\vtixx}) &= \abs{\vtijacp(\hat{\vtixx})}\vtiJacp^{-1}(\hat{\vtixx})\ttiR(\vtixx)^T\vtiAz \ttiR(\vtixx)\vtiJacp^{-T}(\hat{\vtixx}).
  \end{aligned}\right.

As previously, the matrix :math:`\vtiJacp` is the Jacobian matrix of the transformation from the reference hexadron :math:`\hat{K}` to the current element :math:`K_p` and :math:`\vtijacp = \det(\vtiJacp)`. The computation of :math:`\ttiStiffxy` and :math:`\ttiStiffz` should then be about the same complexity as for :math:`\vtiStiffxy` and :math:`\vtiStiffz`. As a result, TTI numerical simulations should be of the same complexity as the VTI ones. Finally, the matrix-matrix products :math:`\ttiR(\vtixx)^T\vtiAxy\ttiR(\vtixx)` and :math:`\ttiR(\vtixx)^T\vtiAz\ttiR(\vtixx)` could be done analytically to reduce the computational cost.

Interior vs Boundary Nodes
--------------------------

See paragraph :ref:`sec-vti-boundary-nodes` for the VTI case.


.. footbibliography::
