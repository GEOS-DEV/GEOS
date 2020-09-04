.. _DruckerPragerPlasticity:

############################################
Drucker-Prager plasticity solid model
############################################

Overview
=========================

This model may be used to represents a solid material with plastic response to loading according to the `Drucker-Prager <https://en.wikipedia.org/wiki/Drucker%E2%80%93Prager_yield_criterion>`__ yield criterion below:

.. math::
   f (p,q) = q + B \, p - A = 0,

where :math:`q = \sqrt{3J_2} = \sqrt{3/2}\norm{\boldsymbol{s}}` is the von Mises stress, :math:`p = I_1/3` is the mean normal stress, and . :math:`I_1 ` and :math:`J_2 ` are  the first invariant of the stress tensor and second invariant of the deviator stress, defined as follows.

.. math::
   I_1 = tr(\tensor{\sigma})/3 \, , \quad J_2 = \frac{1}{2} \|\boldsymbol{s}\|^2 \, , \quad \boldsymbol{s}=\boldsymbol{\sigma}-p \boldsymbol{1}

in which :math:`\boldsymbol{1}` is the identity tensor. Here :math:`f` represents the yield surface, as shown on the figure below.

.. _druckerPragerYield:
.. figure:: ../../../coreComponents/constitutive/docs/DPyield.png
   :align: center
   :width: 500
   :figclass: align-center
   
   Mohr-Coulomb and Drucker-Prager yield surfaces in principal stress axes (Borja, 2002). 
   
The material behavior is linear elastic (see :ref:`LinearElasticIsotropic`) for :math:`f < 0`, and plastic for :math:`f =0`. The two material parameters are :math:`A` and :math:`B` are derived by approximating the Mohr-Coulomb surface with a cone. The Drucker-Prager yield surface has a circular cross-section in deviators plane that passes through the tension or compression corners of the Mohr-Coulomb yield surface, as shown in the Figure 2. The material parameters :math:`A` and :math:`B` are derived as:

.. math::
   A = \frac{6 \, c \, \cos\phi}{3 \pm \sin\phi} \, , \quad B=\frac{6  \, \sin\phi}{3 \pm \sin\phi}
   
where plus signs are for circles passing through the tension corners, and minus signs are for circles passing through tension corners. Also, :math:`\phi` and :math:`c` denote friction angle and cohesion, respectively, as defined by the Mohr-Coulomb failure envelope shown in Figure 3. In this code, :math:`A` and :math:`B` are named ``cohesionParameter`` and ``frictionParameter``, respectively. 

.. _deviatoricView:
.. figure:: ../../../coreComponents/constitutive/docs/DevView.png
   :align: center
   :width: 300
   :figclass: align-center
   
   Mohr-Coulomb and Drucker-Prager yield surfaces on the deviatoric plane (Borja, 2013). 
   
.. _mohrCoulombEnvelope:
.. figure:: ../../../coreComponents/constitutive/docs/MohrCoulomb.png
   :align: center
   :width: 400
   :figclass: align-center
   
   The Mohr-Coulomb failure envelope.
   
We consider a plastic potential in the following form to determine the direction of plastic flow.

.. math::
   g (p,q) = q + b \, p - g_0= 0,

where :math:`g_0` is a constant and :math:`b \leq B` is the dilatancy parameter.  :math:`b = B` corresponds to associative flue rule, while for :math:`b <  B` non-associative flow is obtained. In the code, :math:`b` is named ``dilationParameter`` and is related to `dilation angle <https://en.wikipedia.org/wiki/Dilatancy_(granular_material)>`__  as:

.. math::
   b = \frac{6  \, \sin\psi}{3 \pm \sin\psi},
   
 where :math:`\psi \leq \phi` is the dilation angle. If :math:`\psi > 0`, then the plastic flow is dilative.

The total strain :math:`\boldsymbol{\epsilon}` can be additively split into elastic (:math:`\boldsymbol{\epsilon}^e`) and plastic :math:`\boldsymbol{\epsilon}^p`) strains:

.. math::
   \boldsymbol{\epsilon} = \boldsymbol{\epsilon}^e + \boldsymbol{\epsilon}^p.

The plastic strain tensor is obtained from the flow rule: 

.. math::
   \dot{\boldsymbol{\epsilon}}^p=\dot{\lambda}\frac{\partial g}{\partial\boldsymbol{\sigma}},
   
in which :math:`\dot{\lambda} \geq 0` is the magnitude of plastic strain rate. The elastic strain is related to Cauchy stress tensor in rate form as:

.. math::
  \dot{\boldsymbol\sigma}} = \tensor{C}^e \cdot \dot{\boldsymbol{\epsilon}}^e,

where :math:`\tensor{C}^e` is the fourth order elastic stiffness tensor. The Cauchy stress tensor is related to the total strain as

.. math::
  \dot{\boldsymbol\sigma}} = \boldsymbol{C}^{ep} \cdot \dot{\boldsymbol{\epsilon}},
  
in which :math:`\tensor{C}^{ep}` is the fourth order elasto-plastic stiffness tensor.

A hardening rule is defined which determines how the yield surface will change as a result of plastic deformations. Here we use linear hardening for the cohesion parameter, :math:`A`, 

.. math::
   \dot{A}= h \, \dot{\lambda},
   
where :math:`h` is the hardening parameter. 

Usage
=========================

The following attributes are supported:

.. include:: /coreComponents/fileIO/schema/docs/DruckerPrager.rst

Input example
=========================

.. code-block:: xml

  <Constitutive>
    <DruckerPrager name="shale"
                              defaultDensity="2700"
                              defaultBulkModulus="1.0"
                              defaultShearModulus="1.0"
                              defaultTanFrictionAngle="0.5"
                              defaultTanDilationAngle="0.5"
                              defaultHardeningRate="0.0"
                              defaultCohesion="1.0" />
  </Constitutive>
