.. _CamClayPlasticity:

############################################
Cam-Clay plasticity solid model
############################################

Overview
=========================

This model may be used to represents a solid material with plastic response to loading according to the `Modified Cam-Clay  (MCC) <https://en.wikipedia.org/wiki/Critical_state_soil_mechanics>`__ critical state model .
The implementation also accommodates `Delft-Egg  <https://link.springer.com/chapter/10.1007%2F978-94-011-1046-4_10>`__  constitutive model, which is a generalization of MCC. The yield function for the generalized Cam-Clay model is 

.. math::
  f = q^2 - M^2 \left[ \alpha^2 \, p \left(\frac{2 \alpha}{\alpha+1} p_c -p \right) - \frac{\alpha^2 (\alpha-1)}{\alpha+1} p_c^2 \right] = 0 , 

where :math:`\alpha \geq 1` is the shape parameter. For :math:`\alpha = 1`, the MCC model will be derived with an ellipsoidal yield surface. For :math:`\alpha > 1, the Self-Egg model will be derived with an egg-shaped yield surface. :math:`p_c` is the preconsolidation pressure, and :math:`M` is the slope of the critical state line (CSL).  :math:`M` can be related to the critical state friction angle :math:`\phi_{cs}` as

.. math::
   M = \frac{6 \sin \phi_{cs}}{3-\sin \phi_{cs}}.
   
:math:`q = \sqrt{3J_2} = \sqrt{3/2}\norm{\boldsymbol{s}}` is the von Mises stress, :math:`p = I_1/3` is the mean normal stress, and . :math:`I_1 ` and :math:`J_2 ` are  the first invariant of the stress tensor and second invariant of the deviator stress, defined as follows.

.. math::
   I_1 = tr(\tensor{\sigma})/3 \, , \quad J_2 = \frac{1}{2} \|\boldsymbol{s}\|^2 \, , \quad \boldsymbol{s}=\boldsymbol{\sigma}-p \boldsymbol{1}

in which :math:`\boldsymbol{1}` is the identity tensor. Here :math:`f` represents the yield surface, as shown on the figure below.

.. _CamClaypq:
.. figure:: ../../../coreComponents/constitutive/docs/CamClaypq.png
   :align: center
   :width: 400
   :figclass: align-center
   
   Cam-Clay yield surface in p-q space (Borja, 2013). 
   
The material behavior is linear elastic (see :ref:`LinearElasticIsotropic`) for :math:`f < 0`, and plastic for :math:`f =0`.    
The total strain :math:`\boldsymbol{\epsilon}` can be additively split into elastic (:math:`\boldsymbol{\epsilon}^e`) and plastic :math:`\boldsymbol{\epsilon}^p`) strains:

.. math::
   \boldsymbol{\epsilon} = \boldsymbol{\epsilon}^e + \boldsymbol{\epsilon}^p.

We consider an associative flow rule, where the plastic potential and the yield surface are the same: 

.. math::
   \dot{\boldsymbol{\epsilon}}^p=\dot{\lambda}\frac{\partial f}{\partial\boldsymbol{\sigma}},
   
in which :math:`\dot{\lambda} \geq 0` is the magnitude of plastic strain rate. The elastic strain is related to Cauchy stress tensor in rate form as:

.. math::
  \dot{\boldsymbol\sigma}} = \tensor{C}^e \cdot \dot{\boldsymbol{\epsilon}}^e,

where :math:`\tensor{C}^e` is the fourth order elastic stiffness tensor. The Cauchy stress tensor is related to the total strain as

.. math::
  \dot{\boldsymbol\sigma}} = \boldsymbol{C}^{ep} \cdot \dot{\boldsymbol{\epsilon}},
  
in which :math:`\tensor{C}^{ep}` is the fourth order elasto-plastic stiffness tensor.

Here we use a hyper-elastic constitutive law using the following elastic rate constitutive equation 

.. math::
  \dot{p} = - \frac{p}{c_r} \dot{\epsilon}^e_v,
  
where :math:`c_r > 0` is the elastic compressibility index. The tangential elastic bulk modulus is :math:`K=- \frac{p}{c_r} ` and varies linearly with pressure. We assume a constant shear modulus, and can write stress invariants p and q as

.. math::
  p = p_0 \exp \left( \frac{\epsilon_{v0} - \epsilon_v^e}{c_r}\right) , \quad q = 3 \mu \epsilon_s^e,
  
where :math:`p_0` is the reference pressure and :math:`\epsilon_{v0}` is the reference volumetric strain. The hardening law is derived from the linear relationship between logarithm of specific volume and logarithm of preconsolidation pressure, as show in Figure 2. 

.. _CamClayHardening:
.. figure:: ../../../coreComponents/constitutive/docs/CamClayHardening.png
   :align: center
   :width: 400
   :figclass: align-center
   
   Bilogarithmic hardening law derived from isotropic compression tests  (Borja, 2013). 

The hardening law describes evolution of the preconsolidation pressure :math:`p_c` as

.. math::
  \dot{p_c} = - \frac{tr(\dot{\boldsymbol{\epsilon}}^p)}{c_c-c_r} p_c,

where :math:`c_c` is the virgin compressibility index and we have :math:` 0 < c_r < c_c`.

Usage
=========================

The following attributes are supported:

.. include:: /coreComponents/fileIO/schema/docs/CamClay.rst

Input example
=========================

.. code-block:: xml

  <Constitutive>
    <CamClay name="shale"
                              defaultDensity="2700"
                              defaultRefPressure="-90.0"
                              defaultRefStrainVol="-1e-4"
                              defaultShearModulus="5400.0"
                              defaultPreConsolidationPressure="-90.0"
                              defaultShapeParameter="1.0"
                              defaultCslSlope="1.0"
                              defaultRecompressionIndex="0.018"
                              defaultVirginCompressionIndex="0.13" />
  </Constitutive>
