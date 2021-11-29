.. _ElasticIsotropicPressureDependent:

############################################
Model: Elastic Isotropic Pressure Dependent
############################################

Overview
=========================

This model may be used for solid materials with a pressure-dependent elastic isotropic behavior.
The relationship between stress and strain is given by a `hyperelastic law <https://en.wikipedia.org/wiki/Hyperelastic_material>`__. The elastic constitutive equations for the volumetric and deviatoric stresses and strain are expressed as:

.. math::
   p = p_0 \exp{\left(\frac{\epsilon_{v0}-\epsilon_v^e}{c_r}\right)} \, , \quad q = 3 \mu \epsilon_s^e
   
where :math:`p` and  :math:`q` are the volumetric and deviatoric components of the Cauchy stress tensor.
:math:`\epsilon_{v}^e` and :math:`\epsilon_{s}^e` are the volumetric and deviatoric components of the strain tensor. :math:`\epsilon_{v0}` and :math:`p_0` are the initial volumetric strain and initial pressure. :math:`C_r` denotes the elastic compressibility index,
and :math:`\mu` is the elastic shear modulus. In this model, the shear modulus is constant and the bulk modulus, :math:`K`, varies linearly with pressure as follows: 

.. math::
   K = -\frac{p}{c_r}

Parameters
=========================

The following attributes are supported.  Note that two elastic constants :math:`c_r` and :math:`\mu`, as well as the initial volumetric strain and initial pressure need to be provided.  The "default" keyword in front of certain properties indicates that this
is the default value adopted for a region unless the user separately specifies a heterogeneous field via the
``FieldSpecification`` mechanism. 


Example
=========================

A typical ``Constititutive`` block will look like:

.. code-block:: xml

  <Constitutive>
    <ElasticIsotropicPressureDependent
      name="elasticPressure"
      defaultDensity="2700"
      defaultRefPressure="-1.0"
      defaultRefStrainVol="1"
      defaultRecompressionIndex="0.003"
      defaultShearModulus="200"/>
  </Constitutive>
