.. _LinearElasticIsotropic:

############################################
Linear elastic isotropic solid model
############################################

Overview
=========================

This model represents a linear elastic isotropic solid.

Usage
=========================

The following attributes are supported:

.. include:: /coreComponents/fileIO/schema/docs/LinearElasticIsotropic.rst

Input example
=========================

.. code-block:: xml

  <Constitutive>
    <LinearElasticIsotropic name="shale"
                            density0="2700"
                            BulkModulus0="61.9e6"
                            ShearModulus0="28.57e6"
                            BiotCoefficient="1"
                            referencePressure="2.125e6"
                            compressibility="3e-10"/>
  </Constitutive>