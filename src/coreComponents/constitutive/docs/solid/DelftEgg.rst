.. _DelftEgg:

=========================
Model: Delft Egg
=========================
The `Delft-Egg  <https://link.springer.com/chapter/10.1007%2F978-94-011-1046-4_10>`__  plasticity model is a generalization of the Modified Cam-Clay model, with the yield function defined as 

.. math::
  f = q^2 - M^2 \left[ \alpha^2 \, p \left(\frac{2 \alpha}{\alpha+1} p_c -p \right) - \frac{\alpha^2 (\alpha-1)}{\alpha+1} p_c^2 \right] = 0 , 

where :math:`\alpha \geq 1` is the shape parameter. For :math:`\alpha = 1` leads to a MCC model with an ellipsoidal yield surface. For :math:`\alpha > 1`, the Delft-Egg model with an egg-shaped yield surface will result. 
The rest of the implementation follows the same procedure as the Modified Cam-Clay model. 

In GEOSX we use a single implementation which is called Modified Cam-Clay and can accommodate both MCC and Delft-Egg models. 

Parameters
~~~~~~~~~~~~~~~~~~~~

The following attributes are supported:

.. include:: /coreComponents/fileIO/schema/docs/ModifiedCamClay.rst

Example
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: xml

  <Constitutive>
    <ModifiedCamClay name="delft"
                     defaultDensity="2700"
                     defaultRefPressure="-1.0"
                     defaultRefStrainVol="0"
                     defaultShearModulus="200.0"
                     defaultPreConsolidationPressure="-1.5"
                     defaultCslSlope="1.2"
                     defaultShapeParameter="1.1"
                     defaultRecompressionIndex="0.002"
                     defaultVirginCompressionIndex="0.003" />
  </Constitutive>

