.. _AdvancedExampleNonLinearThermalOedometric:


######################################################################
Non-linear thermal oedometric example with temperature dependent TEC
######################################################################

------------------------------------------------------------------
Problem description
------------------------------------------------------------------

This example uses the temperature dependent Thermal Expansion Coefficient (TEC) constitutive model to model a non-linear thermo-elastic oedometric problem. To mimic this specific problem, thermal convection and fluid flow are neglected by setting fluid pressure and fluid heat capacity to zero. With a uniform temperature applied on the whole medium that increase linearly from 0 to 90 °C during a period of time from 0 to :math:`10^5` (s).

The confining stress variation :math:`dP_{conf}` w.r.t. temperature change :math:`dT` can be defined by the following tangent relationship:

.. math::
   dP_{conf} = \frac{E}{1-\nu} \alpha(T) dT

where :math:`E` and :math:`\nu`are the elastic Young's modulus and Poisson's ratio. :math:`\alpha(T)` is the thermal expansion coefficient (TEC) that is temperature dependent. 

Let us consider the particular case of a linear TEC: :math:`\alpha(T) = \alpha_0 + A (T - T_0)` where :math:`\alpha_0` is the TEC at a reference temperature :math:`T_0` and :math:`A` is the constant gradient of TEC w.r.t. the temperature. The confining stress can be obtained by integrating the incremental equation above such as:

.. math::
   P_{conf} = \frac{E}{1-\nu} \frac{\alpha_0 + \alpha(T}{2} (T - T_0)

This analytical results will be used to validate GEOS results.


**Input file**

This benchmark example uses no external input file and everything required is
contained within two GEOS xml files that are located at:

.. code-block:: console

  inputFiles/thermoPoromechanics/ThermoElasticOedometric_base.xml

and

.. code-block:: console

  inputFiles/thermoPoromechanics/ThermoElasticOedometric_benchmark.xml

The corresponding integrated test is

.. code-block:: console

  inputFiles/thermoPoromechanics/ThermoElasticOedometric_LinearTEC_smoke.xml

In this example, we would focus our attention on the ``Constitutive`` and ``FieldSpecifications`` tags.

-----------------------------------------------------------
Constitutive
-----------------------------------------------------------

The thermo-elastic properties are defined in the following xml block:  

.. literalinclude:: ../../../../../../../inputFiles/thermoPoromechanics/ThermoElasticOedometric_base.xml
  :language: xml
  :start-after: <!-- SPHINX_TemperatureBC -->
  :end-before: <!-- SPHINX_TemperatureBCEnd -->

--------------------------------------------------------------------   
FieldSpecifications
--------------------------------------------------------------------

The applied uniform temperature is defined in following xml block:

.. literalinclude:: ../../../../../../../inputFiles/thermoPoromechanics/ThermoElasticOedometric_base.xml
  :language: xml
  :start-after: <!-- SPHINX_FieldSpecificationImposedTemperature -->
  :end-before: <!-- SPHINX_FieldSpecificationImposedTemperatureEnd -->

in which the temperature table is defined by

.. literalinclude:: ../../../../../../../inputFiles/thermoPoromechanics/ThermoElasticOedometric_base.xml
  :language: xml
  :start-after: <!-- SPHINX_TemperatureTable -->
  :end-before: <!-- SPHINX_TemperatureTableEnd -->

---------------------------------
Results and benchmark
---------------------------------

A good agreement between the GEOS results and analytical results is shown in the figure below:


.. plot:: docs/sphinx/advancedExamples/validationStudies/thermoPoromechanics/nonLinearThermalOedometric/plot.py


------------------------------------------------------------------
To go further
------------------------------------------------------------------

**Feedback on this example**

This concludes the example of pure thermal diffusion problem around a wellbore.
For any feedback on this example, please submit a `GitHub issue on the project's GitHub page <https://github.com/GEOS-DEV/GEOS/issues>`_.
