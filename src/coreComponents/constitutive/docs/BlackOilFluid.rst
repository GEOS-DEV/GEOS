.. _BlackOilFluid:

############################################
Black-oil fluid model
############################################

Overview
=========================

In the black-oil model three pseudo-components, oil (o), gas (g) and water (w)
are considered. These are assumed to be partitioned across three fluid phases,
named liquid (l), vapour (v) and aqueous (a).

* :math:`B_o`: oil formation volume factor
* :math:`B_g`: gas formation volume factor
* :math:`R_s`: gas/oil ratio
* :math:`r_s`: oil/gas ratio

Phase behavior is characterized by black-oil tables, that tabulate saturated and
undersaturated oil and gas properties as functions of pressure and solution ratios.

Dead oil
-------------
In **dead-oil** each component occupies only one phase. Thus,the following partition matrix determines the components distribution within the
three phases:

.. math::
    \begin{bmatrix}
    y_{gv} & y_{gl} & y_{ga}\\
    y_{ov} & y_{ol} & y_{oa}\\
    y_{wv} & y_{wl} & y_{wa}
    \end{bmatrix}
    = \begin{bmatrix}
    1 & 0 & 0 \\
    0  & 1 & 0 \\
    0 & 0 & 1
    \end{bmatrix}

.. math::
      \rho_{l} = & \, \frac{\rho_{o}^{STC}}{B_{o}} \\
      \rho_{v} = & \, \frac{\rho_{g}^{STC}}{B_{g}}

Live oil
-------------
The live oil fluid model make no assumptions about the partitioning of the
hydrocarbon components and the following composition matrix can be used

.. math::
    \begin{bmatrix}
    y_{gv} & y_{gl} & y_{ga}\\
    \\
    y_{ov} & y_{ol} & y_{oa}\\
    \\
    y_{wv} & y_{wl} & y_{wa}
    \end{bmatrix}
    = \begin{bmatrix}
    \frac{\rho_{g}^{STC}}{\rho_{g}^{STC} + \rho_{o}^{STC} r_{s}} & \frac{\rho_{g}^{STC} R_{s}}{\rho_{o}^{STC} + \rho_{g}^{STC} R_{s}} & 0 \\
    \\
    \frac{\rho_{o}^{STC} r_{s}}{\rho_{g}^{STC} + \rho_{o}^{STC} r_{s}} & \frac{\rho_{o}^{STC}}{\rho_{o}^{STC} + \rho_{g}^{STC} R_{s}} & 0 \\
    \\
    0 & 0 & 1
    \end{bmatrix}

.. math::
  \rho_{l} = & \, \frac{\rho_{o}^{STC} + \rho_{g}^{STC} R_{s}}{B_{o}} \\
  \rho_{v} = & \, \frac{\rho_{g}^{STC} + \rho_{o}^{STC} r_{s}}{B_{g}}

See `Petrowiki`_ for more information.

Usage
=========================

Both types are represented by ``<BlackOilFluid>`` node in the input.
Under the hood this is a wrapper around ``PVTPackage`` library, which is included as a submodule.
In order to use the model, GEOSX must be built with ``-DENABLE_PVTPACKAGE=ON`` (default).

The following attributes are supported:

.. include:: ../coreComponents/fileIO/schema/docs/BlackOilFluid.rst

Supported phase names are:

* oil
* gas
* water

Input example
=========================

.. code-block:: xml

  <Constitutive>
    <BlackOilFluid name="fluid1"
                   fluidType="LiveOil"
                   phaseNames="oil gas water"
                   surfaceDensities="800.0 0.9907 1022.0"
                   componentMolarWeight="114e-3 16e-3 18e-3"
                   tableFiles="pvto.txt pvtg.txt pvtw.txt"/>
  </Constitutive>


.. _Petrowiki: https://petrowiki.org/Phase_behavior_in_reservoir_simulation#Black-oil_PVT_models
