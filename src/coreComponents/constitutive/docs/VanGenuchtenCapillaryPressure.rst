Van Genuchten capillary pressure model
######################################

********
Overview
********

In GEOSX, the oil-phase pressure is assumed to be the primary
pressure.
The following paragraphs explain how the Van Genuchten capillary
pressure model is used to compute the water-phase and gas-phase
pressures as:

.. math::
    p_w = p_o - P_{c,w}(S_w)

and

.. math::
    p_g = p_o + P_{c,g}(S_g)


The Van Genuchten model computes the water-phase capillary
pressure as a function of the water-phase saturation as:

.. math::

  P_c(S_w) = \alpha_w  ( S_{w,c}^{-1/m_w} - 1 )^{ (1-m_w)/2 }

where the scaled water-phase saturation is computed as:

.. math::

   S_{\textit{w,c}} = \frac{S_w - S_{\textit{w,min}} }{1 - S_{\textit{w,min}} - S_{\textit{o,min}} - S_{\textit{g,min} }}

The gas-phase capillary pressure is computed analogously.

****************
Model parameters
****************

The capillary pressure constitutive model is listed in
``<Constitutive>`` block of the input XML file.
The capillary pressure model must be assigned a unique name via
``name`` attribute.
This name is used to assign the model to regions of the physical
domain via a ``materialList`` attribute of the ``<ElementRegion>``
node.

Besides ``name``, the model also requires the specification
of five parameters:

   * ``phaseNames`` - The list of phases. The number of phases can be either 2 or 3. The capillary model assumes that oil is always present. Example: "oil water gas".

   * ``phaseMinVolFraction`` - The mimimum volume fraction :math:`S_{\ell,min}` for each phase specified in ``phaseNames``, in the same order. Below this volume fraction, the phase is assumed to be immobile. Default value for each phase: 0.0.

   * ``phaseCapPressureExponentInv`` - The parameter :math:`m_{\ell}` for each phase specified in ``phaseNames``, in the same order. The parameter corresponding to the oil phase is not used. Default value for each phase: 0.5.

   * ``phaseCapPressureMultiplier`` - The parameter :math:`\alpha_{\ell}` for each phase specified in ``phaseNames``, in the same order. The parameter corresponding to the oil phase is not used. Default value for each phase: 1.

   * ``capPressureEpsilon`` - The parameter :math:`\epsilon`. This parameter is used for both the water-phase and gas-phase capillary pressure. To avoid extremely large, or infinite, capillary pressure values, we set :math:`P_{c,w}(S_w) := P_{c,w}(\epsilon)` whenever :math:`S_w < \epsilon`. The gas-phase capillary pressure is treated analogously. Default value: 0.0

**************
Input example
**************

.. code-block:: xml

    <Constitutive>
       ...
       <VanGenuchtenCapillaryPressure name="capPressure"
                                      phaseNames="water oil"
                                      phaseMinVolumeFraction="0.1 0.015"
                                      phaseCapPressureExponentInv="0.55 0"
                                      phaseEntryPressure="1e6 0"
                                      capPressureEpsilon="1e-7"/>
      ...
    </Constitutive>
