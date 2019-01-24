#########################
Capillary pressure models
#########################

In GEOS-X, the oil pressure is assumed to be the primary pressure.
The following paragraphs explain how the capillary pressure models are used to compute the water-phase and gas-phase pressures as:

.. math::
    p_w = p_o - P_{c,w}(S_w)

and

.. math::
    p_g = p_o + P_{c,g}(S_g)


Brooks-Corey model
==================

In this model, the capillary pressure is computed as a function of saturation with the following expression:

.. math::

   P_c(S_w) = p_e S_{\textit{w,c}} S^{-1/\lambda}

where the scaled water-phase saturation is computed as:

.. math::

   S_{\textit{w,c}} = \frac{S_w - S_{\textit{w,min}} }{1 - S_{\textit{w,min}} - S_{\textit{o,min}} - S_{\textit{g,min} }}

The gas capillary pressure is computed analogously.

This model requires the specification of five types of parameters:

* ``phaseNames``: Defines the list of phases.
The number of phases can be either 2 or 3.
The capillary models assume that oil is always present.
Example: "oil water gas".

* ``phaseMinVolFractions``: Defines the mimimum volume fraction :math:`S_{\ell,min}` for each phase specified in ``phaseNames``, in the same order.
Below this volume fraction, the phase is assumed to be immobile. Default for each phase: 0.0.

* ``phaseCapPressureExponentInv``: Defines the parameter :math:`\lambda` for each phase specified in ``phaseNames``, in the same order.
The parameter corresponding to the oil phase is not used. Default for each phase: 2.0.

* ``phaseEntryPressure``: Defines the parameter :math:`p_e` for each phase specified in ``phaseNames``, in the same order.
The parameter corresponding to the oil phase is not used. Default for each phase: 1.0.

* ``capPressureEpsilon``: Defines the parameter :math:`\epsilon`.
This parameter is used for both the water-phase and gas-phase capillary pressure.
Default: 0.0

Input example:
***************************************************

.. code-block:: xml

<BrooksCoreyCapillaryPressure name="capPressure"
                              phaseNames="oil gas"
                              phaseMinVolumeFraction="0.01 0.015"
                              phaseCapPressureExponentInv="0 6"
                              phaseEntryPressure="0 1e8"
                              capPressureEpsilon="1e-8"/>


