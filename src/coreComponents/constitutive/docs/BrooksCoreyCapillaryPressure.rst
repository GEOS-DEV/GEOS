.. _BrooksCoreyCapillaryPressure:

#####################################
Brooks-Corey capillary pressure model
#####################################

********
Overview
********

In GEOSX, the oil-phase pressure is assumed to be the primary pressure.
The following paragraphs explain how the Brooks-Corey capillary pressure
model is used to compute the water-phase and gas-phase pressures as:

.. math::
    p_w = p_o - P_{c,w}(S_w)

and

.. math::
    p_g = p_o + P_{c,g}(S_g)

In the Brooks-Corey model, the water-phase capillary pressure
is computed as a function of the water-phase volume fraction with
the following expression:

.. math::
   P_c(S_w) = p_{e,w} S_{\textit{w,c}} S^{-1/\lambda_w}

where the scaled water-phase volume fraction is computed as:

.. math::
   S_{\textit{w,c}} = \frac{S_w - S_{\textit{w,min}} }{1 - S_{\textit{w,min}} - S_{\textit{o,min}} - S_{\textit{g,min} }}

The gas capillary pressure is computed analogously.

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

Besides ``name``, this model requires the specification of five
types of parameters:

* ``phaseNames`` - The list of phase names. The number of phases can be either 2 or 3. The capillary model assumes that oil is always present. Example: "oil water gas".

* ``phaseMinVolFraction`` - The list of minimum volume fractions :math:`S_{\ell,min}` for each phase specified in ``phaseNames``, in the same order. Below this volume fraction, the phase is assumed to be immobile. Default value for each phase: 0.0.

* ``phaseCapPressureExponentInv`` - The list of exponents :math:`\lambda_{\ell}` for each phase specified in ``phaseNames``, in the same order. The parameter corresponding to the oil phase is not used. Default value for each phase: 2.0.

* ``phaseEntryPressure`` - The list of entry pressures :math:`p_{e,\ell}` for each phase specified in ``phaseNames``, in the same order. The parameter corresponding to the oil phase is not used. Default value for each phase: 1.0.

* ``capPressureEpsilon`` - The parameter :math:`\epsilon`. This parameter is used for both the water-phase and gas-phase capillary pressure. To avoid extremely large, or infinite, capillary pressure values, we set :math:`P_{c,w}(S_w) := P_{c,w}(\epsilon)` whenever :math:`S_w < \epsilon`. The gas-phase capillary pressure is treated analogously. Default value: 0.0

*************
Input example
*************

.. code-block:: xml

   <Constitutive>
      ...
      <BrooksCoreyCapillaryPressure name="capPressure"
                                    phaseNames="oil gas"
                                    phaseMinVolumeFraction="0.01 0.015"
                                    phaseCapPressureExponentInv="0 6"
                                    phaseEntryPressure="0 1e8"
                                    capPressureEpsilon="1e-8"/>
      ...
   </Constitutive>
