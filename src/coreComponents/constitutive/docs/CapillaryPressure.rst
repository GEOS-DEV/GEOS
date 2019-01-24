#########################
Capillary pressure models
#########################

This section describes the use of the capillary pressure models implemented in GEOS-X.

In GEOS-X, the oil pressure is assumed to be the primary pressure. 
The following paragraphs explain how the capillary pressure models are used to compute the water-phase and gas-phase pressures as:

.. math::
    p_w = p_o - P_{c,w}(S_w) 

and

.. math::    
    p_g = p_o + P_{c,g}(S_g)

Two analytical models are currently implemented to compute capillary pressure:

* The Brooks-Corey model
* The Van Genuchten model

Each of these models require the specification of five types of parameters. The first type is common to the Brooks-Corey and Van Genuchten models:

* ``phaseNames``: 
defines the list of phases. 
The number of phases can be 2 or 3. 
The capillary models assume that oil is always present. 
Example: "oil water gas".


Brooks-Corey model
==================

In this model, the capillary pressure is computed as a function of saturation with the following expression:

.. math::

   P_c(S_w) = p_e S_{\textit{w,c}} S^{-1/\lambda}

where the scaled water-phase saturation is computed as:

.. math::

   S_{\textit{w,c}} = \frac{S_w - S_{\textit{w,min}} }{1 - S_{\textit{w,min}} - S_{\textit{o,min}} - S_{\textit{g,min} }}

The gas capillary pressure is computed analogously.

In addition to ``phaseNames``, this model requires the specification of four types of parameters:

* ``phaseMinVolFractions``: 
defines the mimimum volume fraction for each phase specified in ``phaseNames``, in the same order.
Below this volume fraction, the phase is assumed to be immobile. Default: "0.0 0.0 0.0".

* ``phaseCapPressureExponentInv``: 
defines the parameter :math:`\lambda` for each phase specified in ``phaseNames``, in the same order. 
The parameter corresponding to the oil phase is not used. Default: "2.0 0.0 2.0".

* ``phaseEntryPressure``:
defines the parameter :math:`p_e` for each phase specified in ``phaseNames``, in the same order.
The parameter corresponding to the oil phase is not used. Default: "1.0 1.0 1.0".

* ``capPressureEpsilon``:
defines the parameter :math:`\epsilon`. 
This parameter is used for both the water-phase and gas-phase capillary pressure. 
Default: "0.0"

Input example:
***************************************************



Van Genuchten model
===================

The Van Genuchten model computes the capillary pressure as a function of saturation as

.. math::

  P_c(S_w) = \alpha  ( S_{w,c}^{-1/m} - 1 )^{ (1-m)/2 }

where the scaled water-phase saturation is computed as:

.. math::

   S_{\textit{w,c}} = \frac{S_w - S_{\textit{w,min}} }{1 - S_{\textit{w,min}} - S_{\textit{o,min}} - S_{\textit{g,min} }}

This model also requires the specification of four parameters

* ``phaseMinVolFractions``
* ``phaseCapPressureExponentInv``
* Z
* U

Input example:
**************************************************************
