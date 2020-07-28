.. _FiniteVolumeDiscretization:

Discretization of Finite Volumes
---------------------------------

In the current implementation of the :ref:`SinglePhaseFlow` and the :ref:`CompositionalMultiphaseFlow`, the discretization is based on a cell-centered finite-volume scheme.


The flux is computed using a Two-Point Flux Approximation (TPFA). 
Considering the general compositional multiphase case, the numerical flux of component :math:`c` at the interface between elements :math:`i` and :math:`j` reads:

.. math::
  F_{c,(ij)} = T_{(ij)} \sum_{\ell} x^{upw}_{c,\ell} \lambda^{upw}_{\ell} \Delta p_{\ell,(ij)},

where :math:`\Delta p_{\ell,(ij)}` is the phase potential difference at the interface. 
The phase component fraction, :math:`x^{upw}_{c,\ell}`, and the phase mobility, :math:`\lambda^{upw}_{\ell}`, are upwinded using the sign of the phase potential difference at the interface (single-point upwinding). 
Finally, :math:`T_{(ij)}` is the standard TPFA transmissibility at the interface.

The discretization of the flux for single-phase flow is analogous and is omitted for brevity.
