.. _BrooksCoreyRelativePermeability:

############################################
Brooks-Corey relative permeability model
############################################


Overview
======================

The following paragraphs explain how the Brooks-Corey
model is used to compute the phase relative permeabilities as a function
of volume fraction with the expression:

.. math::
    k_{r\ell} = k_{\textit{r}\ell,\textit{max}} S_{\ell,\textit{scaled}}^{\lambda_{\ell}},

where the scaled volume fraction of phase :math:`\ell` is computed as:

.. math::
   S_{\ell,\textit{scaled}} = \frac{S_{\ell} - S_{\ell,\textit{min}} }{1 - S_{\textit{w,min}} - S_{\textit{o,min}} - S_{\textit{g,min}} }.

The minimum phase volume fractions :math:`S_{\ell,\textit{min}}` are model parameters specified by the user.

Parameters
======================

The relative permeability constitutive model is listed in
``<Constitutive>`` block of the input XML file.
The relative permeability model must be assigned a unique name via
``name`` attribute.
This name is used to assign the model to regions of the physical
domain via a ``materialList`` attribute of the ``<ElementRegion>``
node.

The following attributes are supported:

.. include:: /coreComponents/fileIO/schema/docs/BrooksCoreyRelativePermeability.rst

Below are some comments on the model parameters.

* ``phaseNames`` - The number of phases can be either 2 or 3. The capillary pressure model assumes that oil is always present. Supported phase names are:

===== ===========
Value Phase
===== ===========
oil   Oil phase
gas   Gas phase
water Water phase
===== ===========

* ``phaseMinVolFraction`` - The list of minimum volume fractions :math:`S_{\ell,min}` for each phase is specified in the same order as in ``phaseNames``. Below this volume fraction, the phase is assumed to be immobile.

* ``phaseRelPermExponent`` - The list of exponents :math:`\lambda_{\ell}` for each phase is specified in the same order as in ``phaseNames``.

* ``phaseMaxValue`` - The list of maximum values :math:`k_{\textit{r} \ell,\textit{max}}` for each phase is specified in the same order as in ``phaseNames``.


Example
=======================

.. code-block:: xml

   <Constitutive>
      ...
      <BrooksCoreyRelativePermeability name="relPerm"
                                       phaseNames="oil water"
                                       phaseMinVolumeFraction="0.02 0.015"
                                       phaseRelPermExponent="2 2.5"
                                       phaseRelPermMaxValue="0.8 1.0"/>
      ...
   </Constitutive>
