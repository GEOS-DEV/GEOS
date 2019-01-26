.. _BrooksCoreyRelativePermeability:

############################################
Brooks-Corey relative permeability model
############################################

********
Overview
********

The following paragraphs explain how the Brooks-Corey
model is used to compute the phase relative permeabilities as a function
of volume fraction with the expression:

.. math::
    k_{r, m} = k_{\textit{rm,max}} S_{c,m}^{\lambda_{m}}

where the scaled volume fraction of phase :math:`\ell` is computed as:

.. math::

   S_{\textit{m,c}} = \frac{S_{m} - S_{\textit{m,min}} }{1 - S_{\textit{w,min}} - S_{\textit{o,min}} - S_{\textit{g,min}} }

****************
Model parameters
****************

The relative permeability constitutive model is listed in
``<Constitutive>`` block of the input XML file.
The relative permeability model must be assigned a unique name via
``name`` attribute.
This name is used to assign the model to regions of the physical
domain via a ``materialList`` attribute of the ``<ElementRegion>``
node.

Besides ``name``, this model requires the specification of four
types of parameters:

* ``phaseNames`` - The list of phase names. The number of phases can be either 2 or 3. Example: "oil water gas".

* ``phaseMinVolFraction`` - The list of minimum volume fractions :math:`S_{m,min}` for each phase specified in ``phaseNames``, in the same order. Below this volume fraction, the phase is assumed to be immobile. Default value for each phase: 0.0.

* ``phaseRelPermExponent`` - The list of exponents :math:`\lambda_{m}` for each phase specified in ``phaseNames``, in the same order. Default value for each phase: 1.0.

* ``phaseMaxValue`` - The list of maximum values :math:`k_{\textit{rm,max}}` for each phase specified in ``phaseNames``, in the same order. Default value for each phase: 0.0.


**************
Input example
**************

.. code-block:: xml

   <Constitutive>
      ...
      <BrooksCoreyRelativePermeability name="relPerm"
                                       phaseNames="oil water"
                                       phaseMinVolumeFraction="0.02 0.015"
                                       phaseRelPermExponent="2 2.5"
                                       phaseMaxValue="0.8 1.0"/>
      ...
   </Constitutive>
