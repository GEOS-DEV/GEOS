.. _ParallelPlatesPermeability:


####################################################
Parallel Plates Permeability Model
####################################################


Overview
======================

The parallel plates permeability model defines the relationship between the hydraulic
fracture aperture and its corresponding permeability following the classic lubrication model (`Witherspoon et al.`_ ) .
In this model, the two fracture walls are assumed to be smooth and parallel to each other and separated by a uniform aperture.

.. math::
   k =  \frac{a^3}{12}

where :math:`a` denotes the hydraulic fracture aperture.

Remark: :math:`k`, dimensionally, is not a permeability (as it is expressed in :math:`m^3`).


Joint roughness coefficient
----------------------------

Optionally, the effect of the roughness of the fracture walls can be accounted for through the
joint roughness coefficient (JRC), following the empirical relation of `Barton et al. (1985)`_.
The aperture :math:`E` provided to the model is then replaced by a reduced aperture :math:`e_h`:

.. math::
   e_h = \frac{E^2}{\text{JRC}^{2.5}}

where both apertures are expressed in microns (:math:`\mu m`).

- The JRC must be between 1 and 20 (from smooth to very rough walls); non-negative values outside this range raise an input error.
- The reduction is only applied where :math:`e_h \le E`, i.e. where :math:`E^{0.4} \le \text{JRC}` (:math:`E` in :math:`\mu m`).
  For larger apertures, :math:`E` is used directly.
- The correction is enabled by setting the ``jointRoughnessCoefficient`` attribute.
  If it is not specified or negative (default: -1), the correction is disabled.

For example, with :math:`\text{JRC} = 10` and :math:`E = 100 \mu m`, the reduced aperture is
:math:`e_h = 100^2 / 10^{2.5} \approx 31.6 \mu m`. With :math:`\text{JRC} = 5`, the same aperture
would give :math:`e_h \approx 179 \mu m > E`, so :math:`E` is kept.



Parameters
======================

The Parallel Plates Permeability Model can be called in the
``<Constitutive>`` block of the input XML file.
This permeability model must be assigned a unique name via the
``name`` attribute.
This name is used to assign the model to regions of the physical
domain via a ``materialList`` attribute in the ``<ElementRegions>``
block.

The following attributes are supported:

.. include:: /docs/sphinx/datastructure/ParallelPlatesPermeability.rst


.. _Witherspoon et al.: https://agupubs.onlinelibrary.wiley.com/doi/pdf/10.1029/WR016i006p01016

.. _Barton et al. (1985): https://doi.org/10.1016/0148-9062(85)93227-9
