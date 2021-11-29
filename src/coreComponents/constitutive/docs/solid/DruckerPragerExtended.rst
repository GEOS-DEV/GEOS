.. _DruckerPragerExtended:

############################################
Model: Extended Drucker-Prager
############################################

Overview
-------------------

This model implements a more sophisticated version of the Drucker-Prager model (see :ref:`DruckerPrager`) allowing for both
cohesion and friction hardening / softening. 
We implement the specific hardening model reported in `Liu et al. (2020) <https://doi.org/10.1007/s00603-019-01992-5>`_.
The yield surface is given by

.. math::
   f(p,q) = q + b \left( p - \frac{a_i}{b_i} \right) = 0,

where :math:`b` is the current yield surface slope, :math:`b_i` is the initial slope, and :math:`a_i` is the
initial cohesion intercept in p-q space.  
The vertex of the Drucker-Prager cone is fixed at :math:`p=a_i/b_i`.
Let :math:`\lambda` denote the accumulated plastic strain measure.  The current yield surface slope is given by
the hyperbolic relationship

.. math::
   b = b_i + \frac{\lambda}{m+\lambda} \left( b_r - b_i \right)

with :math:`m` a parameter controlling the hardening rate.  Here, :math:`b_r` is the residual yield surface slope.  
If :math:`b_r < b_i`, hardening behavior will be observed, while for :math:`b_r < b_i` softening behavior will occur.

In the resulting model, the yield surface begins at an initial position defined by the initial cohesion and friction angle.
As plastic deformation occurs, the friction angle hardens (or softens) so that it asymptoptically approaches a
residual friction angle.  The vertex of the cone remains fixed in p-q space, but the cohesion intercept evolves in
tandem with the friction angle.  See `Liu et al. (2020) <https://doi.org/10.1007/s00603-019-01992-5>` for complete details.

In order to allow for non-associative behavior, we define a "dilation ratio" parameter :math:`\theta \in [0,1]` such
that :math:`b' = \theta b`, where :math:`b'` is the slope of the plastic potential surface, while :math:`b` is 
the slope of the yield surface.  Choosing :math:`\theta=1` leads to associative behavior, while :math:`\theta=0`
implies zero dilatancy.

Parameters
~~~~~~~~~~~~~~~~~~~~

The supported attributes will be documented soon.

Example
~~~~~~~~~~~~~~~

.. code-block:: xml

  <Constitutive>
    <ExtendedDruckerPrager 
      name="edp"
      defaultDensity="2700"
      defaultBulkModulus="500"
      defaultShearModulus="300"
      defaultCohesion="0.0"
      defaultInitialFrictionAngle="15.0"
      defaultResidualFrictionAngle="23.0"
      defaultDilationRatio="1.0"
      defaultHardening="0.001"
    />
  </Constitutive>
