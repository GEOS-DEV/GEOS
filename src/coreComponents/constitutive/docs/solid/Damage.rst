.. _DamageModel:

############################################
Damage Model (MPM Only)
############################################

Overview
========================
We incorporate brittle failure by using a simple Rankine damage model following the approach 
of `Homel et al, 2015 <https://link.springer.com/article/10.1007/s00466-015-1187-5>`_.  Brittle fracture was observed experimentally in triaxial compression tests 
at low confining pressure.  To account for this, we define a fracture stress (:math:`\sigma_{fracture}`) 
under triaxial compression. The fracture model will increment damage if the fracture shear stress 
is exceeded but only if the pressure is below the specified brittle-ductile transition pressure 
(:math:`\sigma_{BD}`), which we define to correspond with the apex of the initial yield surface. 
The dissipation associated with fracturing in a continuum damage model will be inherently 
length-scale dependent unless a regularization is applied.  In our model, the increment in damage, 
D, is calculated based on the increment in plastic work relative to the 
fracture energy release rate, :math:`G_{f}`, and is normalized by a length scale, L.

.. math::
	\Delta W_p = \sigma : \dot{\epsilon}_p \Delta t

.. math::
	\Delta D = \frac{\Delta W_p L}{G_f}

Finally, we compute an updated value of coherence (:math:`\varsigma=1-D`), 
which is used to scale the yield surface parameters:

.. math::
	\varsigma = \varsigma_{n} - \Delta D.

The yield surface softening has a nonlinear dependence on the damage variable allowing us to 
control the rate of softening which has a competing effect with strain hardening: 

.. math::

   \mathrm{FSLOPE}_{h}
   = \varsigma^{\psi} \times \mathrm{FSLOPE}_{i}
   + \left[ \left(1-\varsigma^{\psi}\right) \times \mathrm{FSLOPE}_{f} \right]



Above, the hardened FSLOPE (:math:`\text{FSLOPE}_h`) is updated according to :math:`\varsigma` 
and :math:`\psi`, a fracture softening exponent. We reduce :math:`\text{YSLOPE}_h` if necessary 
to ensure :math:`\text{YSLOPE}_h < \text{FSLOPE}_h`. The maximum hydrostatic tensile stress 
evolves with both damage and hardening (H) as:

.. math::
   \text{peak}_{\bar{I}_1,h} =
   \varsigma^{\psi} \left( \text{peak}_{\bar{I}_1} + \frac{H}{\text{FSLOPE}_h} \right)

where H is determined as:

.. math::
   H = K_h \left( 1 - e^{\,n_h \gamma^p} \right)



As damage approaches 1, the material loses all cohesive strength
(:math:`\mathrm{peak}_{\bar{I}_1 h} \to 0`), and the initial slope of
the yield surface reduces to :math:`FSLOPE_f`. 
:math:`X` and therefore the cap function :math:`F_c`,
evolve with volumetric plastic strain :math:`\epsilon_v^p` where :math:`\epsilon_v^{p,n}` 
relates the initial and current porosity.


.. math::
   X(\epsilon_v^p) =
   \begin{cases}
     \dfrac{\ln\left( \dfrac{\epsilon_v^p + p_3}{p_3} \right) + p_0 p_1}{p_1},
       & \epsilon_v^p \le 0,\\[6pt]
     p_0 (\epsilon_v^p + 1)^{\dfrac{1}{p_0 p_1 p_3}},
       & \epsilon_v^p > 0
   \end{cases}



.. math::
   \epsilon_v^{p,n}
   = \ln\!\left(\dfrac{\phi_p - 1}{\phi_0 - 1}\right)




Parameters, User Inputs
========================

.. .. include:: /coreComponents/schema/docs/DamageParameterTable.rst
.. include:: /coreComponents/schema/docs/Damage.rst


Where fracture energy release rate is 0 or negative to disable damage. 
The fracture softening exponent (:math:`\psi`) controls how quickly the strength drops to a
residual value as damage approaches 1. The damage evolution criterion is set to either 0, 1, or 2.
If set to 0, This approach is designed to increment damage when there is plastic loading with dilation above some
threshold stress.  


The damage evolution is based on dilational plastic work.
If set to 1 or 2, the increment damage is based on the increment in plastic dilatational work relative to 
the fracture energy release rate, normalized by the length scale, to be per-unit-volume. 
If pressure (p=-I1/3) is below the brittle ductile transition pressure, increment damage based 
on plastic work increment relative to fracture energy regularized by length scale.
If damageEvolutionCriterion is set to 1, then the brittleductiletransition value (:math:`I_{DB}`) is a user-specified 
pressure. 
If damageEvolutionCriterion is set to 2, then :math:`I_{DB}` is calculated based upon the yield surface apex.
