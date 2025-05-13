.. _GeomechanicsModel:

############################################
Geomechanics Models
############################################

Damage Overview
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
	FSLOPE_{h} = \varsigma^{\psi} * FSLOPE_{i} + \left[(1-\varsigma^{\psi}) *  FSLOPE_{f}  \right]

Above, the hardened FSLOPE ( :math:`$FSLOPE_{h}$` ) is updated according to :math:`\varsigma` 
and :math:`\psi`, a fracture softening exponent.  We reduce $\text{YSLOPE}$ if necessary to 
ensure the $\text{YSLOPE}_h < \text{FSLOPE}_h < $.  The maximum hydrostatic tensile stress 
evolves with both damage and hardening as:

.. math::
	peak_{\bar{I}1 h} = \varsigma^{\psi} \left( peak_{\bar{I}1} + H/FSLOPE_h \right)

User Inputs and Code
========================
In the particle file writer (pfw) scripts, the user will define:
.. code-block:: xml

   fractureEnergyReleaseRate =1.5e-8
   fractureStress= 0.0184   #GPa
   fractureSofteningExponent = 0.75  
   damageEvolutionCriterion = 1      
   brittleDuctileTransition = 0.020

Where fracture energy release rate is 0 or negative to disable damage. 
The fracture softening exponent (:math:`\psi`) controls how quickly the strength drops to a
resisdual value as damage approaches 1. The damage evolution criterion is set to either 0, 1, or 2.
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






Plasticity and Hardening Overview
========================
For geomechanics deformation, we typically work from the “Arenisca” constitutive model that was developed originally to predict the effect of explosives used on drained
and saturated sandstones for wellbore completion technology (Homel et al. (`2015 <https://link.springer.com/article/10.1007/s00707-015-1407-2>`_, 2014 (internal manual)). The Arenisca constitutive 
model itself adapted functional forms for strength,compaction, and elasticity from the Kayenta material model based out of Sandia National Labs (Brannon et al., 2009 (internal manual)),
 which was adapted from an earlier Sandia geo-model.

The original yield surface defined by Arenisca is shown here (i) and is allowed to evolve to a hardened yield surface (h)


incorporate brittle failure by using a simple Rankine damage model following the approach 
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
	FSLOPE_{h} = \varsigma^{\psi} * FSLOPE_{i} + \left[(1-\varsigma^{\psi}) *  FSLOPE_{f}  \right]

Above, the hardened FSLOPE ( :math:`$FSLOPE_{h}$` ) is updated according to :math:`\varsigma` 
and :math:`\psi`, a fracture softening exponent.  We reduce $\text{YSLOPE}$ if necessary to 
ensure the $\text{YSLOPE}_h < \text{FSLOPE}_h < $.  The maximum hydrostatic tensile stress 
evolves with both damage and hardening as:

.. math::
	peak_{\bar{I}1 h} = \varsigma^{\psi} \left( peak_{\bar{I}1} + H/FSLOPE_h \right)

User Inputs and Code
========================
In the particle file writer (pfw) scripts, the user will define:
.. code-block:: xml

   fractureEnergyReleaseRate =1.5e-8
   fractureStress= 0.0184   #GPa
   fractureSofteningExponent = 0.75  
   damageEvolutionCriterion = 1      
   brittleDuctileTransition = 0.020

Where fracture energy release rate is 0 or negative to disable damage. 
The fracture softening exponent (:math:`\psi`) controls how quickly the strength drops to a
resisdual value as damage approaches 1. The damage evolution criterion is set to either 0, 1, or 2.
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
