.. _GraphiteModel:

############################################
Graphite Model
############################################

Overview
========================
This is a damage model for modeling material failure in brittle materials. This damage model is intended for use with damage-field partitioning (DFG) within the
MPM solver, but can also be used without DFG by any solver. It is only appropriate for
schemes implementing explicit time integration. The model is really a hybrid plasticity/
damage model in the sense that we assume damaged material behaves like granular material
and hence follows a modified Mohr-Coulomb law. 





The modifications are that at low pressures,
the shape of the yield surface is modified to resemble a maximum principal stress criterion,
and at high pressures the shape converges on the von Mises yield surface. The damage
parameter results in softening of the deviatoric response i.e. causes the yield surface to
contract. Furthermore, damage is used to scale back tensile pressure: p = (1 - d) * pTrial.
pTrial is calculatd as pTrial = -k * log(J), where the Jacobian J of the material motion is
integrated and tracked by this model.






Parameters, User Inputs
========================

.. include:: /coreComponents/schema/docs/DamageParameterTable.rst


