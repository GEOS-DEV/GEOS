.. _ChiumentiModel:

############################################
Chiumenti Damage Model (MPM Only)
############################################

Overview, Summary
========================
The Chiumenti model is a brittle damage formulation designed for materials such as ceramics, 
where failure is governed by tensile cracking rather than ductile plasticity :cite:p:`DFG`. The model builds 
on a hyperelastic constitutive response, typically implemented as a compressible neo-Hookean 
trial stress, and introduces damage evolution based on the maximum principal stress. Damage is 
incremented when the maximum tensile stress exceeds a user-defined failure strength and reduces 
tensile principal stresses accordingly.

The model implements Rankine-type fracture behavior with linear strain-softening, where the 
softening rate is tied to the Mode I fracture energy and a characteristic length, providing 
mesh-independent results. Optional extensions include a time-to-failure model to predict fracture 
propagation speed and stochastic variations of particle strength via a Weibull distribution to 
represent material heterogeneity.

This implementation captures the essential features described by :cite:t:`CerveraAndChiumenti`, 
including tension-driven damage, scalar damage evolution, and fracture-energy regularization. 
While it does not explicitly implement the full effective stress framework, Kuhn–Tucker conditions, or 
the consistent tangent operator, the simplified, solver-friendly formulation reproduces the core 
Rankine damage behavior, allowing verification of fracture response and energy dissipation in finite element simulations.

Key features of the model include:

#. Hyperelastic stress response (elastic predictor)
#. Damage initiation based on maximum principal stress
#. Fracture-energy-controlled softening
#. Length-scale regularization for mesh independence
#. Selective degradation of tensile stresses

The model is particularly well-suited for brittle materials where tensile cracking dominates 
the failure mechanism.


Theory
========================

Elastic Response
------------------------
The model begins with a standard hyperelastic stress update, computed from the deformation gradient providing
a trial stress assuming purely elastic behavior:

(Equation 1)

.. math::

   \boldsymbol{\sigma}^{\text{trial}} = \text{HyperelasticUpdate}(F)


Inelastic Toggle
------------------------
If inelastic effects are disabled, the trial stress is returned directly without further updates.

Material Properties for Fracture-Energy Softening
------------------------
Evaluate element-wise material properties required for fracture-energy-based softening:

(Equation 2)

.. math::

   \sigma_f^{(k)} = \sigma_f \cdot \text{strengthScale}[k]

(Equation 3)

.. math::

   E = E(\lambda, G), \quad
   \ell = \ell_{\text{crit}} \cdot \text{lengthScale}[k]

(Equation 4)

.. math::

   \beta = \frac{\sigma_f^2}{2 E G_f}

(Equation 5)

.. math::

   \beta^* = \frac{\beta \ell}{1 - \beta \ell}

Principal Stress Decomposition
------------------------
Perform spectral decomposition of the stress tensor to compute principal stresses and directions:

(Equation 6)

.. math::

   \boldsymbol{\sigma} = \sum_{i=1}^{3} \sigma_i \, \mathbf{n}_i \otimes \mathbf{n}_i



Damage Criterion
------------------------
Determine if the maximum principal stress exceeds the scaled failure threshold:

(Equation 7)

.. math::

   \sigma_{\max} = \max(\sigma_1, \sigma_2, \sigma_3)

If exceeded, compute provisional new damage using the fracture-energy-based softening law:

(Equation 8)

.. math::

   d_{\text{new}} =
   \begin{cases}
   (1 + \beta^*) \left(1 - \frac{\sigma_f}{\sigma_{\max}}\right), & \sigma_{\max} \le \sigma_f \left(1 + \frac{1}{\beta^*}\right) \\
   1, & \text{otherwise}
   \end{cases}

Damage Update
------------------------
Constrain damage to be monotonically increasing and bounded:

(Equation 9)

.. math::

   d = \max(d_{\text{old}}, d_{\text{new}}), \quad 0 \le d \le 1

Stress Degradation
------------------------
Damage affects only tensile principal stresses:

(Equation 10)

.. math::

   \sigma_i^{\text{damaged}} =
   \begin{cases}
   (1 - d)\,\sigma_i, & \sigma_i > 0 \\
   \sigma_i, & \sigma_i \le 0
   \end{cases}

Reconstruct the full stress tensor from the modified principal stresses:

(Equation 11)

.. math::

   \boldsymbol{\sigma}^{\text{damaged}} =
   \sum_{i=1}^{3} \sigma_i^{\text{damaged}} \, \mathbf{n}_i \otimes \mathbf{n}_i


Implementation
========================

Chiumenti Hyperelastic + Damage Update - Function Overview
========================
This function updates the Cauchy stress using a hyperelastic predictor followed by a 
damage correction based on principal stresses.

.. list-table:: The 7-Step Process with Inline Variable Annotations
   :widths: 5 95
   :header-rows: 1

   * - Step
     - Description
   * - 1
     - **Compute elastic trial stress:** Call hyperelastic update using deformation gradient by computing :math:`\boldsymbol{\sigma}^{\text{trial}}` 
       from :math:`F` **(Equation 1)** .
       
   * - 2
     - **Check inelasticity flag:**  
       If inelasticity is **disabled** (:math:`m_{\text{disableInelasticity}} = \text{true}`), return trial stress :math:`\boldsymbol{\sigma}^{\text{trial}}` immediately.  
       If inelasticity is **enabled** (:math:`m_{\text{disableInelasticity}} = \text{false}`), proceed to evaluate damage using :math:`\sigma_{\max}` and :math:`\sigma_f` **(Equations 2,3)**.
       
   * - 3
     - **Compute material properties:** Apply element-wise strength scaling to obtain :math:`\sigma_f` and compute characteristic length :math:`\ell` for fracture-energy regularization **(Equations 5-7)**.  
       Compute Young’s modulus :math:`E`, brittleness factor :math:`\beta`, and scaled brittleness :math:`\beta^*`.
       
   * - 4
     - **Compute principal stresses:** Perform eigenvalue decomposition of the stress tensor and compute principal stresses :math:`\sigma_i` and principal directions :math:`\mathbf{n}_i` **(Equations 8,9)**.  
              
   * - 5
     - **Evaluate damage criterion:** Determine if :math:`\sigma_{\max}` exceeds the scaled failure threshold :math:`\sigma_f`.  
       Compute provisional new damage :math:`d_{\text{new}}` using fracture-energy-based softening **(Equation 10)**  
       which includes the scaled brittleness factor :math:`\beta^*`.
       
   * - 6
     - **Update damage variable:** Constrain the damage to be monotonically increasing and bounded:  
       :math:`d = \max(d_{\text{old}}, d_{\text{new}}), \quad 0 \le d \le 1` **(Equation 11)**.
       
   * - 7
     - **Apply damage to stress:** Reduce tensile principal stresses and reconstruct full stress tensor. 
       Compute :math:`\sigma_i^{\text{damaged}} = (1-d)\,\sigma_i` for :math:`\sigma_i>0` and reconstruct :math:`\boldsymbol{\sigma}^{\text{damaged}}` **(Equations 12)**  .

Parameters, User Inputs
========================

.. include:: /coreComponents/schema/docs/Chiumenti.rst


References
========================

.. bibliography::
   :style: plain
   :filter: docname in docnames