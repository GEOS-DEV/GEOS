.. __StrainHardening Model:

############################################
Strain Hardening Model
############################################

Overview
========================
This damage model is intended for use with damage-field partitioning (DFG) within the
MPM solver, but can also be used without DFG by any solver. It is only appropriate for
schemes implementing explicit time integration.




Model Formulation
========================
First, we store the trial stress for computing the plastic strain increment. 

.. code-block:: cpp

  real64 trialStress[6] = { 0 };
  LvArray::tensorOps::copy< 6 >( trialStress, stress);


This essentially creates:

.. math::

   \sigma_i^{\text{trial}} = \sigma_i, \quad i = 1,\dots,6



Then we decompose the trial stress into pressure (P) and 
Von Mises (Q) stress invarients.

.. code-block:: cpp

  real64 trialP;
  real64 trialQ;
  real64 deviator[6] = { 0 };
  twoInvariant::stressDecomposition( trialStress,
                                     trialP,
                                     trialQ,
                                     deviator );


twoInvariant::stressDecomposition is defined in InvariantDecompositions.hpp. The function decomposes a stress tensor (Voigt notation) into volumetric stress, 
into deviatoric stress, and the deviator (or the unit orientation of the deviatoric part). To explore that file, you'll want to navigate to it through:

.. code-block:: cpp

  GEOS/src/coreComponents/constitutive/solid/InvariantDecompositions.hpp


The following section computes the unrotated state needed to evaluate the constitutive model
in the material frame. This involves extracting the rotation, unrotating the plastic strain,
and computing the right stretch tensor from the deformation gradient. 


.. code-block:: cpp

   real64 rotationTranspose[3][3];
   LvArray::tensorOps::transpose< 3, 3 >( rotationTranspose, beginningRotation );


The steps below remove the rigid-body rotation 
from both the deformation gradient and the plastic strain tensor so the 
constitutive model operates in the material (unrotated) configuration. 
By computing the unrotated plastic strain, unrotated deformation gradient,
 and finally the principal stretches and directions from the right stretch tensor, 
the model isolates the true material deformation needed 
for subsequent yield, hardening, and damage computations.


**Step 1 — Compute the transpose of the rotation tensor**

The tensor ``beginningRotation`` stores the rotation :math:`R` from the polar decomposition
of the deformation gradient. Its transpose :math:`R^{\mathsf{T}}` is computed and stored in
``rotationTranspose``. This is required because the constitutive model operates on **unrotated** quantities.

.. math::
   R^{\mathsf{T}} = R^{-1}


.. code-block:: cpp

   real64 oldPlasticStrain[6] = { 0 };
   LvArray::tensorOps::copy< 6 >( oldPlasticStrain, m_plasticStrain[k][q] );
   oldPlasticStrain[3] *= 0.5;
   oldPlasticStrain[4] *= 0.5;
   oldPlasticStrain[5] *= 0.5;

**Step 2 — Extract the plastic strain tensor and convert to tensor form**

``m_plasticStrain[k][q]`` is stored in engineering-strain Voigt notation.  
The shear components (indices 3–5) must be halved to convert from engineering shear
to tensorial shear strain. This produces ``oldPlasticStrain`` in true tensor form.

.. math::
   \epsilon_{12} = \tfrac{1}{2}\gamma_{12}



.. code-block:: cpp

   real64 unrotatedOldPlasticStrain[6] = { 0 };
   LvArray::tensorOps::Rij_eq_AikSymBklAjl< 3 >( unrotatedOldPlasticStrain,
                                                 rotationTranspose,
                                                 oldPlasticStrain );


**Step 3 — Unrotate the plastic strain**

The plastic strain tensor is unrotated using the following. 
This produces the plastic strain in the material (lattice) configuration.


.. math::
   \epsilon^{p,\text{unrot}} = R^{\mathsf{T}}\, \epsilon^p \, R


.. code-block:: cpp

   unrotatedOldPlasticStrain[3] *= 2.0;
   unrotatedOldPlasticStrain[4] *= 2.0;
   unrotatedOldPlasticStrain[5] *= 2.0;

**Step 4 — Return to engineering shear strain**

After unrotation, the tensorial shear strains are converted **back** to engineering form:

.. math::
   \gamma_{12} = 2\,\epsilon_{12}


.. code-block:: cpp

   real64 unrotatedDeformationGradient[3][3] = { { 0 } };
   LvArray::tensorOps::Rij_eq_AikBkj< 3, 3, 3 >( unrotatedDeformationGradient,
                                                 rotationTranspose,
                                                 m_deformationGradient[k] );


**Step 5 — Compute the unrotated deformation gradient**

The deformation gradient is unrotated using the following. This isolates the **right stretch tensor** contained in :math:`F^{\text{unrot}}`.


.. math::
   F^{\text{unrot}} = R^{\mathsf{T}} F



.. code-block:: cpp

   real64 U[6] = { 0.0 };
   LvArray::tensorOps::denseToSymmetric< 3 >( U, unrotatedDeformationGradient );


**Step 6 — Extract the symmetric right stretch tensor**

The symmetric tensor :math:`U` is extracted from the unrotated deformation gradient,
representing the **right stretch** from the polar decomposition:

.. math::
   F = R\,U


.. code-block:: cpp

   real64 stretch[3] = { 0 };
   real64 eigenVectors[3][3] = { { 0 } };
   LvArray::tensorOps::symEigenvectors< 3 >( stretch, eigenVectors, U );


**Step 7 — Compute principal stretches and directions**

Eigen-decomposition of :math:`U` gives:

- :math:`\lambda_i` → principal stretches (in ``stretch``)  
- corresponding eigenvectors → principal directions (in ``eigenVectors``)

This completes the calculation of the unrotated, principal-stretch representation of the deformation.
Then we find the largest eigenvalues and compare them to the maximum allowable 
failure stretch, which is temperature dependent via the ``thermalSoftening`` function

.. code-block:: cpp

   real64 failureStretch = m_maximumStretch * StrainHardeningPolymerUpdates::thermalSoftening(m_temperature[k], m_maximumStretchT0, m_maximumStretchA, m_maximumStretchB );     



.. code-block:: cpp

   real64 StrainHardeningPolymerUpdates::thermalSoftening( const real64 & T,
                                                           const real64 & T0,
                                                           const real64 & A,
                                                           const real64 & B
   ) const 
   { 
     if (std::abs(A) > 1.e-16)
     {
       return 1. + A / (1. + std::exp( B * (T-T0) ) );
     }
     else
     {
       return 1.;
     } 
   }
















The yield surface evolves with the yield strength, strain hardening slope, and 
shear softening magnitude all softening with temperature (in the thermalSoftening function).


.. code-block:: cpp

   yieldStrength = m_initialYield + plasticSoftening + stretchHardening;


.. math::

   \sigma_{\mathrm{y}} \;=\; \sigma_{0} \;+\; r_{0} \;+\; G_{r}


Then we incorporate the plasticity integration. We will soon compute the trial state, check the yield, 
return to the yield surface, and update the plastic strain. To achieve this, 
we need a workspace to hold the current unrotated plastic strain and a workspace 
to accumulate how much strain is added during the step.

.. code-block:: cpp
    real64 unrotatedTempPlasticStrain[6] = { 0 };
    real64 plasticStrainIncrement[6] = { 0 };



Next, we will do a  fixed-point iteration to find plastic strain and consistent 
return to updated yield surface. Iteratively, we:


**Step 8 – Relaxtion**

This code computes the magnitude of the plastic strain tensor, :math:`\gamma_p`, 
from the unrotated plastic strain in Voigt notation. The factor 0.5*(1 + (i < 3)) 
weights normal components (indices 0–2) by 0.5 and shear components (3–5) by 1.0 
so that, after summing the squared components,
:math:`\sqrt{\gamma_p}` gives the correct equivalent plastic shear strain.

.. code-block:: cpp

   real64 gamma_p = 0.0;
   for( int i = 0; i < 6; i++ )
   {
     gamma_p += 0.5*( 1 + (i < 3) ) * unrotatedTempPlasticStrain[i] * unrotatedTempPlasticStrain[i];
   }
   gamma_p = sqrt( gamma_p );



We start with :math:`r_0`, the initial relaxation with point of with plastic  strain to give 
plastic softening. Below, S is the shearSfoteningMagnitude 

.. math::

   r_0 = S \, \exp\!\left( \max\!\left( -\,\left( \frac{\gamma_p}{r_1} \right)^{r_2},\; -16 \right) \right)




**Step 1 — Damage check based on maximum principal stretch**

The array ``stretch`` contains the three principal stretches
:math:`\lambda_i`. The loop computes the maximum principal stretch

.. math::

   \lambda_{\max} = \max_{i=1,2,3} \lambda_i.

If this maximum stretch exceeds the prescribed failure stretch
:math:`\lambda_{\mathrm{fail}} = \texttt{failureStretch}`, the material
point is marked as fully damaged:

.. math::

   \texttt{m\_damage}[k][q] = 1
   \quad\text{if}\quad
   \lambda_{\max} > \lambda_{\mathrm{fail}}.

This corresponds to complete loss of load-carrying capacity at that
integration point.

**Step 2 — Parameters for the return-to-yield iteration**

The variables ``tol`` and ``maxEvals`` set the numerical controls for
the iterative return-to-yield algorithm used later in the update:

.. math::

   \texttt{tol} = 10^{-10}, \qquad \texttt{maxEvals} = 100.

A fixed-point iteration is currently used, but the comments note that a
Newton solver may be more efficient and robust.

**Step 3 — Temperature-dependent yield-strength components**

The previous-step yield strength is stored as

.. math::

   \sigma_{y}^{\text{old}} = \texttt{m\_yieldStrength}[k].

The updated yield-strength components are decomposed into three parts:
a base yield strength, a strain-hardening term, and a shear-softening
term, each scaled by the thermal-softening function
:math:`s_\theta(T; T_0, A, B)`:

.. math::

   \sigma_{y,0}
   &= \texttt{m\_defaultYieldStrength} \;
      s_\theta\bigl(T_k; T_{0,y}, A_y, B_y\bigr), \\[4pt]
   G_r
   &= \texttt{m\_strainHardeningSlope} \;
      s_\theta\bigl(T_k; T_{0,Gr}, A_{Gr}, B_{Gr}\bigr), \\[4pt]
   r_0
   &= \texttt{m\_shearSofteningMagnitude} \;
      s_\theta\bigl(T_k; T_{0,r}, A_r, B_r\bigr).

Here :math:`T_k = \texttt{m\_temperature}[k]` is the material-point
temperature, and each parameter triple :math:`(T_0, A, B)` defines the
temperature dependence of one contribution via the thermal-softening
law described in Eq. above.

These quantities (``yield0``, ``strainHardeningSlope``, and
``shearSofteningMagnitude``) will later be combined to form the updated
yield strength used in the plastic return.

**Step 4 — Allocate plastic-strain work arrays**

Finally, two temporary arrays are initialized:

- ``unrotatedTempPlasticStrain[6]`` — to store the unrotated trial
  plastic strain for this update step,
- ``plasticStrainIncrement[6]`` — to store the plastic strain increment
  computed by the return-mapping algorithm.

Both are initialized to zero in preparation for the subsequent
integration of the constitutive response.











.. Damage
========================
Damage is incorporated through the Damage constitutive model package. 
Damage documentation can be found here:   - :ref:`Damage <DamageModel>`



.. Parameters, User Inputs
========================
.. include:: /coreComponents/schema/docs/GeomechanicsParameterTable.rst
References
========================
.. bibliography::
   :cited:
   :style: plain
