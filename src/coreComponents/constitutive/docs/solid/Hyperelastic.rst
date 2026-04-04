.. _HyperelasticModel:

############################################
Hyperelastic Model (MPM Only)
############################################

Overview, Summary
========================
The hyperelastic constitutive model in this code describes the isotropic elastic response of a 
solid :cite:p:`fu2001nonlinear`
capable of handling both small and large deformations. It uses fundamental material properties such as
bulk modulus and shear modulus (or derived constants like Young’s modulus and Poisson ratio) to compute
stresses from strains. For small-strain problems, it applies linear elasticity using a 6×6 stiffness 
matrix derived from the Lamé constants, optionally including rotation effects. For finite strains, it
computes the deformation gradient and the right Cauchy-Green tensor, separates volumetric and 
deviatoric contributions, and calculates the stress accordingly, ensuring accurate response to both 
volume changes and distortions. The model also supports flexible initialization from different pairs of 
elastic constants and provides both stress-only and stress-plus-stiffness updates for use in explicit 
or implicit numerical solvers. The code ensures that exactly two elastic constants are provided by the user and then computes all other 
related elastic constants consistently, or throws an error if the input is invalid.

Theory 
========================

Baseline Stiffness
------------------------
The solver prepares the elastic stiffness matrix, encoding how the material responds linearly to small strains.
Later update steps will compute the stress with this stiffness and the strain for deformation gradient.
First the 6x6 stiffness matrix for linear isotropic elasticity is computed in Voigt notation at a given element, k. 
Material constants including the shear modulus (``G``) and the first Lamé constant (computed from bulk modulus ``K`` and 
shear modulus ``G``) are retrieved and the stiffness matrix is initialized to zero.
Hooke's Law for isotropic Materials is encoded by first setting 
the normal (diagonal) components for the 3 normal strains (xx, yy, zz):
     (Equation 1)

     .. math::

        C_{ii} = \lambda + 2G

And the off-diagonal terms:

     (Equation 2)

     .. math::

        C_{ij} = \lambda \quad \text{for } i \neq j, \; i,j = 1,2,3

Shear components (xy, yx, xz) on the diagonal are set accordinly:
   
   (Equation 3)

   .. math::

      C_{44} = C_{55} = C_{66} = G

.. raw:: html

   <br><br><br><br>

Elastic Strain
------------------------
Hooke’s law is inverted to compute the elastic strain components from a given stress tensor,
using the element’s material constants. Elastic strain at a element ``k`` is and quatrature ``q`` is computed from the stress state 
with **linear isotropic elasticity** in Voigt notation. Material constants ``E`` , Young's Modulus, and 
``nu``, Poisson's ratio, are computed from ``G`` and ``K``. Normal strains in the xx, yy, and zz directions
and shear strains in the xy, yz, and xz directions are computed from respective directions: 

  (Equations 4-9)

   .. math::
     \begin{aligned}
     \varepsilon_{xx} &= \frac{\sigma_{xx} - \nu \, \sigma_{yy} - \nu \, \sigma_{zz}}{E} \\
     \varepsilon_{yy} &= \frac{-\nu \, \sigma_{xx} + \sigma_{yy} - \nu \, \sigma_{zz}}{E} \\
     \varepsilon_{zz} &= \frac{-\nu \, \sigma_{xx} - \nu \, \sigma_{yy} + \sigma_{zz}}{E} \\
     \varepsilon_{xy} &= \frac{\sigma_{xy}}{G} \\
     \varepsilon_{yz} &= \frac{\sigma_{yz}}{G} \\
     \varepsilon_{xz} &= \frac{\sigma_{xz}}{G} 
     \end{aligned}

.. raw:: html

   <br><br><br><br>

Small-Strain Updates
------------------------
The small-strain functions provide the constitutive model with lightweight “no-state” or incremental update 
pathways. These are designed to compute stress and/or stiffness without necessarily updating internal state 
arrays, which can be useful for explicit solvers or solver routines that need stiffness information separately.

1. **smallStrainNoStateUpdate_StressOnly**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This is a placeholder function. It takes in total strain and a stress array but does not modify them. 
Its purpose is to satisfy the solver interface when a stress-only update is required but no computation is needed.

``smallStrainNoStateUpdate_StressOnly(k, q, totalStrain, stress)``



2. smallStrainNoStateUpdate
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
There are two overloaded versions. In both versions, stress is not updated in the constitutive model; 
these functions exist purely to expose stiffness information

**Full 6×6 stiffness version**: Computes the elastic stiffness matrix for the element, using getElasticStiffness.
 Stress is passed through but not modified (no state update).

 ``smallStrainNoStateUpdate(k, q, totalStrain, stress, stiffness[6][6])``


**DiscretizationOps version**: Instead of a full stiffness matrix, it writes the element’s bulk and shear 
moduli into a DiscretizationOps object. This allows the solver to access stiffness in a more abstract, optimized way.

``smallStrainNoStateUpdate(k, q, totalStrain, stress, stiffness)``



3. **smallStrainUpdate_StressOnly**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
These functions compute stress from a strain increment:

**Non-rotation version**: Computes incremental stress, adds it to the previous stress, and saves the result as 
the new stress.

``smallStrainUpdate_StressOnly(k, q, timeIncrement, strainIncrement, stress)``


**Rotation-aware version**: Accepts rigid-body rotation tensors, but currently ignores them. Internally, 
it calls the non-rotation version.

``smallStrainUpdate_StressOnly(k, q, timeIncrement, beginningRotation, endRotation, strainIncrement, stress)``



4. **smallStrainUpdate**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
These are overloaded to either return a full stiffness matrix or a DiscretizationOps object along with updated stress.
In both cases, stress is updated using the incremental strain, and the corresponding stiffness 
(matrix or moduli) is returned for solver use.

**Full 6×6 stiffness version**:

``smallStrainUpdate(k, q, timeIncrement, strainIncrement, stress, stiffness[6][6])``

**DiscretizationOps version**:

``smallStrainUpdate(k, q, timeIncrement, strainIncrement, stress, stiffness)``

.. raw:: html

   <br><br><br><br>


Hyperelastic Stress
------------------------
A full finite-strain hyperelastic stress is updated using the deformation gradient, 
computing both the volumetric and deviatoric contributions, and stores the result in Voigt notation.
The Cauchy stress for a hyperelastic (Neo-Hookean-like) material is calculated given a deformation 
gradient increment F - I. It is part of a finite-strain update, so unlike the small-strain functions, 
it works directly with the full deformation rather than incremental linearized strains.

The total deformation gradient, F, is computed from the input increment F-I whereI is the identity matrix 
added back to the establish the absolute deformation gradient.

(Equation 10)

.. math::

   F = (F - I) + I

The Cauchy-Green deformation tensor, C, is computed and implemented via a double loop over the matrix entries, and it's trace
is determined:

(Equations 11-13)

.. math::

   C = F^T F

.. math::
   C_{ij} = \sum_{k=0}^{2} F_{ik} F_{jk}

.. math::

   \text{tr}(C) = C_{00} + C_{11} + C_{22}


Volume change is determined through calculating the determinant of F, the Jacobian, J:

(Equation 14)

.. math::

   J = \det(F) = F_{00}(F_{11} F_{22} - F_{12} F_{21}) - F_{01}(F_{10} F_{22} - F_{12} F_{20}) + F_{02}(F_{10} F_{21} - F_{11} F_{20})

The 6-component stress vector is assembled in Voigt notiation using computed stress scaling factros for the 
Neo-Hookean formulation.

(Equations 15-22)

.. math::

   x_1 = \frac{G}{J^{5/3}}, \quad
   x_2 = K (J - 1) - \frac{G}{3 J^{5/3}} \text{tr}(C)


.. math::
   \begin{aligned}
   \sigma_{xx} &= x_1 C_{00} + x_2, \\
   \sigma_{yy} &= x_1 C_{11} + x_2, \\
   \sigma_{zz} &= x_1 C_{22} + x_2, \\
   \sigma_{yz} &= x_1 C_{12}, \\
   \sigma_{xz} &= x_1 C_{02}, \\
   \sigma_{xy} &= x_1 C_{01}
   \end{aligned}



The **HyperelasticUpdates::hyperUpdate** function is overloaded: one version computes only the 
hyperelastic stress (Equations 10–22), while the other version takes an additional stiffness argument and 
computes both the stress and the elastic stiffness matrix (Equations 1–3 and 10–22), 
providing the solver with both the current stress state and the tangent stiffness when needed.
The Hyperelastic class wraps the HyperelasticUpdates kernel, providing element-wise material properties and
factory methods to instantiate stress and stiffness update kernels for use by the solver.



Parameters, User Inputs
========================

.. include:: /coreComponents/schema/docs/Hyperelastic.rst




References
========================

.. bibliography::
   :style: plain
   :filter: docname in docnames
