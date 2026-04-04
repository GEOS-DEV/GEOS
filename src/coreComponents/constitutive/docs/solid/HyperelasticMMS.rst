.. _HyperelasticMMSModel:

############################################
HyperelasticMMS Model (MPM Only)
############################################

Overview, Summary
========================
The HyperelasticMMS model differs from the standard Hyperelastic model primarily in its 
purpose and implementation for verification: it is designed specifically for the Method of 
Manufactured Solutions (MMS), providing “update” kernels (HyperelasticMMSUpdates) that allow stress and 
stiffness computations to be tested against analytically prescribed solutions, rather than modeling real material 
response. Structurally, it inherits the same small-strain update and hyperelastic stress routines as 
the standard hyperelastic model, but it exposes additional flexibility for stress-only or no-state updates, 
supports explicit specification of Lame constants and shear modulus per element, and includes infrastructure 
to separate kernel updates from the main state data. Unlike the standard hyperelastic model, its primary focus 
is on verification and testing, rather than predictive simulation. Key differences between the two models are
summarized:



+--------------------+---------------------------------------+--------------------------------------------+
| Feature            | Hyperelastic                          | HyperelasticMMS                            |
+====================+=======================================+============================================+
| Purpose            | Predictive simulation of isotropic    | Verification using Method of Manufactured  |
|                    | hyperelastic materials                | Solutions (MMS)                            |
+--------------------+---------------------------------------+--------------------------------------------+
| Kernel Updates     | Standard small-strain and hyperelastic| Small-strain and hyperelastic stress       |
|                    | stress updates                        | updates designed for MMS; supports stress- |
|                    |                                       | only or no-state updates                   |
+--------------------+---------------------------------------+--------------------------------------------+
| State Management   | Maintains full stress and stiffness   | Can optionally bypass full state (no-state |
|                    | state                                 | updates) for testing                       |
+--------------------+---------------------------------------+--------------------------------------------+
| Input Parameters   | Typically derived from bulk/shear     | Direct per-element Lame constants          |
|                    | modulus, Young’s modulus, or Poisson  | (`lambda`)and shear modulus (`G`); supports|
|                    | ratio                                 | default values for MMS                     |
+--------------------+---------------------------------------+--------------------------------------------+
| Stiffness/Stress   | Integrated with material state        | Separate `HyperelasticMMSUpdates` kernels; |
| Computation        |                                       | flexible for verification                  |
+--------------------+---------------------------------------+--------------------------------------------+
| Focus              | Accurate material response under      | Analytical solution verification; stress   |
|                    | small or finite strains               | and stiffness comparisons for testing      |
+--------------------+---------------------------------------+--------------------------------------------+
| Post-Processing    | Standard update routines              | Includes post-input initialization and     |
|                    |                                       | flexible kernel creation for test cases    |
+--------------------+---------------------------------------+--------------------------------------------+

Theory
======

The model assumes **small-strain isotropic hyperelasticity** with per-element **first Lame constant** (`lambda`) and **shear modulus** (`G`). The stress and stiffness follow standard linear elasticity relations:

1. **Elastic stiffness tensor** `C` (Voigt notation):

   .. math::

      C_{ijkl} =
      \lambda \delta_{ij} \delta_{kl} + G (\delta_{ik} \delta_{jl} + \delta_{il} \delta_{jk})

   In Voigt 6x6 form:

   .. code-block:: text

      stiffness[0][0] = lambda + 2*G
      stiffness[0][1] = lambda
      stiffness[0][2] = lambda
      stiffness[1][0] = lambda
      stiffness[1][1] = lambda + 2*G
      stiffness[1][2] = lambda
      stiffness[2][0] = lambda
      stiffness[2][1] = lambda
      stiffness[2][2] = lambda + 2*G
      stiffness[3][3] = G
      stiffness[4][4] = G
      stiffness[5][5] = G

2. **Elastic strain** from stress using Young's modulus `E` and Poisson ratio `nu`:

   .. math::

      \epsilon_{xx} = \frac{\sigma_{xx} - \nu (\sigma_{yy} + \sigma_{zz})}{E}, \\
      \epsilon_{yy} = \frac{\sigma_{yy} - \nu (\sigma_{xx} + \sigma_{zz})}{E}, \\
      \epsilon_{zz} = \frac{\sigma_{zz} - \nu (\sigma_{xx} + \sigma_{yy})}{E}, \\
      \epsilon_{yz} = \frac{\sigma_{yz}}{G}, \quad
      \epsilon_{xz} = \frac{\sigma_{xz}}{G}, \quad
      \epsilon_{xy} = \frac{\sigma_{xy}}{G}

3. **Hyperelastic stress update** (Green-Lagrange formulation for MMS):

   - Compute deformation gradient `F` from `F - I`:

     .. code-block:: text

        F[i][j] = FminusI[i][j]
        F[i][i] += 1  # Add identity

   - Right Cauchy-Green tensor:

     .. math::

        C = F^T F

   - Jacobian determinant:

     .. math::

        J = det(F)

   - Stress:

     .. math::

        \sigma_{ii} = \frac{\lambda \ln J}{J} + \frac{G}{J} (C_{ii} - 1), \quad
        \sigma_{ij} = \frac{G}{J} C_{ij} \quad (i \neq j)

Implementation
==============

The implementation separates the **data container** from the **update kernels**:

- **HyperelasticMMS**: main constitutive class  
  - Holds per-element constants: `lambda`, `G`  
  - Supports default values for Lame constants  
  - Provides factory methods for kernel updates  

- **HyperelasticMMSUpdates**: kernel class  
  - Performs stress-only or full stress + stiffness updates  
  - Provides `smallStrainUpdate` and `hyperUpdate` methods  
  - Can bypass material state arrays for MMS verification  

**Example usage**:

.. code-block:: cpp

    HyperelasticMMS material("TestMMS", parentGroup);

    // Create kernel updates including state
    auto updates = material.createKernelUpdates(true);

    // Access per-element Lame constants
    arrayView1d<real64 const> lambda = material.getLambda();
    arrayView1d<real64 const> G = material.getShearModulus();

    // Update stress at element k, quadrature q
    real64 stress[6];
    real64 FminusI[3][3] = { ... }; // deformation gradient minus identity
    updates.hyperUpdate(k, q, FminusI, stress);

Notes
=====

- HyperelasticMMS is intended for **verification and testing** rather than production predictive simulations.  
- Stress-only or no-state updates allow comparison against analytical MMS solutions without maintaining full history arrays.  
- All stiffness and stress computations are **isotropic and linear in small-strain regime**, consistent with MMS requirements.  

Parameters, User Inputs
========================

.. include:: /coreComponents/schema/docs/Graphite.rst




References
========================

.. bibliography::
   :style: plain
   :filter: docname in docnames
