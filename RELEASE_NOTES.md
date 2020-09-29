[comment]: # (-----------------------------------------------------------------)
[comment]: # (SPDX-License-Identifier: LGPL-2.1-only)
[comment]: # 
[comment]: # (Copyright 2018-2020 Lawrence Livermore National Security LLC)
[comment]: # (Copyright 2018-2020 The Board of Trustees of the Leland Stanford)
[comment]: # (                    Junior University)
[comment]: # (Copyright 2018-2020 Total, S.A)
[comment]: # (Copyright 2019-     GEOSX Contributors)
[comment]: # (All right reserved)
[comment]: # 
[comment]: # (For more details see:)
[comment]: # (  https://github.com/GEOSX/GEOSX/LICENSE)
[comment]: # (  https://github.com/GEOSX/GEOSX/COPYRIGHT)
[comment]: # (  https://github.com/GEOSX/GEOSX/CONTRIBUTORS)
[comment]: # (  https://github.com/GEOSX/GEOSX/NOTICE)
[comment]: # (  https://github.com/GEOSX/GEOSX/ACKNOWLEDGEMENTS)
[comment]: # (  https://github.com/GEOSX/GEOSX/RELEASE)


Version v0.2.0 -- Release date 2020-06-20
==========================================
* Extensive updates to Data Repository.
* Discretization:
  * Finite Element Interface
    * Kernel launching interface looping abstraction.
    * Element Formulations for 1st order:
      * 8 node hexahedron
      * 6 node wedge
      * 5 node pyramid
      * 4 node tetrahedron
  * Finite Volume Stencils for Classical (cell-centered) FVM formulation (TPFA, infrastructure support for MPFA, but currently 
    not implemented)
  * Inner Products for Hybrid FVM formulation (current support for TPFA and quasi-TPFA) 
* Physics Solvers
  * Solid Mechanics Explicit on GPU, Implicit Assembly on GPU
  * Single Phase Flow (Assembly on GPU)
    * Classical FVM and Hybrid FVM formulations
    * Porous matrix and DFM fracture flow
  * Compositional Multiphase Flow (Assembly on GPU)
  * Multi-segmented Wells for Single Phase and Compositional Multiphase Flow (Assembly on GPU)
    * Fully Implicit Isothermal Overall Composition Formulation
    * A library of fluid constitutive models:
      * Equation-of-State based hydrocarbon compositional
      * Three-phase extended Black-oil
      * Two-phase CO2-brine  
  * Surface Generation
    * Topology change (Legacy GEOS approach)
* Embedded Discrete Fractures
  * Enriched Finite Element Method for the discretization of the mechanics
  * Piecewise constant displacement jump enrichment
  * Hydrofracture Solver (Legacy GEOS approach)
  * Small Strain aligned Contact using Lagrange Multipliers
    * Discrete Fracture Model using a low order stabilized mixed finite element
  * Proppant Transport Solver
    * FVM formulation
    * Modeling the following major physical processes:
      * Proppant-fluid slurry flow and multicomponent transport in fractures
      * Proppant gravitational settling
      * Proppant bed build-up and development
      space
* Mesh Structure
  * Introduced concept of extrinsic mesh data
  * Fracture Elements to represent FV cells in fracture
* VTK output
* Linear Algebra Interface layers for Hypre, Trilinos, Petsc
  * Common interface for supported linear algebra packages
  * Krylov solvers (CG, GMRES, BiCGSTAB)
  * Preconditioners (Algebraic Multigrid, incomplete factorizations)
  * Block matrix and vector support
  * Serial and parallel direct solvers


Version v0.1.0 -- Release date 2018-02-15
==========================================
Initial Code Release containing:
* Data Repository
  * Group
  * Wrapper
  * Input processing 
* Physics Solver hierarchy
  * Solid Mechanics
* Mesh data structure
  * NodeManager, EdgeManager, FaceManager, ElementManager
  * Silo Output
