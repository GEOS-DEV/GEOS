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
* Discretizations
  * Finite Element Interface
  * Finite Volume Stencil
* Physics Solvers
  * Solid Mechanics Explcit on GPU, Implicit Assembly on GPU
  * Single Phase Flow
  * Composiitonal Multiphase Flow
  * Surface Generation
    * Topology change (Legacy GEOSX approach)
    * Embedded Enhanced Assumed Strain
  * Hydrofracture Solver
  * Small Strain aligned Contact using Lagrange Multipliers
* Mesh Structure
  * Introduced concept of extrinsic mesh data
  * Fracture Elements to represent FV cells in fracture
* Proppant modeling 

Version v0.1.0 -- Release date 2018-02-15
==========================================
Initial Code Release containing:
* Data Repository
  * Group
  * Wrapper
  * Input processing 
* Physics Solver hierarcy
  * Solid Mechanics
* Mesh data structure
  * NodeManager, EdgeManager, FaceManager, ElementManager
  * Silo Output
