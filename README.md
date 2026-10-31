[![DOI](https://zenodo.org/badge/131810578.svg)](https://zenodo.org/badge/latestdoi/131810578)
[![codecov](https://codecov.io/github/GEOS-DEV/GEOS/graph/badge.svg?token=0VTEHPQG58)](https://codecov.io/github/GEOS-DEV/GEOS)
[![CI](https://github.com/GEOS-DEV/GEOS/actions/workflows/ci_tests.yml/badge.svg)](https://github.com/GEOS-DEV/GEOS/actions?query=branch%3Adevelop)
[![docs](https://readthedocs.com/projects/geosx-geosx/badge/?version=latest)](https://geosx-geosx.readthedocs-hosted.com/en/latest/)

Welcome to the GEOS project!
-----------------------------
GEOS is a simulation framework for modeling coupled flow, transport, and geomechanics
in the subsurface.  The code provides advanced solvers for a number of target applications,
including
  - carbon sequestration,
  - geothermal energy,
  - and similar systems.  

Its capabilities include single-phase, multiphase and compositional flow (isothermal or thermal), wells,
solid mechanics and poromechanics, fracture and fault mechanics (contact, embedded fractures, hydraulic fracturing),
wave propagation and induced seismicity.

A key focus of the project is achieving scalable performance on current and next-generation
high performance computing systems.  We do this through a portable programming model and research into scalable algorithms.

You may want to browse our
[publications](https://geosx-geosx.readthedocs-hosted.com/en/latest/docs/sphinx/Publications.html)
page for more details on the HPC, numerics,
and applied engineering components of this effort.

Getting started
---------------

GEOS is distributed as source code. To download, build and run a first simulation, follow the
[Quick Start Guide](https://geosx-geosx.readthedocs-hosted.com/en/latest/docs/sphinx/QuickStart.html),
then work through the [tutorials](https://geosx-geosx.readthedocs-hosted.com/en/latest/docs/sphinx/tutorials/Index.html).
Docker images with all third-party libraries pre-built are available (see the Quick Start Guide), and
this repository can be opened directly in a GitHub Codespace or a VS Code Dev Container.

Documentation
-------------

Please visit the [Main documentation for GEOS](https://geosx-geosx.readthedocs-hosted.com/en/latest/).

If you would like to contribute to GEOS, please see the [developer guide](https://geosx-geosx.readthedocs-hosted.com/en/latest/docs/sphinx/developerGuide/Index.html)

If you would like to report a bug, please submit an [issue](https://github.com/GEOS-DEV/GEOS/issues/new/choose). 

Who develops GEOS?
-------------------
GEOS is an open source project and is developed by a community of researchers at
several institutions.  The bulk of the code has been written by contributors from
four main organizations:
  - Lawrence Livermore National Laboratory,
  - Stanford University,
  - TotalEnergies,
  - Chevron

See our
[authors](https://geosx-geosx.readthedocs-hosted.com/en/latest/docs/sphinx/Contributors.html)
and
[acknowledgements](https://geosx-geosx.readthedocs-hosted.com/en/latest/docs/sphinx/Acknowledgements.html)
page for more details.  

History
-------
GEOS was known as GEOSX until 2023; some file names, library names and URLs still use the old name.
It succeeds an earlier code developed at LLNL, also called GEOS. The current code differs from that earlier effort in two important ways:
  - It uses a fundamentally different programming model to achieve
    high performance on the complicated chip architectures common on today's
    HPC systems, including GPU-based exascale-class systems.
  - It is developed as an open-source effort to encourage collaboration
    within the research and industrial community, under the [LGPL 2.1 License](./LICENSE).


Release
-------

For release details and restrictions, please read the [LICENSE](./LICENSE) file.

For copyrights, please read the [COPYRIGHT](./COPYRIGHT ) file.

For contributors, please read the [CONTRIBUTORS](./CONTRIBUTORS ) file.

For acknowledgements, please read the [ACKNOWLEDGEMENTS](./ACKNOWLEDGEMENTS ) file.

For notice, please read the [NOTICE](./NOTICE ) file.

`LLNL-CODE-812638`  `OCEC-18-021`
