Welcome to the GEOSX project!
-----------------------------
GEOSX is a simulation framework for modeling coupled flow, transport, and geomechanics
in the subsurface.  The code provides advanced solvers for a number of target applications,
including
  - carbon sequestration,
  - geothermal energy,
  - unconventional oil and gas,
  - and similar systems.  

A key focus of the project is achieving scalable performance on current and next-generation
high performance computing systems.  We do this through a portable programming model and research into scalable algorithms.

You may want to browse our
[publications](https://geosx-geosx.readthedocs-hosted.com/en/latest/docs/sphinx/Publications.html)
page for more details on the HPC, numerics,
and applied engineering components of this effort.

Who develops GEOSX?
-------------------
GEOSX is an open source project and is developed by a community of researchers at
several institutions.  The bulk of the code has been written by contributors from
three main organizations:
  - Lawrence Livermore National Laboratory
  - Stanford University
  - Total, S.A.

See our
[authors](https://geosx-geosx.readthedocs-hosted.com/en/latest/docs/sphinx/Contributors.html)
and
[acknowledgements](https://geosx-geosx.readthedocs-hosted.com/en/latest/docs/sphinx/Acknowledgements.html)
page for more details.  

How does GEOSX relate to GEOS?
------------------------------
GEOSX is the offshoot of an earlier code developed at LLNL called GEOS.  The new
code differs from our previous efforts in two important ways:
  - GEOSX uses a fundamentally different programming model to achieve
    high performance on the complicated chip architectures common on today's
    HPC systems.  The "X" in GEOSX emphasizes our goal of being ready for exascale-class systems as they are delivered.
  - The new code has been released as an open-source effort to encourage collaboration
    within the research and industrial community.  See the release notes below
    for details of the [LGPL 2.1 License](./LICENSE) that has been adopted.

Documentation
---------------------

Documentation is hosted [here](https://geosx-geosx.readthedocs-hosted.com/en/latest/?)

We recommend newcomers to the code read one or both of the following introductory
guides:

- [Getting Started Guide](https://geosx-geosx.readthedocs-hosted.com/en/latest/docs/sphinx/buildingGuide/DownloadingAndCompiling.html)
: Basic instructions for downloading, building, and running the code.

- [Developer Guide](https://geosx-geosx.readthedocs-hosted.com/en/latest/docs/sphinx/developerGuide/Index.html)
: Basic instructions for those interested in contributing code.

Once you are familiar with the basics, you may want to explore our tutorials page:
- [Tutorials](https://geosx-geosx.readthedocs-hosted.com/en/latest/docs/sphinx/tutorials/Index.html) 

We have two categories of tutorials.  User-focused tutorials for those primarily interested
in running existing models, and developer-focused tutorials for those primarily interested
in adding new capabilities to the code.

Release
-------

For release details and restrictions, please read the [LICENSE](./LICENSE) file.

For copyrights, please read the [COPYRIGHT](./COPYRIGHT ) file.

For contributors, please read the [CONTRIBUTORS](./CONTRIBUTORS ) file.

For acknowledgements, please read the [ACKNOWLEDGEMENTS](./ACKNOWLEDGEMENTS ) file.

For notice, please read the [NOTICE](./NOTICE ) file.

`LLNL-Code-746361`  `OCEC-18-021`
