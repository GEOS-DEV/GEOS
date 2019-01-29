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
[publications](https://github.com/GEOSX/GEOSX/blob/develop/src/docs/sphinx/publications.rst)
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
[authors](https://github.com/GEOSX/GEOSX/blob/develop/src/docs/sphinx/authors.rst)
and
[acknowledgements](https://github.com/GEOSX/GEOSX/blob/develop/src/docs/sphinx/acknowledgements.rst)
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
We recommend newcomers to the code read one or both of the following introductory
guides:

- [Getting Started Guide](https://github.com/GEOSX/GEOSX/blob/develop/src/docs/sphinx/getting_started.rst)
: Basic instructions for downloading, building, and running the code.

- [Developer Guide](https://github.com/GEOSX/GEOSX/blob/develop/src/docs/sphinx/developer_guidelines.rst)
: Basic instructions for those interested in contributing code.

Once you are familiar with the basics, you may want to explore our tutorials page:
- [Tutorials](https://github.com/GEOSX/GEOSX/blob/develop/src/docs/sphinx/tutorials.rst) 

We have two categories of tutorials.  User-focused tutorials for those primarily interested
in running existing models, and developer-focused tutorials for those primarily interested
in adding new capabilities to the code.

Release
-------
Copyright (c) 2018, Lawrence Livermore National Security, LLC.

Produced at the Lawrence Livermore National Laboratory.

All rights reserved.

Unlimited Open Source - LGPL 2.1 License

For release details and restrictions, please read the [LICENSE](./LICENSE) file.

`LLNL-Code-746361`  `OCEC-18-021`
