---
title: 'GEOS-2023: A portable multi-physics simulation framework'
tags:
  - reservoir simulations
  - computational mechanics
  - multiphase flow
  - c++
authors:
  - name: Randolph R. Settgast
    orcid: 0000-0002-2536-7867
    corresponding: true
    affiliation: 1
  - name: Benjamin C. Corbett
    affiliation: 1
  - name: James Corbett
    affiliation: 1
  - name: Francois Hamon
    affiliation: 2
  - name: Thomas Gazzola
    affiliation: 2
  - name: Matteo Cusini
    affiliation: 1
  - name: Chris S. Sherman
    affiliation: 1
  - name: Sergey Klevzoff
    affiliation: 3
  - name: Nicola Castelletto
    affiliation: 1
  - name: Victor Paludetto Magri
    affiliation: 1
  - name: William R. Tobin
    affiliation: 1
  - name: Joshua White
    affiliation: 1
  - name: Brian M. Han
    affiliation: 1
  - name: Herve Gross
    affiliation: 2
  - name: Stefano Framba
    affiliation: 2
  - name: Aurilian Citrain
    affiliation: 2
  - name: Andrea Franceschini
    affiliation: 2
  - name: Andrea Borio
    affiliation: 4
  - name: Jian Huang
    affiliation: 2
affiliations:
 - name: Lawrence Livermore National Laboratory, USA
   index: 1
 - name: TotalEnergies E&P Research & Technology, USA
   index: 2
 - name: Stanford University, USA
   index: 3
 - name: Politecnico di Torino, Italy
date: 15 December 2023
bibliography: paper.bib

---

# Summary

GEOS is a simulation framework focused on implementing tightly-coupled multi-physics problems with an initial emphasis subsurface reservoir applications.
Specifically, GEOS provides implementations for studying carbon sequestration, geothermal energy, hydrogen storage, and similar problems, and allows developers to easily extend or add new formulations to the suite of capabilities.
The unique aspect of GEOS that differentiates it from existing reservoir simulators is the ability to provide tightly-coupled compositional flow, poromechanics, faults and fractures, and thermal effects.
Extensive documentation for GEOS is available at https://geosx-geosx.readthedocs-hosted.com/en/latest. 

# Statement of need

The increasing threat of climate change has resulted in an increased focus on mitigating carbon emissions into the atmosphere.
Carbon Capture and Storage (CCS) of CO2 in subsurface reservoirs and saline aquifers is one of the most important technologies required to meet global climate goals. 
Given the 2050 net-zero GHG goals, CO2 storage capacities required to offset emissions is orders of magnitude greater than current levels.(reference needed)
One factor in the evaluation of CO2 storage sites are the containment risks associated with the injection of liquefied CO2 in the subsurface.
The primary goal of GEOS is to provide the global community with an open-source tool that is capable of simulating the complex coupled physics that occurs when liquefied CO2 is injected into a subsurface reservoir. 
Thus, GEOS is a freely available tool that is focused on the simulation of reservoir integrity through various failure mechanisms such as caprock failure, fault leakage, and wellbore failure.
Additionally GEOS provides the potential to estimate seismic events induced by CO2 injection.

# C++ Infrastructure Components 

The core c++17 infrastructure provides components to perform common computer science tasks that are required in a simulation code. 
The components of the infrastructure include a data hierarchy, a discrete mesh data structure, a physics package interface, MPI communications tools, degree-of-freedom management, IO facilities, and an event manager.

The data repository defines a `Wrapper` class to hold anything from data arrays to arbitrary objects, and a `Group` class that serves as a container to form a hierarchy.
Drawing an analogy with a standard folder/file hierarchy, the `Group` class can be thought of as a "Folder" as it holds other `Group`'s as well as a collection of `Wrapper` objects. 
The `Wrapper` can be thought of as a "File" as it contains the relevant data that is stored in the repository.
The mesh interface is built on top of the data repository as a collection of managers for each mesh object type as shown in Figure \autoref{fig:meshHierarchy}.
On each MPI rank there is a `MeshBody` object that represents a physical body.
Each `MeshBody` contains a collection of `MeshLevel` objects that represent a discretization of the `MeshBody`.
Each `MeshLevel` holds a collection of managers objects that contain data on each type of discrete mesh object (i.e. nodes, edges, faces, elements).
The role of each mesh object manager is to hold maps between the mesh objects, and to hold field/dof data.

![UML diagram of the mesh interface hierarchy.\label{fig:meshHierarchy}](MeshHierarchy.svg){ width=40% }

The performance portability strategy utilized by GEOS applies LLNL's suite of portability tools RAJA[@Beckingsale:2019], CHAI[@CHAI:2023], and Umpire[@Beckingsale:2020].
The RAJA performance portability layer provides portable kernel launching and wrappers for reductions, atomics, and local/shared memory to achieve performance on both CPU and GPU hardware.
The combination of CHAI/Umpire provides memory motion management for platforms with heterogeneous memory spaces (i.e. host memory and device memory).
Through this strategy GEOS has been successfully run on platforms ranging from GPU-based Exa-scale systems to CPU-based laptops.


In addition to the c++ core, GEOS maintains a Python3 interface that allows for the integration of the simulation capabilities into complex python workflows involving components unrelated to GEOS.
The Python3 interface provides data exchange between GEOS simulations and the Python driver, as well as allowing the Python layer to call specific GEOS packages outside of standard GEOS event manager workflow.

# Applications
The development of GEOS specifically targets simulation of subsurfaces reservoirs.
To date GEOS has been used to simulate problems relevant to CO2 storage, enhanced geothermal systems, hydrogen storage, and both conventional and unconventional oil and gas extraction.

However, GEOS is intended to be a generic multi-physics simulation platform.

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Acknowledgements


This work performed under the auspices of the U.S. Department of Energy by Lawrence Livermore National Laboratory under Contract DE-AC52-07NA27344

This research was supported by the Exascale Computing Project (ECP), Project Number: 17-SC-20-SC, a collaborative effort of two DOE organizations - the Office of Science and the National Nuclear Security Administration, responsible for the planning and preparation of a capable exascale ecosystem, including software, applications, hardware, advanced system engineering and early testbed platforms, to support the nation's exascale computing imperative.

FC-MAELSTROM statement...

# References