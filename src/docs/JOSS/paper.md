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
  - name: Stefano Frambati
    affiliation: 2
  - name: Aurélien Citrain
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
Thus, GEOS is freely available and focused on the simulation of reservoir integrity through various failure mechanisms such as caprock failure, fault leakage, and wellbore failure.
Additionally GEOS provides the potential to estimate seismic events induced by CO2 injection.

# C++ Infrastructure Components 

The core c++17 infrastructure provides components to perform common computer science tasks typical of a discrete numerical simulation. 
The components of the infrastructure provided by GEOS include a data hierarchy, a discrete mesh data structure, a physics package interface, mesh based MPI communications interface, degree-of-freedom management, IO facilities, and an event manager.

The data repository forms an object hierarchy through the definition of a `Wrapper` class and a `Group` class. 
Drawing an analogy with a standard folder/file hierarchy, the `Group` class can be thought of as a "Folder" as it holds other `Group`'s as well as a collection of `Wrapper` objects, while the `Wrapper` is a container for any arbitrary object.
The mesh interface is built on top of the data repository as a collection of object managers for each mesh object type (e.g. node, edge, face, element).
The management of distributed memory parallelism is done through MPI, and the execution of distributed memory parallelism typically requires little intervention from the physics developer.

The performance portability strategy utilized by GEOS applies LLNL's suite of portability tools RAJA[@Beckingsale:2019], CHAI[@CHAI:2023], and Umpire[@Beckingsale:2020].
The RAJA performance portability layer provides portable kernel launching and wrappers for reductions, atomics, and local/shared memory to achieve performance on both CPU and GPU hardware.
The combination of CHAI/Umpire provides memory motion management for platforms with heterogeneous memory spaces (i.e. host memory and device memory).
Through this strategy GEOS has been successfully run on platforms ranging from GPU-based Exa-scale systems to CPU-based laptops with minimal loss of performance due to platform changes.

GEOS maintains a generic linear algebra interface (LAI) capable of wrapping various linear algebra packages.
However as a matter of practice, the primary linear algebra package used for the great majority of GEOS simulations is LLNL's hypre[@hypre].

In addition to its c++ core, the the GEOS team provides a Python3 interface that allows for the integration of the simulation capabilities into complex python workflows involving components unrelated to GEOS.
The Python3 interface provides data exchange between GEOS simulations and the Python driver, as well as allowing the Python layer to call specific GEOS packages outside of standard GEOS c++ workflow.

GEOS is intended to be a generic multi-physics simulation platform.
As such, single physics capabilities are developed and tested independently. 
When coupling one or more single physics capabilities together to create a couple physics package, the strategy can be expressed as a monolithic linear system with an underlying block structure corresponding where the row/col of the block corresponds with a set of constraint equations/degrees-of-freedom associated with a physics package.
In figure \autoref{fig:matrix}, a system matrix ($A_{coupled}$) that couples 3 distinct physics packages (1,2,3) is shown.

the diagonal blocks result from each single physics package.
The off-diagonal blocks represent the coupling between physics packages and are typically filled through various options, such as through the coupled physics package, or through a callback mechanism in the single physics package.

\begin{equation}
A_{coupled} = \left(\begin{array}{@{}c|c|c@{}}
    A_{11} & A_{12} & A_{13}  \\ \hline
    A_{21} & A_{22} & A_{23}  \\ \hline
    A_{31} & A_{32} & A_{33}  \\ 
\end{array}\right)
\end{equation}\label{fig:matrix}

# Applications
The development of GEOS targets multi-physics simulations of subsurfaces reservoirs.
To date GEOS has been used to simulate problems relevant to CO2 storage, enhanced geothermal systems, hydrogen storage, and both conventional and unconventional oil and gas extraction.
Often these simulations involve coupling between compositional multiphase flow and transport, poroelasticity, thermal transport, and interactions with faults and fractures.

The coupling strategy applied in GEOS is to require the capability of a tightly coupled monolithic system as a baseline capability.
In cases where such tight coupling is not required, one may decompose the monolithic system into blocks and apply a sequential coupling approach.


# Acknowledgements


This work performed under the auspices of the U.S. Department of Energy by Lawrence Livermore National Laboratory under Contract DE-AC52-07NA27344

This research was supported by the Exascale Computing Project (ECP), Project Number: 17-SC-20-SC, a collaborative effort of two DOE organizations - the Office of Science and the National Nuclear Security Administration, responsible for the planning and preparation of a capable exascale ecosystem, including software, applications, hardware, advanced system engineering and early testbed platforms, to support the nation's exascale computing imperative.

FC-MAELSTROM statement...

# References