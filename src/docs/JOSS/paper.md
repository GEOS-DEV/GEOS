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
  - name: Ryan M. Aronson
    affiliation: 3
  - name: Julien R. Besset
    affiliation: 2
  - name: Andrea Borio
    affiliation: 5
  - name: Thomas J. Byer
    affiliation: 1
  - name: Nicola Castelletto
    affiliation: 1
  - name: Aurélien Citrain
    affiliation: 2
  - name: Benjamin C. Corbett
    affiliation: 1
  - name: James Corbett
    affiliation: 1
  - name: Philippe Cordier
    affiliation: 2
  - name: Matthias A. Cremon
    affiliation: 1
  - name: Cameron M. Crook
    affiliation: 1
  - name: Matteo Cusini
    affiliation: 1
  - name: Fan Fei
    affiliation: 1
  - name: Stefano Frambati
    affiliation: 2
  - name: Andrea Franceschini
    affiliation: 3
  - name: Matteo Frigo
    affiliation: 3
  - name: Thomas Gazzola
    affiliation: 2
  - name: Herve Gross
    affiliation: 2
  - name: Francois Hamon
    affiliation: 2
  - name: Brian M. Han
    affiliation: 1
  - name: Michael Homel
    affiliation: 1
  - name: Jian Huang
    affiliation: 2
  - name: Tao Jin
    affiliation: 1
  - name: Dickson Kachuma
    affiliation: 2
  - name: Mohammad Karimi-Fard
    affiliation: 3
  - name: Sergey Klevtsov
    affiliation: 3
  - name: Alexandre Lapene
    affiliation: 2
  - name: Victor A. P. Magri
    affiliation: 1
  - name: Daniel Osei-Kuffuor
    affiliation: 1
  - name: Stefan Povolny
    affiliation: 1
  - name: Shabnam J. Semnani
  - name: Chris S. Sherman
    affiliation: 1
  - name: Melvin Rey
    affiliation: 2
  - name: Hamdi A. Tchelepi
    affiliation: 3
  - name: William R. Tobin
    affiliation: 1
  - name: Pavel Tomin
    affiliation: 4
  - name: Lionel Untereiner
    orcid: 0000-0002-8025-2616
  - name: Joshua A. White
    affiliation: 1
  - name: Hui Wu
    affiliation: 1
affiliations:
 - name: Lawrence Livermore National Laboratory, USA
   index: 1
 - name: TotalEnergies E&P Research & Technology, USA
   index: 2
 - name: Stanford University, USA
   index: 3
 - name: Chevron Technical Center, USA
   index: 4
 - name: Politecnico di Torino, Italy
   index: 5
date: 15 December 2023
bibliography: paper.bib

---

# Summary

GEOS is a simulation framework focused solving tightly-coupled multi-physics problems with an initial emphasis subsurface reservoir applications.
Currently GEOS actively supports implementations for studying carbon sequestration, geothermal energy, hydrogen storage, and similar subsurface applications.
The unique aspect of GEOS that differentiates it from existing reservoir simulators is the ability to provide tightly-coupled compositional flow, poromechanics, faults and fractures slip, and thermal effects, etc.
The GEOS repository is located at https://github.com/GEOS-DEV/GEOS, and extensive documentation for GEOS is available at https://geosx-geosx.readthedocs-hosted.com/en/latest. 
Note that the version of GEOS described here should be considered a separate work form the previous incarnation of GEOS referred to in [@Settgast:2017].

# Statement of need

The increasing threat of climate change has resulted in an increased focus on mitigating carbon emissions into the atmosphere.
Carbon Capture and Storage (CCS) of CO2 in subsurface reservoirs and saline aquifers is an important component in the strategy to meet global climate goals. 
Given the 2050 net-zero GHG goals, CO2 storage capacities required to offset emissions is orders of magnitude greater than current levels [@IPCC_2023].
The ability to evaluate the reservoir performance and containment risks associated with the injection of liquefied CO2 in the subsurface in a reproducible and transparent manner is an important consideration when developing new storage sites.
The primary goal of GEOS is to provide the global community with an open-source tool that is capable of simulating the complex coupled physics that occurs when liquefied CO2 is injected into a subsurface reservoir. 
Thus, GEOS is freely available and focused on the simulation of reservoir integrity through various failure mechanisms such as caprock failure, fault leakage, and wellbore failure.

# GEOS Components 

The core c++17 infrastructure provides common computer science capabilities typically required for solving differential equations using a spatially discrete method. 
The components of the infrastructure provided by GEOS include a data hierarchy, a discrete mesh data structure, a mesh based MPI communications interface, degree-of-freedom management, IO services, and a physics package interface.

By design, GEOS is intended to be a generic multi-physics simulation platform.
The physics package interface in GEOS is intended to encapsulate the development of numerical methods applied to the solution of governing equations relevant to a problem.
When implementing a physics package for a set of coupled physics equations, each individual physics package is first developed as a stand-alone capability. 
The single physics capabilities are then applied together in a coupled physics package and solved through a flexible strategy ranging from solving the fully monolithic system, to a split operator approach. 

To solve the linear systems that arise from the boundary value problem, GEOS maintains a generic linear algebra interface (LAI) capable of wrapping various linear algebra packages such as hypre [@hypre], PETSc[@petsc-web-page], and Trilinos[@trilinos-website].
Currently, in GEOS only the hypre interaface is actively maintained.
For every multi-physics problems involving the solution of a coupled linear system, GEOS currently relies on a multigrid reduction preconditioning strategy available in hypre as presented by [@BUI:2020;@BUI:2021114111].

The performance portability strategy utilized by GEOS applies LLNL's suite of portability tools RAJA[@Beckingsale:2019], CHAI[@CHAI:2023], and Umpire[@Beckingsale:2020].
The RAJA performance portability layer provides portable kernel launching and wrappers for reductions, atomics, and local/shared memory to achieve performance on both CPU and GPU hardware.
The combination of CHAI/Umpire provides memory motion management for platforms with heterogeneous memory spaces (i.e. host memory and device memory).
Through this strategy GEOS has been successfully run on platforms ranging from GPU-based Exa-scale systems to CPU-based laptops with near optimal of performance.

In addition to its c++ core, the GEOS project provides a Python3 interface that allows for the integration of the simulation capabilities into complex python workflows involving components unrelated to GEOS.

# Applications
To date GEOS has been used to simulate problems relevant to CO2 storage, enhanced geothermal systems, hydrogen storage, and both conventional and unconventional oil and gas extraction.
Often these simulations involve coupling between compositional multiphase flow and transport, poroelasticity, thermal transport, and interactions with faults and fractures.

As an example of a field case where GEOS has been applied, we present a coupled compositional flow/mechanics simulation of CO2 injection and storage at a large real-world storage site.
Figure \ref{RW_final}a illustrates the computational mesh and Figure \ref{RW_final}b shows results after 25 years of injection.
Simulations such as this will play a critical role in predicting the viability of potential CO2 storage sites.

![Real world CO2 storage site: (a) discrete mesh, transparency is used for the overburden region to reveal the complex faulted structure of the storage reservoir; (b) results of a compositional flow simulation after 25 years of CO2 injection. The CO2 plume is shown in white near the bottom of the well. Colors in the reservoir layer indicate changes in fluid pressure, and the colors in the overburden indicate vertical displacement resulting from the injection. Note that color scales have been removed intentionally.\label{RW_results}](RW_final.pdf){ width=100% }

As an example of the weak scalability of GEOS on exascale systems, we present two weak scaling studies on a simple wellbore geometry using the exascale Frontier supercomputer located at Oak Ridge National Laboratory (ORNL).
The results from the weak scaling study (Figure \ref{fig:Frontier_scaling}a) shows flat scaling of the GEOS processes (assembly/field synchronization) up to 16,384 MPI ranks and 81.3e9 degrees-of-freedom (1/4 of Frontier).
There is a moderate decrease in efficiency with the application of the hypre preconditioner setup and solve, but given the complexity of those algorithms this level of scaling efficiency is excellent.
The compositional flow study presented in Figure \ref{fig:Frontier_scaling}b shows similarly good weak scaling. 

![Weak scaling results on ORNL/Frontier: execution time per timestep vs number of cluster ranks for a mechanics (a) and a compositional flow (b) simulation, respectively.\label{fig:Frontier_scaling}](GEOS_Frontier_scaling.pdf){ width=100% }

# Acknowledgements
This work performed under the auspices of the U.S. Department of Energy by Lawrence Livermore National Laboratory under Contract DE-AC52-07NA27344

This research was supported by the Exascale Computing Project (ECP), Project Number: 17-SC-20-SC, a collaborative effort of two DOE organizations - the Office of Science and the National Nuclear Security Administration, responsible for the planning and preparation of a capable exascale ecosystem, including software, applications, hardware, advanced system engineering and early testbed platforms, to support the nation's exascale computing imperative.

Partial funding was provided by TotalEnergies and Chevron through the FC-MAELSTROM project, a collaborative effort between Lawrence Livermore National Laboratory, TotalEnergies, Chevron, and Stanford University, aiming to develop an exascale compatible, multiscale, research-oriented simulator for modeling fully coupled flow, transport and geomechanics in geological formations.

# References
