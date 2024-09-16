---
title: 'GEOS: A performance portable multi-physics simulation framework for subsurface applications'
tags:
  - reservoir simulations
  - computational mechanics
  - multiphase flow
  - C++
authors:
  - name: Randolph R. Settgast
    orcid: 0000-0002-2536-7867
    corresponding: true
    affiliation: 1
  - name: Ryan M. Aronson
    orcid: 0009-0004-0785-5084
    affiliation: "2,3"
  - name: Julien R. Besset
    affiliation: 7
  - name: Andrea Borio
    orcid: 0000-0003-2016-5403
    affiliation: 5
  - name: Quan M. Bui
    orcid: 0000-0003-2648-0586
    affiliation: 1
  - name: Thomas J. Byer
    affiliation: 1
  - name: Nicola Castelletto
    orcid: 0000-0001-6816-6769
    affiliation: 1
  - name: Aurélien Citrain
    affiliation: 7
  - name: Benjamin C. Corbett
    orcid: 0009-0008-7108-9651
    affiliation: 1
  - name: James Corbett
    affiliation: 1
  - name: Philippe Cordier
    orcid: 0000-0002-6439-9263
    affiliation: 2
  - name: Matthias A. Cremon
    orcid: 0000-0001-7458-6401
    affiliation: 1
  - name: Cameron M. Crook
    orcid: 0000-0002-5366-6418
    affiliation: 1
  - name: Matteo Cusini
    orcid: 0000-0002-6024-861X
    affiliation: 1
  - name: Fan Fei
    orcid: 0000-0001-7273-4458
    affiliation: 1
  - name: Stefano Frambati
    orcid: 0000-0003-0683-1203
    affiliation: 7
  - name: Jacques Franc
    orcid: 0000-0002-8833-9425
    affiliation: 3
  - name: Andrea Franceschini
    orcid: 0000-0003-4395-5125
    affiliation: 3
  - name: Matteo Frigo
    orcid: 0000-0001-8150-1090
    affiliation: 3
  - name: Pengcheng Fu
    orcid: 0000-0002-7408-3350
    affiliation: 1
  - name: Thomas Gazzola
    orcid: 0000-0002-6103-4605
    affiliation: 2
  - name: Herve Gross
    orcid: 0000-0002-1747-2018
    affiliation: 2
  - name: Francois Hamon
    orcid: 0000-0001-8229-963X
    affiliation: 2
  - name: Brian M. Han
    orcid: 0009-0002-8549-7644
    affiliation: 1
  - name: Rasim Hasanzade
    affiliation: "3,4"
  - name: Michael Homel
    orcid: 0000-0002-0399-0092
    affiliation: 1
  - name: Jian Huang
    orcid: 0000-0002-5380-2563
    affiliation: 2
  - name: Tao Jin
    orcid: 0000-0001-6658-8941
    affiliation: 1
  - name: Isaac Ju
    orcid: 0000-0003-4110-7472
    affiliation: 3
  - name: Dickson Kachuma
    affiliation: 2
  - name: Mohammad Karimi-Fard
    orcid: 0000-0001-5707-165X
    affiliation: 3
  - name: Taeho Kim
    affiliation: 2
  - name: Sergey Klevtsov
    orcid: 0000-0001-9044-1827
    affiliation: 3
  - name: Alexandre Lapene
    affiliation: 2
  - name: Victor A. P. Magri
    orcid: 0000-0002-3389-523X
    affiliation: 1
  - name: Antoine Mazuyer
    orcid: 0000-0002-0329-3385
    affiliation: "2,3"
  - name: Mamadou N'diaye
    affiliation: 3
  - name: Daniel Osei-Kuffuor
    orcid: 0000-0002-6111-6205
    affiliation: 1
  - name: Stefan Povolny
    affiliation: 1
  - name: Guotong Ren
    orcid: 0000-0002-5821-9158
    affiliation: 4
  - name: Shabnam J. Semnani
    affiliation: 6
  - name: Chris S. Sherman
    orcid: 0000-0003-3550-0657
    affiliation: 1
  - name: Melvin Rey
    affiliation: 8
  - name: Hamdi A. Tchelepi
    orcid: 0000-0002-3084-6635
    affiliation: 3
  - name: William R. Tobin
    orcid: 0009-0001-3960-6064
    affiliation: 1
  - name: Pavel Tomin
    affiliation: 4
    orchid: 0000-0003-4862-4288
  - name: Lionel Untereiner
    affiliation: 8
  - name: Arturo Vargas
    affiliation: 1
  - name: Sohail Waziri
    affiliation: "3,4"
  - name: Xianhuan Wen
    affiliation: 4
  - name: Joshua A. White
    orcid: 0000-0003-3491-142X
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
  - name: University of California San Diego
    index: 6
  - name: Inria, Universite de Pau et des Pays de l’Adour
    index: 7
  - name: Independent
    index: 8

date: 28 May 2024
bibliography: paper.bib

---

# Summary

GEOS is a simulation framework focused on solving tightly coupled multi-physics problems with an initial emphasis on subsurface reservoir applications.
Currently, GEOS supports capabilities for studying carbon sequestration, geothermal energy, hydrogen storage, and related subsurface applications.
The unique aspect of GEOS that differentiates it from existing reservoir simulators is the ability to simulate tightly coupled compositional flow, poromechanics, fault slip, fracture propagation, and thermal effects, etc.
Extensive documentation is available on the [GEOS documentation pages](https://geosx-geosx.readthedocs-hosted.com/en/latest) [@GEOS_RTD].
Note that GEOS, as presented here, is a complete rewrite of the previous incarnation of the GEOS referred to in [@Settgast:2017].
<!--Additionally, this project has no relation to the computational geometry tool bearing the same name [@libgeos].-->

# Statement of need

The threat of climate change has resulted in an increased focus on mitigating carbon emissions into the atmosphere.
Carbon Capture and Storage (CCS) of CO~2~ in subsurface reservoirs and saline aquifers is an important component in the strategy to meet global climate goals.
Given the 2050 net-zero emissions goals, CO~2~ storage capacities required to offset emissions is orders of magnitude greater than current levels [@IPCC_2023].
Evaluation of reservoir performance and containment risks associated with the injection of liquefied CO~2~ in the subsurface in a reproducible and transparent manner is an important consideration when evaluating new storage sites.
As an example of typical complexities in carbon storage reservoirs, the 11th Society of Petroleum Engineers Comparative Solution Project (SPE11) [@Nordbotten2024] provides a benchmark example for evaluating the the predictions of carbon storage simulators.
The goal of GEOS is to provide the global community with an exascale capable open-source tool that is capable of simulating the complex coupled physics that occurs when liquefied CO~2~ is injected into a subsurface reservoir.
To this end, GEOS is freely available and focused on the simulation of reservoir integrity through various failure mechanisms such as caprock failure, fault leakage, and wellbore failure.
Open source projects such as OPM [@RASMUSSEN2021159], OpenGeoSys [@ogs:6.5.2], DuMux [@Kochetal2020Dumux] and Darts [@Voskov2024] are example efforts that share similar objectives.
However, GEOS stands out in two key areas: the explicit fault modeling coupled with flow and mechanical deformation, and the focus on  performance portability on platforms ranging from workstations to exascale supercomputers.




# GEOS Components

The core C++17 infrastructure provides common computer science capabilities typically required for solving differential equations using a spatially discrete method.
The components of the infrastructure provided by GEOS include a data hierarchy, a discrete mesh data structure, a mesh-based MPI communications interface, degree-of-freedom management, IO services, and a physics package interface.

By design, GEOS is intended to be a generic multi-physics simulation platform.
The physics package interface in GEOS is intended to encapsulate the development of numerical methods applied to the solution of governing equations relevant to a problem.
When implementing a physics package for a set of coupled physics equations, each individual physics package is first developed as a stand-alone capability.
The single physics capabilities are then applied together in a coupled physics package and solved through a flexible strategy ranging from solving the fully monolithic system to a split operator approach.

To solve the linear systems that arise from the boundary value problem, GEOS maintains a generic linear algebra interface (LAI) capable of wrapping various linear algebra packages such as hypre [@hypre], PETSc [@petsc-web-page], and Trilinos [@Her_etal05].
Currently, in GEOS only the hypre interface is actively maintained.
For every multi-physics problem involving the solution of a coupled linear system, GEOS currently relies on a multigrid reduction preconditioning strategy available in hypre [@BUI:2020;@BUI:2021114111].

The performance portability strategy utilized by GEOS applies LLNL's suite of portability tools RAJA [@Beckingsale:2019], CHAI [@CHAI:2023], and Umpire [@Beckingsale:2020].
The RAJA performance portability layer provides [performance portable](https://performanceportability.org) kernel launching and wrappers for reductions, atomics, and local/shared memory to achieve performance on both CPU and GPU hardware.
The combination of CHAI/Umpire provides memory motion management for platforms with heterogeneous memory spaces (i.e., host and device memory).
Through this strategy, GEOS has been successfully run on platforms ranging from GPU-based Exa-scale systems to CPU-based laptops with near-optimal of performance.

In addition to its C++ core, the GEOS project provides a Python3 interface that allows for the integration of the simulation capabilities into complex Python workflows involving components unrelated to GEOS.

# Applications
To date GEOS has been used to simulate problems relevant to CO~2~ storage, enhanced geothermal systems, hydrogen storage, and both conventional and unconventional oil and gas extraction.
Often these simulations involve coupling between compositional multiphase flow and transport, poroelasticity, thermal transport, and interactions with faults and fractures.

As an example of a field case where GEOS has been applied, we present a coupled compositional flow/mechanics simulation of CO~2~ injection and storage at a large real-world storage site.
Figure \ref{RW_results}a illustrates the computational mesh and Figure \ref{RW_results}b shows results after 25 years of injection.
Simulations such as this will play a critical role in predicting the viability of potential CO~2~ storage sites.

![Real world CO~2~ storage site: (a) discrete mesh, transparency is used for the overburden region to reveal the complex faulted structure of the storage reservoir; (b) results of a compositional flow simulation after 25 years of CO~2~ injection. The CO~2~ plume is shown in white near the bottom of the well. Colors in the reservoir layer indicate changes in fluid pressure, and the colors in the overburden indicate vertical displacement resulting from the injection. Note that color scales have been removed intentionally.\label{RW_results}](RW_final.pdf){ width=100% }

As an example of the weak scalability of GEOS on an exascale class system, we present two weak scaling studies on a simple wellbore geometry run on the Frontier supercomputer at Oak Ridge National Laboratory.
Frontier is comprised of 9,472 Cray EX235a nodes that each contain a single AMD EPYC 7A53 CPU and four AMD MI250X GPUs [@frontier].
Note that each MI250X is comprised of two Graphics Compute Dies (GCD), with each GCD appearing as a GPU to the operating system. 
A more detailed discussion and instructions to reproduce the results are available in the [Performance Benchmarks](https://geosx-geosx.readthedocs-hosted.com/en/latest/docs/sphinx/advancedExamples/performanceBenchmarks/Index.html) of the GEOS documentation.

The weak scaling results for mechanics are presented in (Figure \ref{fig:Frontier_scaling}a) and shows nearly flat scaling of the GEOS processes (assembly/field synchronization) up to 32,768 GPUs ($81.3 \times 10^{9}$ degrees-of-freedom).
There is a moderate decrease in efficiency with the application of the hypre preconditioner setup and solve, but given the complexity of those algorithms, this level of scaling efficiency is excellent.
The weak scaling results of compositional flow are presented in Figure \ref{fig:Frontier_scaling}b shows excellent scaling up to 2,048 GPUs. 

![Weak scaling results on ORNL/Frontier: average execution time per newton iteration vs number of GPUs for a mechanics (a) and a compositional flow (b) simulation, respectively.\label{fig:Frontier_scaling}](nearwell_scaling_frontier.pdf){ width=100% }

# Acknowledgements
This work was performed under the auspices of the U.S. Department of Energy by Lawrence Livermore National Laboratory under Contract DE-AC52-07NA27344. LLNL release number LLNL-JRNL-864747.

This research was supported by the Exascale Computing Project (ECP), Project Number: 17-SC-20-SC, a collaborative effort of two DOE organizations - the Office of Science and the National Nuclear Security Administration, responsible for the planning and preparation of a capable exascale ecosystem, including software, applications, hardware, advanced system engineering and early testbed platforms, to support the nation's exascale computing imperative.

Support was provided by TotalEnergies and Chevron through the FC-MAELSTROM project, a collaborative effort between Lawrence Livermore National Laboratory, TotalEnergies, Chevron, and Stanford University, aiming to develop an exascale compatible, multiscale, research-oriented simulator for modeling fully coupled flow, transport and geomechanics in geological formations.

# References
