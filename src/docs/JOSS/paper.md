---
title: 'GEOS: A portable multi-physics simulation framework'
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
  - name: Arturo Vargas
    affiliation: 1
  - name: Joshua White
    affiliation: 1
  - name: William R. Tobin
    affiliation: 1
  - name: Brian M. Han
    affiliation: 1
  - name: Herve Gross
    affiliation: 2
  - name: Stefan Framba
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

GEOS is an open-source simulation framework focused on implementing tightly-coupled multi-physics problems with an emphasis subsurface reservoir applications.
The C++ infrastructure of GEOS provides facilities to assist in the implementation of constraint equations such as a discrete mesh data structure, MPI communications tools, degree-of-freedom management, IO facilities, etc.
The performance portability strategy applies LLNL's suite of portability tools RAJA[@Beckingsale:2019], CHAI[@CHAI:2023], and Umpire[@Beckingsale:2020].
Through this strategy GEOS has been successfully run on platforms ranging from GPU-based Exa-scale systems down to CPU-based laptops.
The code provides implementations for a number of target applications, including, carbon sequestration, geothermal energy, hydrogen storage, and similar systems.
The Python interface of GEOS allows for the integration of the computational simulation capabilities in GEOS into complex python workflows involving other tools.

# Statement of need

`Gala` is an Astropy-affiliated Python package for galactic dynamics. Python
enables wrapping low-level languages (e.g., C) for speed without losing
flexibility or ease-of-use in the user-interface. The API for `Gala` was
designed to provide a class-based and user-friendly interface to fast (C or
Cython-optimized) implementations of common operations such as gravitational
potential and force evaluation, orbit integration, dynamical transformations,
and chaos indicators for nonlinear dynamics. `Gala` also relies heavily on and
interfaces well with the implementations of physical units and astronomical
coordinate systems in the `Astropy` package [@astropy] (`astropy.units` and
`astropy.coordinates`).

`Gala` was designed to be used by both astronomical researchers and by
students in courses on gravitational dynamics or astronomy. It has already been
used in a number of scientific publications [@Pearson:2017] and has also been
used in graduate courses on Galactic dynamics to, e.g., provide interactive
visualizations of textbook material [@Binney:2008]. The combination of speed,
design, and support for Astropy functionality in `Gala` will enable exciting
scientific explorations of forthcoming data releases from the *Gaia* mission
[@gaia] by students and experts alike.

# Infrastructure Components 

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

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements


This work performed under the auspices of the U.S. Department of Energy by Lawrence Livermore National Laboratory under Contract DE-AC52-07NA27344

This research was supported by the Exascale Computing Project (ECP), Project Number: 17-SC-20-SC, a collaborative effort of two DOE organizations - the Office of Science and the National Nuclear Security Administration, responsible for the planning and preparation of a capable exascale ecosystem, including software, applications, hardware, advanced system engineering and early testbed platforms, to support the nation's exascale computing imperative.

FC-MAELSTROM statement...

# References