/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file EigenSolverParameters.hpp
 */

#ifndef GEOS_LINEARALGEBRA_UTILITIES_EIGENSOLVERPARAMETERS_HPP_
#define GEOS_LINEARALGEBRA_UTILITIES_EIGENSOLVERPARAMETERS_HPP_

#include "common/DataTypes.hpp"
#include "common/format/EnumStrings.hpp"

namespace geos
{

/**
 * @brief Parameters of a generalized symmetric eigensolver for the pencil (K, M), K x = lambda M x.
 */
struct EigenSolverParameters
{
  /// Available eigensolver algorithms
  enum class SolverType : integer
  {
    arnoldi,  ///< Block Krylov-Schur (thick-restart Arnoldi/Lanczos) with spectral transformation (K - shift M)^{-1} M
    lobpcg    ///< Locally optimal block preconditioned conjugate gradient (smallest eigenvalues, needs a preconditioner only)
  };

  /// Eigensolver algorithm
  SolverType solverType = SolverType::arnoldi;

  /// Number of eigenpairs requested (those closest to the shift)
  integer numEigenvalues = 10;

  /// Spectral shift sigma, in the units of the eigenvalue lambda (lambda = omega^2 for vibration modes)
  real64 shift = 0.0;

  /// Relative convergence tolerance on the Ritz residual estimate
  real64 tolerance = 1.0e-8;

  /// Maximum number of restarts (Arnoldi) or iterations (LOBPCG)
  integer maxIterations = 300;

  /// Maximum dimension of the Krylov basis (Arnoldi), or block size of the iterates (LOBPCG, if larger than
  /// numEigenvalues, the extra vectors accelerate convergence); 0 selects a default
  integer subspaceSize = 0;

  /// Number of vectors expanded at once (Arnoldi block size). Needs to be at least the multiplicity of the
  /// eigenvalues to find all of them without the completeness check.
  integer blockSize = 1;

  /// After convergence, inject a fresh random block orthogonal to the converged vectors to detect missed
  /// copies of repeated eigenvalues (e.g. the six rigid-body modes of a free structure)
  integer completenessCheck = 1;

  /// Seed of the random starting vectors
  integer seed = 1;

  /// Verbosity (0 = silent, 1 = one line per restart)
  integer logLevel = 0;
};

ENUM_STRINGS( EigenSolverParameters::SolverType,
              "arnoldi",
              "lobpcg" );

} // namespace geos

#endif /* GEOS_LINEARALGEBRA_UTILITIES_EIGENSOLVERPARAMETERS_HPP_ */
