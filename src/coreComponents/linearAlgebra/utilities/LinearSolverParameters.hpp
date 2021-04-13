/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file LinearSolverParameters.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_UTILITIES_LINEARSOLVERPARAMETERS_HPP_
#define GEOSX_LINEARALGEBRA_UTILITIES_LINEARSOLVERPARAMETERS_HPP_

#include "codingUtilities/EnumStrings.hpp"

namespace geosx
{

/**
 * @brief Set of parameters for a linear solver or preconditioner.
 *
 * This class holds a simple tree of linear solver options.
 * They are set to default values, but can be overwritten as needed.
 */
struct LinearSolverParameters
{
  /**
   * @brief Linear solver type.
   */
  enum class SolverType : integer
  {
    direct,        ///< Direct solver
    cg,            ///< CG
    gmres,         ///< GMRES
    fgmres,        ///< Flexible GMRES
    bicgstab,      ///< BiCGStab
    preconditioner ///< Preconditioner only
  };

  /**
   * @brief Preconditioner type.
   */
  enum class PreconditionerType : integer
  {
    none,   ///< No preconditioner
    jacobi, ///< Jacobi smoothing
    gs,     ///< Gauss-Seidel smoothing
    sgs,    ///< Symmetric Gauss-Seidel smoothing
    iluk,   ///< Incomplete LU with k-level of fill
    ilut,   ///< Incomplete LU with thresholding
    icc,    ///< Incomplete Cholesky
    ict,    ///< Incomplete Cholesky with thresholding
    amg,    ///< Algebraic Multigrid
    mgr,    ///< Multigrid reduction (Hypre only)
    block   ///< Block preconditioner
  };

  integer logLevel = 0;     ///< Output level [0=none, 1=basic, 2=everything]
  integer dofsPerNode = 1;  ///< Dofs per node (or support location) for non-scalar problems
  bool isSymmetric = false; ///< Whether input matrix is symmetric (may affect choice of scheme)
  integer stopIfError = 1;  ///< Whether to stop the simulation if the linear solver reports an error

  SolverType solverType = SolverType::direct;                        ///< Solver type
  PreconditionerType preconditionerType = PreconditionerType::iluk;  ///< Preconditioner type

  /// Direct solver parameters: used for SuperLU_Dist interface through hypre and PETSc
  struct Direct
  {
    /**
     * @brief How to permute the columns
     */
    enum class ColPerm : integer
    {
      none,        ///< natural
      MMD_AtplusA, ///< multiple minimum degree on At+A
      MMD_AtA,     ///< multiple minimum degree on At*A (heavy)
      colAMD,      ///< approximate minimum degree on columns
      metis,       ///< using METIS
      parmetis     ///< using ParMETIS
    };

    /**
     * @brief How to permute the rows
     */
    enum class RowPerm : integer
    {
      none, ///< natural
      mc64  ///< using HSL routine MC64
    };

    integer checkResidual = 0;        ///< Whether to check the linear system solution residual
    integer equilibrate = 1;          ///< Whether to scale the rows and columns of the matrix
    ColPerm colPerm = ColPerm::metis; ///< Columns permutation
    RowPerm rowPerm = RowPerm::mc64;  ///< Rows permutation
    integer replaceTinyPivot = 1;     ///< Whether to replace tiny pivots by sqrt(epsilon)*norm(A)
    integer iterativeRefine = 1;      ///< Whether to perform iterative refinement
    integer parallel = 1;             ///< Whether to use a parallel solver (instead of a serial one)
  }
  direct;                             ///< direct solver parameter struct

  /// Krylov-method parameters
  struct Krylov
  {
    real64 relTolerance = 1e-6;       ///< Relative convergence tolerance for iterative solvers
    integer maxIterations = 200;      ///< Max iterations before declaring convergence failure
    integer maxRestart = 200;         ///< Max number of vectors in Krylov basis before restarting
    integer useAdaptiveTol = false;   ///< Use Eisenstat-Walker adaptive tolerance
    real64 weakestTol = 1e-3;         ///< Weakest allowed tolerance when using adaptive method
  }
  krylov;                             ///< Krylov-method parameter struct

  /// Matrix-scaling parameters
  struct Scaling
  {
    integer useRowScaling = false;      ///< Apply row scaling
    integer useRowColScaling = false;   ///< Apply row and column scaling (not yet implemented)
  }
  scaling;                              ///< Matrix-scaling parameter struct

  /// Algebraic multigrid parameters
  struct AMG
  {
    integer maxLevels = 20;                  ///< Maximum number of coarsening levels
    string cycleType = "V";                  ///< AMG cycle type
    string smootherType = "gaussSeidel";     ///< Smoother type
    string coarseType = "direct";            ///< Coarse-level solver/smoother
    string coarseningType = "HMIS";          ///< Coarsening algorithm
    integer interpolationType = 6;           ///< Coarsening algorithm
    integer numSweeps = 2;                   ///< Number of smoother sweeps
    integer numFunctions = 1;                ///< Number of amg functions
    integer aggresiveNumLevels = 0;          ///< Number of levels for aggressive coarsening.
    string preOrPostSmoothing = "both";      ///< Pre and/or post smoothing [pre,post,both]
    real64 threshold = 0.0;                  ///< Threshold for "strong connections" (for classical and
                                             ///< smoothed-aggregation AMG)
    integer separateComponents = false;      ///< Apply a separate component filter before AMG construction
    string nullSpaceType = "constantModes";  ///< Null space type [constantModes,rigidBodyModes]
  }
  amg;                                       ///< Algebraic Multigrid (AMG) parameters

  /// Multigrid reduction parameters
  struct MGR
  {
    /**
     * @brief MGR available strategies
     */
    enum class StrategyType : integer
    {
      compositionalMultiphaseFVM,       ///< finite volume compositional muliphase flow
      compositionalMultiphaseHybridFVM, ///< hybrid finite volume compositional muliphase flow
      compositionalMultiphaseReservoir, ///< reservoir with finite volume compositional multiphase flow
      hydrofracture,                    ///< hydrofracture
      lagrangianContactMechanics,       ///< Lagrangian contact mechanics
      singlePhasePoroelastic,           ///< single phase poroelastic with finite volume single phase flow
      hybridSinglePhasePoroelastic      ///< single phase poroelastic with hybrid finite volume single phase flow
    };

    StrategyType strategy;              ///< Predefined MGR solution strategy (solver specific)
    integer separateComponents = false; ///< Apply a separate displacement component (SDC) filter before AMG construction
    string displacementFieldName;       ///< Displacement field name need for SDC filter
  }
  mgr;                                  ///< Multigrid reduction (MGR) parameters

  /// Incomplete factorization parameters
  struct ILU
  {
    integer fill = 0;        ///< Fill level
    real64 threshold = 0.0;  ///< Dropping threshold
  }
  ilu;                       ///< Incomplete factorization parameter struct

  /// Domain decomposition parameters
  struct DD
  {
    integer overlap = 0;   ///< Ghost overlap
  }
  dd;                      ///< Domain decomposition parameter struct
};

ENUM_STRINGS( LinearSolverParameters::SolverType,
              "direct",
              "cg",
              "gmres",
              "fgmres",
              "bicgstab",
              "preconditioner" )

ENUM_STRINGS( LinearSolverParameters::PreconditionerType,
              "none",
              "jacobi",
              "gs",
              "sgs",
              "iluk",
              "ilut",
              "icc",
              "ict",
              "amg",
              "mgr",
              "block" )

ENUM_STRINGS( LinearSolverParameters::Direct::ColPerm,
              "none",
              "MMD_AtplusA",
              "MMD_AtA",
              "colAMD",
              "metis",
              "parmetis" )

ENUM_STRINGS( LinearSolverParameters::Direct::RowPerm,
              "none",
              "mc64" )

ENUM_STRINGS( LinearSolverParameters::MGR::StrategyType,
              "compositionalMultiphaseFVM",
              "compositionalMultiphaseHybridFVM",
              "compositionalMultiphaseReservoir",
              "hydrofracture",
              "lagrangianContactMechanics",
              "singlePhasePoroelastic",
              "hybridSinglePhasePoroelastic" )

} /* namespace geosx */

#endif /*GEOSX_LINEARALGEBRA_UTILITIES_LINEARSOLVERPARAMETERS_HPP_ */
