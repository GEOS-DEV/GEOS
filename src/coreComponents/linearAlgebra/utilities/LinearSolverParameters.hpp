/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
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
    richardson,    ///< Richardson iteration
    preconditioner ///< Preconditioner only
  };

  /**
   * @brief Preconditioner type.
   */
  enum class PreconditionerType : integer
  {
    none,      ///< No preconditioner
    jacobi,    ///< Jacobi smoothing
    l1jacobi,  ///< l1-Jacobi smoothing
    fgs,       ///< Gauss-Seidel smoothing (forward sweep)
    sgs,       ///< Symmetric Gauss-Seidel smoothing
    l1sgs,     ///< l1-Symmetric Gauss-Seidel smoothing
    chebyshev, ///< Chebyshev polynomial smoothing
    ilu,       ///< Incomplete LU with k-level of fill
    ilut,      ///< Incomplete LU with thresholding
    ic,        ///< Incomplete Cholesky
    ict,       ///< Incomplete Cholesky with thresholding
    amg,       ///< Algebraic Multigrid
    mgr,       ///< Multigrid reduction (Hypre only)
    block,     ///< Block preconditioner
    direct,    ///< Direct solver as preconditioner
    bgs,       ///< Gauss-Seidel smoothing (backward sweep)
    multiscale,///< Multiscale preconditioner
  };

  integer logLevel = 0;     ///< Output level [0=none, 1=basic, 2=everything]
  integer dofsPerNode = 1;  ///< Dofs per node (or support location) for non-scalar problems
  bool isSymmetric = false; ///< Whether input matrix is symmetric (may affect choice of scheme)
  integer stopIfError = 1;  ///< Whether to stop the simulation if the linear solver reports an error

  SolverType solverType = SolverType::direct;          ///< Solver type
  PreconditionerType preconditionerType = PreconditionerType::ilu;  ///< Preconditioner type

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

  /// Relaxation/stationary iteration parameters (Richardson, damped Jacobi, etc.)
  struct Relaxation
  {
    real64 weight = 2.0 / 3.0;        ///< Relaxation weight (omega) for stationary iterations
  } relaxation;

  /// Chebyshev iteration/smoothing parameters
  struct Chebyshev
  {
    integer order = 2;                ///< Chebyshev order
    integer eigNumIter = 10;          ///< Number of eigenvalue estimation CG iterations
  } chebyshev;

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
    /// AMG cycle type
    enum class CycleType : integer
    {
      V, ///< V-cycle
      W, ///< W-cycle
    };

    /// AMG pre/post smoothing option
    enum class PreOrPost : integer
    {
      pre,  ///< pre-smoothing only
      post, ///< post-smoothing only
      both  ///< pre- and post-smoothing
    };

    /// AMG smoother type
    enum class SmootherType : integer
    {
      default_,  ///< Use LAI's default option
      jacobi,    ///< Jacobi smoothing
      l1jacobi,  ///< l1-Jacobi smoothing
      fgs,       ///< Gauss-Seidel smoothing (forward sweep)
      bgs,       ///< Gauss-Seidel smoothing (backward sweep)
      sgs,       ///< Symmetric Gauss-Seidel smoothing
      l1sgs,     ///< l1-Symmetric Gauss-Seidel smoothing
      chebyshev, ///< Chebyshev polynomial smoothing
      ilu,       ///< Incomplete LU
      ilut,      ///< Incomplete LU with thresholding
      ic,        ///< Incomplete Cholesky
      ict,       ///< Incomplete Cholesky with thresholding
    };

    /// AMG coarse solver type
    enum class CoarseType : integer
    {
      default_,  ///< Use LAI's default option
      jacobi,    ///< Jacobi
      l1jacobi,  ///< l1-Jacobi
      fgs,       ///< Gauss-Seidel (forward sweep)
      sgs,       ///< Symmetric Gauss-Seidel
      l1sgs,     ///< l1-Symmetric Gauss-Seidel
      chebyshev, ///< Chebyshev polynomial
      direct,    ///< Direct solver as preconditioner
      bgs,       ///< Gauss-Seidel smoothing (backward sweep)
    };

    /// Null space type
    enum class NullSpaceType : integer
    {
      constantModes,  ///< Constant modes
      rigidBodyModes, ///< Rigid body modes
    };

    integer maxLevels = 20;                         ///< Maximum number of coarsening levels
    integer numCycles = 1;                          ///< Number of multigrid cycles
    CycleType cycleType = CycleType::V;             ///< AMG cycle type
    SmootherType smootherType = SmootherType::fgs;  ///< Smoother type
    CoarseType coarseType = CoarseType::direct;     ///< Coarse-level solver/smoother
    string coarseningType = "HMIS";                 ///< Coarsening algorithm
    integer interpolationType = 6;                  ///< Coarsening algorithm
    integer numSweeps = 2;                          ///< Number of smoother sweeps
    integer numFunctions = 1;                       ///< Number of amg functions
    integer aggresiveNumLevels = 0;                 ///< Number of levels for aggressive coarsening.
    integer maxCoarseSize = 9;                      ///< Threshold for coarse grid size
    PreOrPost preOrPostSmoothing = PreOrPost::both; ///< Pre and/or post smoothing
    real64 threshold = 0.0;                         ///< Threshold for "strong connections" (for classical and smoothed-aggregation AMG)
    integer separateComponents = false;             ///< Apply a separate component filter before AMG construction
    NullSpaceType nullSpaceType = NullSpaceType::constantModes; ///< Null space type [constantModes,rigidBodyModes]
  }
  amg;                                              ///< Algebraic Multigrid (AMG) parameters

  /// Multigrid reduction parameters
  struct MGR
  {
    /**
     * @brief MGR available strategies
     */
    enum class StrategyType : integer
    {
      invalid,                                   ///< default value, to ensure solver sets something
      singlePhaseReservoirFVM,                   ///< finite volume single-phase flow with wells
      singlePhaseHybridFVM,                      ///< hybrid finite volume single-phase flow
      singlePhaseReservoirHybridFVM,             ///< hybrid finite volume single-phase flow with wells
      singlePhasePoromechanics,                  ///< single phase poromechanics with finite volume single phase flow
      hybridSinglePhasePoromechanics,            ///< single phase poromechanics with hybrid finite volume single phase flow
      singlePhasePoromechanicsEmbeddedFractures, ///< single phase poromechanics with finite volume single phase flow and embedded fractures
      compositionalMultiphaseFVM,                ///< finite volume compositional multiphase flow
      compositionalMultiphaseHybridFVM,          ///< hybrid finite volume compositional multiphase flow
      compositionalMultiphaseReservoirFVM,       ///< finite volume compositional multiphase flow with wells
      compositionalMultiphaseReservoirHybridFVM, ///< hybrid finite volume compositional multiphase flow with wells
      thermalCompositionalMultiphaseFVM,         ///< finite volume thermal compositional multiphase flow
      multiphasePoromechanics,                   ///< multiphase poromechanics with finite volume compositional multiphase flow
      hydrofracture,                             ///< hydrofracture
      lagrangianContactMechanics,                ///< Lagrangian contact mechanics
      solidMechanicsEmbeddedFractures            ///< Embedded fractures mechanics
    };

    StrategyType strategy = StrategyType::invalid; ///< Predefined MGR solution strategy (solver specific)
    integer separateComponents = false;            ///< Apply a separate displacement component (SDC) filter before AMG construction
    string displacementFieldName;                  ///< Displacement field name need for SDC filter
  }
  mgr;                                             ///< Multigrid reduction (MGR) parameters

  /// Incomplete factorization parameters
  struct IFact
  {
    integer fill = 0;        ///< Fill level
    real64 threshold = 0.0;  ///< Dropping threshold
  }
  ifact;                       ///< Incomplete factorization parameter struct

  /// Domain decomposition parameters
  struct DD
  {
    integer overlap = 0;   ///< Ghost overlap
  }
  dd;                      ///< Domain decomposition parameter struct

  /// Block preconditioner parameters
  struct Block
  {
    /// Shape of the block preconditioner
    enum class Shape
    {
      Diagonal,            //!< (D)^{-1}
      UpperTriangular,     //!< (DU)^{-1}
      LowerTriangular,     //!< (LD)^{-1}
      LowerUpperTriangular //!< (LDU)^{-1}
    };

    /// Type of Schur complement approximation used
    enum class SchurType
    {
      None,                  //!< No Schur complement - just block-GS/block-Jacobi preconditioner
      FirstBlockDiagonal,    //!< Approximate first block with its diagonal
      RowsumDiagonalProbing, //!< Rowsum-preserving diagonal approximation constructed with probing
      FirstBlockUserDefined  //!< User defined preconditioner for the first block
    };

    /// Type of block row scaling to apply
    enum class Scaling
    {
      None,          //!< No scaling
      FrobeniusNorm, //!< Equilibrate Frobenius norm of the diagonal blocks
      UserProvided   //!< User-provided scaling
    };

    Shape shape = Shape::UpperTriangular;
    SchurType schurType = SchurType::RowsumDiagonalProbing;
    Scaling scaling = Scaling::FrobeniusNorm;

    array1d< LinearSolverParameters const * > subParams; ///< Pointers to parameters for subproblems
  } block;

  struct Multiscale
  {
    enum class BasisType
    {
      msrsb,   ///< Restricted Smoothing Basis Multiscale
    };

    // Core parameters
    BasisType basisType = BasisType::msrsb;                     ///< Type of basis functions
    string fieldName;                                           ///< DofManager field name, populated by the physics solver
    integer maxLevels = 5;                                      ///< Limit on total number of grid levels
    integer galerkin = 1;                                       ///< Whether to use Galerkin definition R = P^T (otherwise R = P_0^T)
    real64 droptol = 0.0;                                       ///< Dropping tolerance for coarse matrix values (relative to row max)
    PreconditionerType coarseType = PreconditionerType::direct; ///< Coarse solver type
    integer separateComponents = false;                         ///< Apply a separate component filter before multiscale
    array1d< string > boundarySets;                             ///< List of boundary node set names (needed for coarse node detection)

    // Debugging/user-display settings
    string label;                                               ///< User-displayed label of the scheme
    integer debugLevel = 0;                                     ///< Flag for enabling addition debugging output

    struct Smoother
    {
      PreconditionerType type = PreconditionerType::sgs; ///< Smoother type
      AMG::PreOrPost preOrPost = AMG::PreOrPost::both;   ///< Pre and/or post smoothing [pre,post,both]
      integer numSweeps = 1;                             ///< Number of smoother sweeps
    } smoother;

    struct Coupled
    {
      integer useBlockSmoother = 1; ///< Whether to use block smoother
    } coupled;

    struct Coarsening
    {
      enum class PartitionType
      {
        graph,          ///< Graph-based
        cartesian,      ///< Cartesian (only with internal mesh)
        semistructured, ///< For 2.5D (extruded) meshes
      };

      PartitionType partitionType = PartitionType::graph; ///< Partitioning/coarsening method
      array1d< real64 > ratio;                            ///< Coarsening ratio (cell-based), total or per-dimension
      globalIndex maxCoarseDof = 0;                       ///< Limit of coarsening globally (trims the grid hierarchy)

      struct Structured
      {
        integer semicoarsening = 0; ///< Use semicoarsenining in Z-direction
      } structured;

      struct Graph
      {
        enum class Method
        {
          metis,
          scotch,
        };

        Method method = Method::metis; ///< Graph partitioning method (library)
        integer minCommonNodes = 3;    ///< Min number of common nodes when constructing a cell connectivity graph
        integer preserveRegions = 0;   ///< Attempt to keep cells from the same region in one aggregate
        integer matrixWeights = 0;     ///< If >0, specifies matrix weight multiplier
        CRSMatrix< real64, localIndex > const * localMatrix{}; ///< Local matrix to use for weights

        struct Metis
        {
          enum class Method
          {
            kway,
            recursive
          };

          Method method = Method::kway;    ///< METIS partitioning algorithm
          integer ufactor = 30;            ///< METIS UFACTOR option (allowed partition imbalance)
          integer seed = 2020;             ///< Random number generator seed
        } metis;

        struct Scotch
        {
          // TODO
        } scotch;

      } graph;

    } coarsening;

    struct MsRSB
    {
      enum class SupportType : integer
      {
        layers,
        matching
      };

      SupportType support     = SupportType::matching; ///< algorithm used to construct support regions
      integer numLayers       = 3;          ///< Number of layers to grow (for support = layers)
      integer maxIter         = 100;        ///< Max number of smoothing iterations
      real64 tolerance        = 1e-3;       ///< Smoothing iteration convergence tolerance
      real64 relaxation       = 2.0 / 3.0;  ///< Relaxation parameter for Jacobi smoothing
      integer checkFrequency  = 10;         ///< Convergence check frequency (num smoothing iterations)
      integer updateFrequency = 10;         ///< Coarse operator update frequency (num smoothing iterations)
    } msrsb;
  }
  multiscale;
};

/// Declare strings associated with enumeration values.
ENUM_STRINGS( LinearSolverParameters::SolverType,
              "direct",
              "cg",
              "gmres",
              "fgmres",
              "bicgstab",
              "richardson",
              "preconditioner" );

/// Declare strings associated with enumeration values.
ENUM_STRINGS( LinearSolverParameters::PreconditionerType,
              "none",
              "jacobi",
              "l1jacobi",
              "fgs",
              "sgs",
              "l1sgs",
              "chebyshev",
              "ilu",
              "ilut",
              "ic",
              "ict",
              "amg",
              "mgr",
              "block",
              "direct",
              "bgs",
              "multiscale" );

/// Declare strings associated with enumeration values.
ENUM_STRINGS( LinearSolverParameters::Direct::ColPerm,
              "none",
              "MMD_AtplusA",
              "MMD_AtA",
              "colAMD",
              "metis",
              "parmetis" );

/// Declare strings associated with enumeration values.
ENUM_STRINGS( LinearSolverParameters::Direct::RowPerm,
              "none",
              "mc64" );

/// Declare strings associated with enumeration values.
ENUM_STRINGS( LinearSolverParameters::MGR::StrategyType,
              "invalid",
              "singlePhaseReservoirFVM",
              "singlePhaseHybridFVM",
              "singlePhaseReservoirHybridFVM",
              "singlePhasePoromechanics",
              "hybridSinglePhasePoromechanics",
              "singlePhasePoromechanicsEmbeddedFractures",
              "compositionalMultiphaseFVM",
              "compositionalMultiphaseHybridFVM",
              "compositionalMultiphaseReservoirFVM",
              "compositionalMultiphaseReservoirHybridFVM",
              "thermalCompositionalMultiphaseFVM",
              "multiphasePoromechanics",
              "hydrofracture",
              "lagrangianContactMechanics",
              "solidMechanicsEmbeddedFractures" );

/// Declare strings associated with enumeration values.
ENUM_STRINGS( LinearSolverParameters::AMG::CycleType,
              "V",
              "W" );

/// Declare strings associated with enumeration values.
ENUM_STRINGS( LinearSolverParameters::AMG::PreOrPost,
              "pre",
              "post",
              "both" );

/// Declare strings associated with enumeration values.
ENUM_STRINGS( LinearSolverParameters::AMG::SmootherType,
              "default",
              "jacobi",
              "l1jacobi",
              "fgs",
              "bgs",
              "sgs",
              "l1sgs",
              "chebyshev",
              "ilu",
              "ilut",
              "ic",
              "ict" );

/// Declare strings associated with enumeration values.
ENUM_STRINGS( LinearSolverParameters::AMG::CoarseType,
              "default",
              "jacobi",
              "l1jacobi",
              "fgs",
              "sgs",
              "l1sgs",
              "chebyshev",
              "direct",
              "bgs", );

/// Declare strings associated with enumeration values.
ENUM_STRINGS( LinearSolverParameters::AMG::NullSpaceType,
              "constantModes",
              "rigidBodyModes" );

/// Declare strings associated with enumeration values.
ENUM_STRINGS( LinearSolverParameters::Block::Shape,
              "D",
              "DU",
              "LD",
              "LDU" );

/// Declare strings associated with enumeration values.
ENUM_STRINGS( LinearSolverParameters::Block::SchurType,
              "none",
              "diagonal",
              "probing",
              "user" );

/// Declare strings associated with enumeration values.
ENUM_STRINGS( LinearSolverParameters::Block::Scaling,
              "none",
              "frobenius",
              "user" );

/// Declare strings associated with enumeration values.
ENUM_STRINGS( LinearSolverParameters::Multiscale::BasisType,
              "msrsb" );

/// Declare strings associated with enumeration values.
ENUM_STRINGS( LinearSolverParameters::Multiscale::Coarsening::PartitionType,
              "graph",
              "cartesian",
              "semistructured" );

/// Declare strings associated with enumeration values.
ENUM_STRINGS( LinearSolverParameters::Multiscale::Coarsening::Graph::Method,
              "metis",
              "scotch" );

/// Declare strings associated with enumeration values.
ENUM_STRINGS( LinearSolverParameters::Multiscale::Coarsening::Graph::Metis::Method,
              "kway",
              "recursive" );

/// Declare strings associated with enumeration values.
ENUM_STRINGS( LinearSolverParameters::Multiscale::MsRSB::SupportType,
              "layers",
              "matching" );

} /* namespace geosx */

#endif /*GEOSX_LINEARALGEBRA_UTILITIES_LINEARSOLVERPARAMETERS_HPP_ */
