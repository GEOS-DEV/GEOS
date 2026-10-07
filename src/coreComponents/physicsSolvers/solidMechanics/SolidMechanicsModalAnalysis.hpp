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
 * @file SolidMechanicsModalAnalysis.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSMODALANALYSIS_HPP_
#define GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSMODALANALYSIS_HPP_

#include "common/format/EnumStrings.hpp"
#include "linearAlgebra/utilities/EigenSolverParameters.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"

namespace geos
{

/**
 * @class SolidMechanicsModalAnalysis
 * @brief Vibration modes of a structure: the generalized eigenproblem K phi = lambda M phi.
 *
 * The solver assembles the tangent stiffness K of the quasi-static solid mechanics solver at the current state
 * and a mass M, and computes the modes closest to a spectral shift with the eigensolvers of the linear algebra
 * layer (shift-and-invert Arnoldi, or LOBPCG). It does not integrate in time: one execution of the solver
 * computes the modes. Displacement boundary conditions are taken as homogeneous and remove the constrained
 * degrees of freedom. A free structure is the main use case: its six rigid-body modes have a zero eigenvalue.
 */
class SolidMechanicsModalAnalysis : public SolidMechanicsLagrangianFEM
{
public:

  /**
   * @enum MassType
   * @brief Discretization of the mass matrix.
   */
  enum class MassType : integer
  {
    lumped,      //!< Row-sum lumped mass (diagonal)
    consistent   //!< Exact consistent mass, available for first-order tetrahedra only
  };

  /**
   * @brief Constructor.
   * @param name the name of the solver instance
   * @param parent the parent group of the solver
   */
  SolidMechanicsModalAnalysis( string const & name,
                               Group * const parent );

  /**
   * @return The string that may be used to generate a new instance from the PhysicsSolverBase::CatalogInterface::CatalogType
   */
  static string catalogName() { return "SolidMechanicsModalAnalysis"; }

  /**
   * @copydoc PhysicsSolverBase::getCatalogName()
   */
  string getCatalogName() const override { return catalogName(); }

  virtual void registerDataOnMesh( Group & meshBodies ) override;

  virtual real64 solverStep( real64 const & time_n,
                             real64 const & dt,
                             integer const cycleNumber,
                             DomainPartition & domain ) override;

  /**
   * @brief Compute the vibration modes of the structure.
   * @param time_n time at the beginning of the step
   * @param dt time step (unused except for the evaluation of time-dependent coefficients)
   * @param cycleNumber cycle number
   * @param domain the domain
   * @return the time step (unchanged)
   *
   * Assembles the tangent stiffness K and the selected modal mass M, and solves the generalized
   * eigenproblem K phi = lambda M phi for the modes closest to the shift @p modalShiftFrequency, with the
   * eigensolver selected by @p modalSolverType. Displacement boundary conditions (taken as homogeneous) remove
   * the constrained degrees of freedom. Frequencies, residuals and participation factors are stored in the
   * solver, and the mode shapes in the nodal fields `modeShape1`, `modeShape2`, ...
   * @note This function is public because it launches device kernels (extended lambdas of CUDA cannot be
   *       defined in protected or private member functions).
   */
  real64 modalAnalysisStep( real64 const & time_n,
                            real64 const & dt,
                            integer const cycleNumber,
                            DomainPartition & domain );

  /**
   * @brief Build the lumped mass vector and the mask of free degrees of freedom of the modal analysis.
   * @param[in] time time at which the displacement boundary conditions are evaluated
   * @param[in] domain the domain
   * @param[out] freeMask vector with 1 on free degrees of freedom and 0 on constrained ones
   * @param[out] massDiag lumped mass, zero on constrained degrees of freedom
   * @return the global number of constrained degrees of freedom
   * @note This function is public for the same reason as modalAnalysisStep().
   */
  globalIndex computeModalDiagonals( real64 const time,
                                     DomainPartition & domain,
                                     ParallelVector & freeMask,
                                     ParallelVector & massDiag );

  /**
   * @brief Assemble the exact consistent mass for first-order tetrahedra, using the stiffness sparsity.
   * @param[in] domain the domain
   * @param[out] massMatrix the mass matrix
   * @note This function is public for the same reason as modalAnalysisStep().
   */
  void assembleModalConsistentMass( DomainPartition & domain, ParallelMatrix & massMatrix );

  /**
   * @struct viewKeyStruct
   * @brief Keys of the input attributes and of the results.
   */
  struct viewKeyStruct : SolidMechanicsLagrangianFEM::viewKeyStruct
  {
    /// @return key of the number of modes
    static constexpr char const * modalNumModesString() { return "modalNumModes"; }
    /// @return key of the shift frequency
    static constexpr char const * modalShiftFrequencyString() { return "modalShiftFrequency"; }
    /// @return key of the eigensolver type
    static constexpr char const * modalSolverTypeString() { return "modalSolverType"; }
    /// @return key of the eigensolver tolerance
    static constexpr char const * modalToleranceString() { return "modalTolerance"; }
    /// @return key of the maximum number of restarts or iterations
    static constexpr char const * modalMaxIterationsString() { return "modalMaxIterations"; }
    /// @return key of the subspace size
    static constexpr char const * modalSubspaceSizeString() { return "modalSubspaceSize"; }
    /// @return key of the block size
    static constexpr char const * modalBlockSizeString() { return "modalBlockSize"; }
    /// @return key of the completeness check flag
    static constexpr char const * modalCompletenessCheckString() { return "modalCompletenessCheck"; }
    /// @return key of the seed of the starting vectors
    static constexpr char const * modalSeedString() { return "modalSeed"; }
    /// @return key of the rigid-body mode deflation flag
    static constexpr char const * modalDeflateRigidBodyModesString() { return "modalDeflateRigidBodyModes"; }
    /// @return key of the mass type
    static constexpr char const * modalMassTypeString() { return "modalMassType"; }
    /// @return key of the free-body verification flag
    static constexpr char const * modalVerifyFreeBodyString() { return "modalVerifyFreeBody"; }
    /// @return key of the eigenvalues
    static constexpr char const * modalEigenvaluesString() { return "modalEigenvalues"; }
    /// @return key of the frequencies
    static constexpr char const * modalFrequenciesString() { return "modalFrequencies"; }
    /// @return key of the residuals
    static constexpr char const * modalResidualsString() { return "modalResiduals"; }
    /// @return key of the participation factors
    static constexpr char const * modalParticipationFactorsString() { return "modalParticipationFactors"; }
  };

  /**
   * @brief Name of the nodal field holding a mode shape.
   * @param mode mode number, starting at 1
   * @return the field name
   */
  static string modeShapeFieldName( integer const mode ) { return GEOS_FMT( "modeShape{}", mode ); }

  /// @return eigenvalues lambda = omega^2 of the last modal analysis, ascending
  arrayView1d< real64 const > modalEigenvalues() const { return m_modalEigenvalues.toViewConst(); }

  /// @return signed frequencies sign(lambda) sqrt(|lambda|) / (2 pi) of the last modal analysis
  arrayView1d< real64 const > modalFrequencies() const { return m_modalFrequencies.toViewConst(); }

  /// @return relative residuals of the eigenpairs of the last modal analysis
  arrayView1d< real64 const > modalResiduals() const { return m_modalResiduals.toViewConst(); }

  /// @return participation factors (mode, direction) of the last modal analysis
  arrayView2d< real64 const > modalParticipationFactors() const { return m_modalParticipationFactors.toViewConst(); }

protected:

  virtual void postInputInitialization() override;

private:

  /// Number of modes requested by the modal analysis
  integer m_modalNumModes;
  /// Spectral shift of the modal analysis, in Hz (signed: the shift in eigenvalue units is sign(f) (2 pi f)^2)
  real64 m_modalShiftFrequency;
  /// Eigensolver used by the modal analysis
  EigenSolverParameters::SolverType m_modalSolverType;
  /// Convergence tolerance of the eigensolver
  real64 m_modalTolerance;
  /// Maximum number of eigensolver restarts/iterations
  integer m_modalMaxIterations;
  /// Maximum dimension of the Krylov basis (0 = default)
  integer m_modalSubspaceSize;
  /// Block size of the Arnoldi eigensolver
  integer m_modalBlockSize;
  /// Whether the Arnoldi eigensolver verifies that no repeated eigenvalue copy was missed
  integer m_modalCompletenessCheck;
  /// Seed of the starting vectors of the eigensolver
  integer m_modalSeed;
  /// Whether the six rigid-body modes of a free structure are deflated from the eigensolve
  integer m_modalDeflateRigidBodyModes;
  /// Mass discretization
  MassType m_modalMassType;
  /// Whether the free-body checks run after the eigensolve
  integer m_modalVerifyFreeBody;

  /// Eigenvalues lambda = omega^2 of the last modal analysis
  array1d< real64 > m_modalEigenvalues;
  /// Signed frequencies (Hz) of the last modal analysis
  array1d< real64 > m_modalFrequencies;
  /// Relative residuals of the last modal analysis
  array1d< real64 > m_modalResiduals;
  /// Participation factors of the last modal analysis: ( mode, direction )
  array2d< real64 > m_modalParticipationFactors;
};

ENUM_STRINGS( SolidMechanicsModalAnalysis::MassType,
              "lumped",
              "consistent" );

} // namespace geos

#endif /* GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSMODALANALYSIS_HPP_ */
