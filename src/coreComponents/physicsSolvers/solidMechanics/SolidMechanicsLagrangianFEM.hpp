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
 * @file SolidMechanicsLagrangianFEM.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSLAGRANGIANFEM_HPP_
#define GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSLAGRANGIANFEM_HPP_

#include "common/format/EnumStrings.hpp"
#include "common/TimingMacros.hpp"
#include "kernels/SolidMechanicsLagrangianFEMKernels.hpp"
#include "kernels/StressStrainAverageKernels.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "mesh/mpiCommunications/MPI_iCommData.hpp"
#include "linearAlgebra/utilities/EigenSolverParameters.hpp"
#include "physicsSolvers/PhysicsSolverBase.hpp"

#include "physicsSolvers/solidMechanics/SolidMechanicsFields.hpp"

namespace geos
{

/**
 * @class SolidMechanicsLagrangianFEM
 *
 * This class implements a finite element solution to the equations of motion.
 */
class SolidMechanicsLagrangianFEM : public PhysicsSolverBase
{
public:

  /// String used to form the solverName used to register single-physics solvers in CoupledSolver
  static string coupledSolverAttributePrefix() { return "solid"; }

  /**
   * @enum TimeIntegrationOption
   *
   * The options for time integration
   */
  enum class TimeIntegrationOption : integer
  {
    QuasiStatic,      //!< QuasiStatic
    ImplicitDynamic,  //!< ImplicitDynamic
    ExplicitDynamic,  //!< ExplicitDynamic
    Modal             //!< Modal analysis: generalized eigenproblem K x = lambda M x (no time integration)
  };

  /**
   * Constructor
   * @param name The name of the solver instance
   * @param parent the parent group of the solver
   */
  SolidMechanicsLagrangianFEM( const string & name,
                               Group * const parent );

  /**
   * @return The string that may be used to generate a new instance from the PhysicsSolverBase::CatalogInterface::CatalogType
   */
  static string catalogName() { return "SolidMechanicsLagrangianFEM"; }
  /**
   * @copydoc PhysicsSolverBase::getCatalogName()
   */
  string getCatalogName() const override { return catalogName(); }

  virtual void initializePreSubGroups() override;

  virtual void registerDataOnMesh( Group & meshBodies ) override;

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/
  virtual
  real64 solverStep( real64 const & time_n,
                     real64 const & dt,
                     integer const cycleNumber,
                     DomainPartition & domain ) override;

  virtual
  real64 explicitStep( real64 const & time_n,
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
   * Assembles the tangent stiffness K at the current state and the lumped mass M, and solves the generalized
   * eigenproblem K phi = lambda M phi for the modes closest to the shift @p modalShiftFrequency, with the
   * eigensolver selected by @p modalSolverType. Displacement boundary conditions (taken as homogeneous) remove
   * the constrained degrees of freedom. Frequencies, residuals and participation factors are stored in the
   * solver, and the mode shapes in the nodal fields `modeShape1`, `modeShape2`, ...
   */
  real64 modalAnalysisStep( real64 const & time_n,
                            real64 const & dt,
                            integer const cycleNumber,
                            DomainPartition & domain );

  virtual void
  implicitStepSetup( real64 const & time_n,
                     real64 const & dt,
                     DomainPartition & domain ) override;

  virtual void
  setupDofs( DomainPartition const & domain,
             DofManager & dofManager ) const override;

  virtual void
  setupSystem( DomainPartition & domain,
               DofManager & dofManager,
               CRSMatrix< real64, globalIndex > & localMatrix,
               ParallelVector & rhs,
               ParallelVector & solution,
               bool setSparsity = true ) override;

  virtual void
  setSparsityPattern( DomainPartition & domain,
                      DofManager & dofManager,
                      CRSMatrix< real64, globalIndex > & localMatrix,
                      SparsityPattern< globalIndex > & pattern ) override;

  virtual std::unique_ptr< PreconditionerBase< LAInterface > >
  createPreconditioner( DomainPartition & domain ) const override;

  virtual void
  assembleSystem( real64 const time,
                  real64 const dt,
                  DomainPartition & domain,
                  DofManager const & dofManager,
                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                  arrayView1d< real64 > const & localRhs ) override;

  virtual void solveLinearSystem( DofManager const & dofManager,
                                  ParallelMatrix & matrix,
                                  ParallelVector & rhs,
                                  ParallelVector & solution,
                                  integer const cycleNumber,
                                  integer const nonlinearIteration ) override;

  virtual void
  applySystemSolution( DofManager const & dofManager,
                       arrayView1d< real64 const > const & localSolution,
                       real64 const scalingFactor,
                       real64 const dt,
                       DomainPartition & domain ) override;

  virtual void updateState( DomainPartition & domain ) override
  {
    // There should be nothing to update
    GEOS_UNUSED_VAR( domain );
  };

  virtual void applyBoundaryConditions( real64 const time,
                                        real64 const dt,
                                        DomainPartition & domain,
                                        DofManager const & dofManager,
                                        CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                        arrayView1d< real64 > const & localRhs ) override;

  virtual real64
  calculateResidualNorm( real64 const & time_n,
                         real64 const & dt,
                         DomainPartition const & domain,
                         DofManager const & dofManager,
                         arrayView1d< real64 const > const & localRhs ) override;

  virtual void resetStateToBeginningOfStep( DomainPartition & domain ) override;

  virtual void implicitStepComplete( real64 const & time,
                                     real64 const & dt,
                                     DomainPartition & domain ) override;

  /**@}*/


  template< typename TYPE_LIST,
            typename KERNEL_WRAPPER,
            typename ... PARAMS >
  real64 assemblyLaunch( MeshLevel & mesh,
                         DofManager const & dofManager,
                         string_array const & regionNames,
                         string const & materialNamesString,
                         CRSMatrixView< real64, globalIndex const > const & localMatrix,
                         arrayView1d< real64 > const & localRhs,
                         real64 const dt,
                         PARAMS && ... params );

  real64 explicitKernelDispatch( MeshLevel & mesh,
                                 string_array const & targetRegions,
                                 string const & finiteElementName,
                                 real64 const dt,
                                 std::string const & elementListName );

  /**
   * Applies displacement boundary conditions to the system for implicit time integration
   * @param time The time to use for any lookups associated with this BC
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param domain The DomainPartition.
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   * @param solution the solution vector
   */
  void applyDisplacementBCImplicit( real64 const time,
                                    DofManager const & dofManager,
                                    DomainPartition & domain,
                                    CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                    arrayView1d< real64 > const & localRhs );

  void applyTractionBC( real64 const time,
                        DofManager const & dofManager,
                        DomainPartition & domain,
                        arrayView1d< real64 > const & localRhs );

  void applyChomboPressure( DofManager const & dofManager,
                            DomainPartition & domain,
                            arrayView1d< real64 > const & localRhs );


  void applyContactConstraint( DofManager const & dofManager,
                               DomainPartition & domain,
                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs );

  virtual real64
  scalingForSystemSolution( DomainPartition & domain,
                            DofManager const & dofManager,
                            arrayView1d< real64 const > const & localSolution ) override;

  void enableFixedStressPoromechanicsUpdate();

  virtual void saveSequentialIterationState( DomainPartition & domain ) override;

  struct viewKeyStruct : PhysicsSolverBase::viewKeyStruct
  {
    static constexpr char const * newmarkGammaString() { return "newmarkGamma"; }
    static constexpr char const * newmarkBetaString() { return "newmarkBeta"; }
    static constexpr char const * massDampingString() { return "massDamping"; }
    static constexpr char const * stiffnessDampingString() { return "stiffnessDamping"; }
    static constexpr char const * timeIntegrationOptionString() { return "timeIntegrationOption"; }
    static constexpr char const * maxNumResolvesString() { return "maxNumResolves"; }
    static constexpr char const * strainTheoryString() { return "strainTheory"; }
    static constexpr char const * solidMaterialNamesString() { return "solidMaterialNames"; }
    static constexpr char const * contactRelationNameString() { return "contactRelationName"; }
    static constexpr char const * noContactRelationNameString() { return "NOCONTACT"; }
    static constexpr char const * maxForceString() { return "maxForce"; }
    static constexpr char const * elemsAttachedToSendOrReceiveNodesString() { return "elemsAttachedToSendOrReceiveNodes"; }
    static constexpr char const * elemsNotAttachedToSendOrReceiveNodesString() { return "elemsNotAttachedToSendOrReceiveNodes"; }
    static constexpr char const * surfaceGeneratorNameString() { return "surfaceGeneratorName"; }

    static constexpr char const * sendOrReceiveNodesString() { return "sendOrReceiveNodes";}
    static constexpr char const * nonSendOrReceiveNodesString() { return "nonSendOrReceiveNodes";}
    static constexpr char const * targetNodesString() { return "targetNodes";}
    static constexpr char const * forceString() { return "Force";}

    static constexpr char const * contactPenaltyStiffnessString() { return "contactPenaltyStiffness"; }

    static constexpr char const * modalNumModesString() { return "modalNumModes"; }
    static constexpr char const * modalShiftFrequencyString() { return "modalShiftFrequency"; }
    static constexpr char const * modalSolverTypeString() { return "modalSolverType"; }
    static constexpr char const * modalToleranceString() { return "modalTolerance"; }
    static constexpr char const * modalMaxIterationsString() { return "modalMaxIterations"; }
    static constexpr char const * modalSubspaceSizeString() { return "modalSubspaceSize"; }
    static constexpr char const * modalBlockSizeString() { return "modalBlockSize"; }
    static constexpr char const * modalCompletenessCheckString() { return "modalCompletenessCheck"; }
    static constexpr char const * modalSeedString() { return "modalSeed"; }
    static constexpr char const * modalDeflateRigidBodyModesString() { return "modalDeflateRigidBodyModes"; }
    static constexpr char const * modalEigenvaluesString() { return "modalEigenvalues"; }
    static constexpr char const * modalFrequenciesString() { return "modalFrequencies"; }
    static constexpr char const * modalResidualsString() { return "modalResiduals"; }
    static constexpr char const * modalParticipationFactorsString() { return "modalParticipationFactors"; }

  };

  SortedArray< localIndex > & getElemsAttachedToSendOrReceiveNodes( ElementSubRegionBase & subRegion )
  {
    return subRegion.getReference< SortedArray< localIndex > >( viewKeyStruct::elemsAttachedToSendOrReceiveNodesString() );
  }

  SortedArray< localIndex > & getElemsNotAttachedToSendOrReceiveNodes( ElementSubRegionBase & subRegion )
  {
    return subRegion.getReference< SortedArray< localIndex > >( viewKeyStruct::elemsNotAttachedToSendOrReceiveNodesString() );
  }

  real64 & getMaxForce() { return m_maxForce; }
  real64 const & getMaxForce() const { return m_maxForce; }

  void computeRigidBodyModes( DomainPartition & domain ) const;

  arrayView1d< ParallelVector > const & getRigidBodyModes( DomainPartition & domain ) const
  {
    computeRigidBodyModes( domain );
    return m_rigidBodyModes;
  }

  /*
   * @brief Utility function to set the stress initialization flag
   * @param[in] performStressInitialization true if the solver has to initialize stress, false otherwise
   */
  void setStressInitialization( bool const performStressInitialization )
  {
    m_performStressInitialization = performStressInitialization;
  }

  TimeIntegrationOption timeIntegrationOption() const { return m_timeIntegrationOption; }

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

  void initializeMass( MeshLevel & mesh, CellElementSubRegion & subRegion );

  /**
   * @brief Build the lumped mass vector and the mask of free degrees of freedom of the modal analysis.
   * @param[in] time time at which the displacement boundary conditions are evaluated
   * @param[in] domain the domain
   * @param[out] freeMask vector with 1 on free degrees of freedom and 0 on constrained ones
   * @param[out] massDiag lumped mass, zero on constrained degrees of freedom
   * @return the global number of constrained degrees of freedom
   */
  globalIndex computeModalDiagonals( real64 const time,
                                     DomainPartition & domain,
                                     ParallelVector & freeMask,
                                     ParallelVector & massDiag );

  virtual void initializePostInitialConditionsPreSubGroups() override;

  virtual void setConstitutiveNamesCallSuper( ElementSubRegionBase & subRegion ) const override;

  real64 m_newmarkGamma;
  real64 m_newmarkBeta;
  real64 m_massDamping;
  real64 m_stiffnessDamping;
  TimeIntegrationOption m_timeIntegrationOption;
  real64 m_maxForce = 0.0;
  integer m_maxNumResolves;
  integer m_strainTheory;

  /// Flag to indicate that the solver is running in fixed stress (sequential) mode
  bool m_isFixedStressPoromechanicsUpdate;
  /// Flag to indicate that the solver is going to perform stress initialization
  bool m_performStressInitialization;

  /// Rigid body modes; TODO remove mutable hack
  mutable array1d< ParallelVector > m_rigidBodyModes;

  real64 m_contactPenaltyStiffness;

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
  /// True while the modal analysis assembles its operators: the Modal option is only valid for a standalone
  /// solver that runs modalAnalysisStep(), not when the solver is driven by a coupled solver
  bool m_isModalAssembly = false;

  /// Eigenvalues lambda = omega^2 of the last modal analysis
  array1d< real64 > m_modalEigenvalues;
  /// Signed frequencies (Hz) of the last modal analysis
  array1d< real64 > m_modalFrequencies;
  /// Relative residuals of the last modal analysis
  array1d< real64 > m_modalResiduals;
  /// Participation factors of the last modal analysis: ( mode, direction )
  array2d< real64 > m_modalParticipationFactors;

private:

  string m_contactRelationName;

  PhysicsSolverBase *m_surfaceGenerator;
  string m_surfaceGeneratorName;
};

ENUM_STRINGS( SolidMechanicsLagrangianFEM::TimeIntegrationOption,
              "QuasiStatic",
              "ImplicitDynamic",
              "ExplicitDynamic",
              "Modal" );

//**********************************************************************************************************************
//**********************************************************************************************************************
//**********************************************************************************************************************


template< typename TYPE_LIST,
          typename KERNEL_WRAPPER,
          typename ... PARAMS >
real64 SolidMechanicsLagrangianFEM::assemblyLaunch( MeshLevel & mesh,
                                                    DofManager const & dofManager,
                                                    string_array const & regionNames,
                                                    string const & materialNamesString,
                                                    CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                    arrayView1d< real64 > const & localRhs,
                                                    real64 const dt,
                                                    PARAMS && ... params )
{
  GEOS_MARK_FUNCTION;

  NodeManager const & nodeManager = mesh.getNodeManager();

  string const dofKey = dofManager.getKey( fields::solidMechanics::totalDisplacement::key() );
  arrayView1d< globalIndex const > const & dofNumber = nodeManager.getReference< globalIndex_array >( dofKey );

  real64 const gravityVectorData[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( gravityVector() );

  KERNEL_WRAPPER kernelWrapper( dofNumber,
                                dofManager.rankOffset(),
                                localMatrix,
                                localRhs,
                                dt,
                                gravityVectorData,
                                std::forward< PARAMS >( params )... );

  return finiteElement::
           regionBasedKernelApplication< parallelDevicePolicy< >,
                                         TYPE_LIST >( mesh,
                                                      regionNames,
                                                      this->getDiscretizationName(),
                                                      materialNamesString,
                                                      kernelWrapper );

}

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSLAGRANGIANFEM_HPP_ */
