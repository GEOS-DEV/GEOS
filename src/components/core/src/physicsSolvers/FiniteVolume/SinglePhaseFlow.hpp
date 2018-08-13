/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file SinglePhaseFlow_TPFA.hpp
 */

#ifndef SINGLE_PHASE_FLOW_TPFA_HPP_
#define SINGLE_PHASE_FLOW_TPFA_HPP_

#include "physicsSolvers/SolverBase.hpp"

class Epetra_FECrsGraph;

namespace geosx
{
namespace dataRepository
{
class ManagedGroup;
}
class BoundaryConditionBase;
class FiniteElementBase;
class DomainPartition;

/**
 * @class SinglePhaseFlow_TPFA
 *
 * class to perform a single phase finite volume solve.
 */
class SinglePhaseFlow : public SolverBase
{
public:
  /**
   * @brief main constructor for NodeManager Objects
   * @param name the name of this instantiation of NodeManager in the repository
   * @param parent the parent group of this instantiation of NodeManager
   */
  SinglePhaseFlow( const std::string& name,
                        ManagedGroup * const parent );


  /// deleted default constructor
  SinglePhaseFlow() = delete;

  /// deleted copy constructor
  SinglePhaseFlow( SinglePhaseFlow const & ) = delete;

  /// default move constructor
  SinglePhaseFlow( SinglePhaseFlow && ) = default;

  /// deleted assignment operator
  SinglePhaseFlow & operator=( SinglePhaseFlow const & ) = delete;

  /// deleted move operator
  SinglePhaseFlow & operator=( SinglePhaseFlow && ) = delete;

  /**
   * @brief default destructor
   */
  virtual ~SinglePhaseFlow() override = default;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string CatalogName() { return "SinglePhaseFlow"; }


  virtual void FillDocumentationNode() override final;

  virtual void FillOtherDocumentationNodes( dataRepository::ManagedGroup * const group ) override final;

  virtual void FinalInitialization( dataRepository::ManagedGroup * const problemManager ) override final;

  virtual real64 SolverStep( real64 const& time_n,
                             real64 const& dt,
                             integer const cycleNumber,
                             DomainPartition * domain ) override;

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/

  virtual void ImplicitStepSetup( real64 const& time_n,
                              real64 const& dt,
                              DomainPartition * const domain,
                              systemSolverInterface::EpetraBlockSystem * const blockSystem ) override;


  virtual void AssembleSystem( DomainPartition * const domain,
                               systemSolverInterface::EpetraBlockSystem * const blockSystem,
                               real64 const time,
                               real64 const dt ) override;

  virtual void ApplyBoundaryConditions( DomainPartition * const domain,
                                        systemSolverInterface::EpetraBlockSystem * const blockSystem,
                                        real64 const time,
                                        real64 const dt ) override;

  virtual real64
  CalculateResidualNorm(systemSolverInterface::EpetraBlockSystem const *const blockSystem,
                        DomainPartition *const domain) override;

  virtual void SolveSystem( systemSolverInterface::EpetraBlockSystem * const blockSystem,
                            SystemSolverParameters const * const params ) override;

  virtual void
  ApplySystemSolution( systemSolverInterface::EpetraBlockSystem const * const blockSystem,
                       real64 const scalingFactor,
                       DomainPartition * const domain ) override;

  virtual void ResetStateToBeginningOfStep( DomainPartition * const domain ) override;

  virtual  void ImplicitStepComplete( real64 const & time,
                                      real64 const & dt,
                                      DomainPartition * const domain ) override;
  /**@}*/

  void SetupSystem ( DomainPartition * const domain,
                     systemSolverInterface::EpetraBlockSystem * const blockSystem );

  /**
   * @brief set the sparsity pattern for the linear system
   * @param domain the domain partition
   * @param sparsity the sparsity pattern matrix
   */
  void SetSparsityPattern( DomainPartition const * const domain,
                           Epetra_FECrsGraph * const sparsity );

  /**
   * @brief sets the dof indices for this solver
   * @param meshLevel the mesh object (single level only)
   * @param numLocalRows the number of local rows on this partition
   * @param numGlobalRows the number of global rows in the problem
   * @param localIndices unused TODO delete
   * @param offset the DOF offset for this solver in the case of a non-block system
   *
   * This function sets the number of global rows, and sets the dof numbers for
   * this solver. dof numbers are referred to trilinosIndices currently.
   */
  void SetNumRowsAndTrilinosIndices( MeshLevel * const meshLevel,
                                     localIndex & numLocalRows,
                                     globalIndex & numGlobalRows,
                                     localIndex_array& localIndices,
                                     localIndex offset );



  /**
   * @brief Function to perform the Application of Dirichlet type BC's
   * @param object the target ObjectManager for the application of the BC.
   * @param time current time
   * @param blockSystem the entire block system
   */
  void ApplyDirichletBC_implicit( DomainPartition * object,
                                  real64 const time, real64 const dt,
                                  systemSolverInterface::EpetraBlockSystem * const blockSystem);

  /**
   * @brief Function to perform the application of Dirichlet BCs on faces
   * @param domain the domain
   * @param time current time
   * @param blockSystem the entire block system
   */
  void ApplyFaceDirichletBC_implicit(DomainPartition * domain,
                                     real64 const time, real64 const dt,
                                     systemSolverInterface::EpetraBlockSystem * const blockSystem);

  /**
   * @enum an enum to lay out the time integration options.
   */
  enum class timeIntegrationOption
  {
    SteadyState,      //!< SteadyState
    ImplicitTransient,//!< ImplicitTransient
    ExplicitTransient //!< ExplicitTransient
  };

  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    constexpr static auto blockLocalDofNumberString = "blockLocalDofNumber_SinglePhaseFlow";

    constexpr static auto fluidPressureString = "fluidPressure";
    constexpr static auto deltaFluidPressureString = "deltaFluidPressure";

    constexpr static auto fluidDensityString = "fluidDensity";
    constexpr static auto deltaFluidDensityString = "deltaFluidDensity";

    constexpr static auto fluidViscosityString = "fluidViscosity";
    constexpr static auto deltaFluidViscosityString = "deltaFluidViscosity";

    constexpr static auto porosityString = "porosity";
    constexpr static auto deltaPorosityString = "deltaPorosity";
    constexpr static auto referencePorosityString = "referencePorosity";

    constexpr static auto gravityFlagString = "gravityFlag";
    constexpr static auto gravityDepthString = "gravityDepth";

    constexpr static auto permeabilityString = "permeability";

    constexpr static auto discretizationString = "discretization";

    dataRepository::ViewKey blockLocalDofNumber = { blockLocalDofNumberString };
    dataRepository::ViewKey functionalSpace = { "functionalSpace" };
  } viewKeys;

  struct groupKeyStruct : SolverBase::groupKeyStruct
  {
  } groupKeys;

private:

  /**
   * @brief This function generates various discretization information for later use.
   * @param domain the domain parition
   */
  void PrecomputeData(DomainPartition *const domain);

  /**
   * @brief This function allocates additional storage (e.g. for derivatives)
   * @param domain the domain partition
   */
  void AllocateAuxStorage(DomainPartition *const domain);

  /// flag indicating whether FV precompute has been performed
  bool m_precomputeDone;

  /// flag to determine whether or not to apply gravity
  integer m_gravityFlag;

  /// name of the FV discretization object in the data repository
  std::string m_discretizationName;

  /// temp storage for derivatives of density w.r.t. pressure
  array<array<array<real64>>> m_dDens_dPres;

  /// temp storage for derivatives of porosity w.r.t. pressure
  array<array<array<real64>>> m_dPoro_dPres;

  /// temp storage for derivatives of porosity w.r.t. pressure
  array<array<array<real64>>> m_dVisc_dPres;

};


} /* namespace geosx */

#endif /*  */
