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

#ifndef SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_COMPOSITIONALMULTIPHASEFLOW_HPP_
#define SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_COMPOSITIONALMULTIPHASEFLOW_HPP_

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
 * @class CompositionalMultiphaseFlow
 *
 * A compositional multiphase solver
 */
class CompositionalMultiphaseFlow : public SolverBase
{
public:

  /**
   * @brief main constructor for ManagedGroup Objects
   * @param name the name of this instantiation of ManagedGroup in the repository
   * @param parent the parent group of this instantiation of ManagedGroup
   */
  CompositionalMultiphaseFlow( const std::string& name,
                               ManagedGroup * const parent );

  /// deleted default constructor
  CompositionalMultiphaseFlow() = delete;

  /// deleted copy constructor
  CompositionalMultiphaseFlow( CompositionalMultiphaseFlow const & ) = delete;

  /// default move constructor
  CompositionalMultiphaseFlow( CompositionalMultiphaseFlow && ) = default;

  /// deleted assignment operator
  CompositionalMultiphaseFlow & operator=( CompositionalMultiphaseFlow const & ) = delete;

  /// deleted move operator
  CompositionalMultiphaseFlow & operator=( CompositionalMultiphaseFlow && ) = delete;

  /**
   * @brief default destructor
   */
  virtual ~CompositionalMultiphaseFlow() override = default;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string CatalogName() { return "CompositionalMultiphaseFlow"; }

  virtual void FillDocumentationNode() override final;

  virtual void FillOtherDocumentationNodes( dataRepository::ManagedGroup * const rootGroup ) override final;

  virtual void FinalInitialization( dataRepository::ManagedGroup * const rootGroup ) override final;

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

  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    constexpr static auto blockLocalDofNumberString = "blockLocalDofNumber_CompositionalMultiphaseFlow";

    /* primary variables */

    constexpr static auto fluidPressureString = "fluidPressure";
    constexpr static auto deltaFluidPressureString = "deltaFluidPressure";

    constexpr static auto componentMassDensityString = "componentMassDensity";
    constexpr static auto deltaComponentMasslDensityString = "deltaComponentMassDensity";

    constexpr static auto totalMassDensityString = "totalMassDensity";
    constexpr static auto deltaTotalMassDensityString = "deltaTotalMassDensity";

    /* dependent variables */

    constexpr static auto componentMoleFractionString = "componentMoleFraction";

    constexpr static auto phaseComponentMassFractionString = "phaseComponentMassFraction";
    constexpr static auto deltaPhaseComponentMassFractionString = "deltaPhaseComponentMassFraction";

    constexpr static auto phaseVolumeFractionString = "phaseVolumeFraction";
    constexpr static auto deltaPhaseVolumeFractionString = "deltaPhaseVolumeFraction";

    constexpr static auto phaseMassDensityString = "phaseMassDensityString";
    constexpr static auto deltaPhaseMassDensityString = "deltaPhaseMassDensityString";

    constexpr static auto phaseViscosityString = "phaseViscosity";
    constexpr static auto deltaPhaseViscosityString = "deltaPhaseViscosity";

    constexpr static auto phaseRelativePermeabilityString = "phaseRelativePermeability";
    constexpr static auto deltaPhaseRelativePermeabilityString = "deltaPhaseRelativePermeability";

    constexpr static auto porosityString = "porosity";
    constexpr static auto deltaPorosityString = "deltaPorosity";

    /* other data */

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

};

} // namespace geosx


#endif //SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_COMPOSITIONALMULTIPHASEFLOW_HPP_
