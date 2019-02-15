/*

 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
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
 * @file SinglePhaseWell.hpp
 */

#ifndef SRC_COMPONENTS_CORE_SRC_WELLS_SINGLEPHASEWELLSOLVER_HPP_
#define SRC_COMPONENTS_CORE_SRC_WELLS_SINGLEPHASEWELLSOLVER_HPP_

#include "WellSolverBase.hpp"

class Epetra_FECrsGraph;

namespace geosx
{

namespace dataRepository
{
class ManagedGroup;

namespace keys
{
string const singlePhaseWell = "SinglePhaseWell";
}
}

namespace constitutive
{
class SingleFluidBase;
}

/**
 * @class SinglePhaseWell
 *
 * A single-phase well solver
 */
class SinglePhaseWell : public WellSolverBase
{
public:

  /**
   * @brief main constructor for ManagedGroup Objects
   * @param name the name of this instantiation of ManagedGroup in the repository
   * @param parent the parent group of this instantiation of ManagedGroup
   */
  SinglePhaseWell( const string& name,
                   ManagedGroup * const parent );

  /// deleted default constructor
  SinglePhaseWell() = delete;

  /// deleted copy constructor
  SinglePhaseWell( SinglePhaseWell const & ) = delete;

  /// default move constructor
  SinglePhaseWell( SinglePhaseWell && ) = default;

  /// deleted assignment operator
  SinglePhaseWell & operator=( SinglePhaseWell const & ) = delete;

  /// deleted move operator
  SinglePhaseWell & operator=( SinglePhaseWell && ) = delete;

  /**
   * @brief default destructor
   */
  virtual ~SinglePhaseWell() override = default;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string CatalogName() { return dataRepository::keys::singlePhaseWell; }

  virtual void InitializePreSubGroups( ManagedGroup * const rootGroup ) override;

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
                               real64 const time_n,
                               real64 const dt ) override;

  virtual real64
  CalculateResidualNorm(systemSolverInterface::EpetraBlockSystem const *const blockSystem,
                        DomainPartition *const domain) override;
  
  virtual void SolveSystem( systemSolverInterface::EpetraBlockSystem * const blockSystem,
                            SystemSolverParameters const * const params ) override;
  
  virtual bool
  CheckSystemSolution(systemSolverInterface::EpetraBlockSystem const * const blockSystem, real64 const scalingFactor,
                      DomainPartition * const domain) override;

  virtual void
  ApplySystemSolution( systemSolverInterface::EpetraBlockSystem const * const blockSystem,
                       real64 const scalingFactor,
                       DomainPartition * const domain ) override;
  
  virtual void ResetStateToBeginningOfStep( DomainPartition * const domain ) override;

  virtual void ImplicitStepComplete( real64 const & time,
                                     real64 const & dt,
                                     DomainPartition * const domain ) override;

  /**
   * @brief Update all relevant fluid models using current values of pressure and composition
   * @param domain the domain containing the mesh and fields
   */
  void UpdateFluidModelAll( DomainPartition * domain );

  /**
   * @brief Recompute all dependent quantities from primary variables (including constitutive models)
   * @param domain the domain containing the mesh and fields
   */
  void UpdateStateAll( DomainPartition * domain );

  /**
   * @brief assembles the accumulation terms for all segments
   * @param domain the physical domain object
   * @param blockSystem the entire block system
   * @param time_n previous time value
   * @param dt time step
   */
  void AssembleAccumulationTerms( DomainPartition * const domain,
                                  Epetra_FECrsMatrix * const jacobian,
                                  Epetra_FEVector * const residual,
                                  real64 const time_n,
                                  real64 const dt );

  /**
   * @brief assembles the flux terms for all connections
   * @param domain the physical domain object
   * @param blockSystem the entire block system
   * @param time_n previous time value
   * @param dt time step
   */
  void AssembleFluxTerms( DomainPartition * const domain,
                          Epetra_FECrsMatrix * const jacobian,
                          Epetra_FEVector * const residual,
                          real64 const time_n,
                          real64 const dt );


  /**
   * @brief assembles the perforation rate terms 
   * @param domain the physical domain object
   * @param blockSystem the entire block system
   * @param time_n previous time value
   * @param dt time step
   */

  void AssembleSourceTerms( DomainPartition * const domain,
                            systemSolverInterface::EpetraBlockSystem * const blockSystem,
                            real64 const time_n,
                            real64 const dt );

  void InitializeWellState( DomainPartition * const domain );
  
  void CheckWellControlSwitch( DomainPartition * const domain );
  
  /**@}*/

  struct viewKeyStruct : WellSolverBase::viewKeyStruct
  {
    // primary solution field
    static constexpr auto pressureString      = "pressure";
    static constexpr auto deltaPressureString = "deltaPressure";
    static constexpr auto velocityString      = "velocity";
    static constexpr auto deltaVelocityString = "deltaVelocity";

    static constexpr auto densityString          = "density";
    static constexpr auto dDensity_dPresString   = "dDensity_dPres";
    static constexpr auto viscosityString        = "viscosity";
    static constexpr auto dViscosity_dPresString = "dViscosity_dPres";
    
    static constexpr auto flowRateString = "flowRate";

    using ViewKey = dataRepository::ViewKey;

    // primary solution field
    ViewKey pressure      = { pressureString };
    ViewKey deltaPressure = { deltaPressureString };
    ViewKey velocity      = { velocityString };
    ViewKey deltaVelovity = { deltaVelocityString };

    // density and viscosity
    ViewKey density          = { densityString };
    ViewKey dDensity_dPres   = { dDensity_dPresString };
    ViewKey viscosity        = { viscosityString };
    ViewKey dViscosity_dPres = { dViscosity_dPresString };
    
    // well controls
    ViewKey flowRate = { flowRateString };

  } viewKeysSinglePhaseWell;

  struct groupKeyStruct : SolverBase::groupKeyStruct
  {
  } groupKeysSinglePhaseWell;

protected:

  virtual void InitializePostInitialConditions_PreSubGroups( dataRepository::ManagedGroup * const rootGroup ) override;


private:

  /**
   * @brief Extract the fluid model used by this solver from a group
   * @param dataGroup target group (e.g. subregion, face/edge/node manager, etc.)
   * @return
   */
  constitutive::SingleFluidBase * GetFluidModel( ManagedGroup * const dataGroup ) const;

  /**
   * @brief Extract the fluid model used by this solver from a group (const version)
   * @param dataGroup target group (e.g. subregion, face/edge/node manager, etc.)
   * @return
   */
  constitutive::SingleFluidBase const * GetFluidModel( ManagedGroup const * const dataGroup ) const;
  
  /**
   * @brief Backup current values of all constitutive fields that participate in the accumulation term
   * @param domain the domain containing the mesh and fields
   */
  void BackupFields( DomainPartition * const domain );

  /**
   * @brief Set up the linear system (DOF indices and sparsity patterns)
   * @param domain the domain containing the mesh and fields
   * @param blockSystem the linear system object
   */
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
   * @param offset the DOF offset for this solver in the case of a non-block system
   *
   * This function sets the number of global rows, and sets the dof numbers for
   * this solver. dof numbers are referred to trilinosIndices currently.
   */
  void SetNumRowsAndTrilinosIndices( MeshLevel * const meshLevel,
                                     localIndex & numLocalRows,
                                     globalIndex & numGlobalRows,
                                     localIndex offset );

  /**
   * @brief Setup stored views into domain data for the current step
   */
  void ResetViews( DomainPartition * const domain ) override;

};

} // namespace geosx


#endif //SRC_COMPONENTS_CORE_SRC_WELLS_SINGLEPHASEWELLS_HPP_
