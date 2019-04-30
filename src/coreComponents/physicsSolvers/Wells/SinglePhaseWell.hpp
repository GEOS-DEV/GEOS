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
class SinglePhaseFlow;
class Well;
  
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

  // define the column offset of the derivatives
  struct ColOffset
  {
    static constexpr integer DPRES = 0;
    static constexpr integer DRATE = 1;
  };

  // define the row offset of the residual equations
  struct RowOffset
  {
    static constexpr integer CONTROL = 0;
    static constexpr integer MASSBAL = 1;
  };
  
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

  virtual void RegisterDataOnMesh(ManagedGroup * const meshBodies) override;

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/

  virtual real64
  CalculateResidualNorm( systemSolverInterface::EpetraBlockSystem const *const blockSystem,
                         DomainPartition *const domain) override;
  
  virtual bool
  CheckSystemSolution( systemSolverInterface::EpetraBlockSystem const * const blockSystem, 
                       real64 const scalingFactor,
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
   * @brief Update all relevant fluid models using current values of pressure
   * @param well the well containing the well elements and their associated fields
   */
  void UpdateFluidModel( Well * well );

  /**
   * @brief Recompute all dependent quantities from primary variables (including constitutive models) on the well
   * @param well the well containing the well elements and their associated fields
   */
  virtual void UpdateState( Well * well ) override;
  
  /**
   * @brief assembles the flux terms for all connections between well elements
   * @param domain the physical domain object
   * @param jacobian the entire jacobian matrix of the system
   * @param residual the entire residual of the system
   * @param time_n previous time value
   * @param dt time step
   */
  virtual void AssembleFluxTerms( DomainPartition * const domain,
                                  Epetra_FECrsMatrix * const jacobian,
                                  Epetra_FEVector * const residual,
                                  real64 const time_n,
                                  real64 const dt ) override;


  /**
   * @brief assembles the perforation rate terms 
   * @param domain the physical domain object
   * @param jacobian the entire jacobian matrix of the system
   * @param residual the entire residual of the system
   * @param time_n previous time value
   * @param dt time step
   */
  virtual void AssemblePerforationTerms( DomainPartition * const domain,
                                         Epetra_FECrsMatrix * const jacobian,
                                         Epetra_FEVector * const residual,
                                         real64 const time_n,
                                         real64 const dt ) override;
  
  /**
   * @brief assembles the volume balance terms for all well elements
   * @param domain the physical domain object
   * @param jacobian the entire jacobian matrix of the system
   * @param residual the entire residual of the system
   * @param time_n previous time value
   * @param dt time step
   */
  virtual void AssembleVolumeBalanceTerms( DomainPartition * const domain,
                                           Epetra_FECrsMatrix * const jacobian,
                                           Epetra_FEVector * const residual,
                                           real64 const time_n,
                                           real64 const dt ) override;

  /**
   * @brief assembles the pressure relations at all connections between well elements except at the well head (first connection)
   * @param domain the physical domain object
   * @param jacobian the entire jacobian matrix of the system
   * @param residual the entire residual of the system
   */
  virtual void FormPressureRelations( DomainPartition * const domain,
                                      Epetra_FECrsMatrix * const jacobian,
                                      Epetra_FEVector * const residual ) override;

  /**
   * @brief assembles the control equation for the well head (first connection)
   * @param domain the physical domain object
   * @param jacobian the entire jacobian matrix of the system
   * @param residual the entire residual of the system
   */
  virtual void FormControlEquation( DomainPartition * const domain,
                                    Epetra_FECrsMatrix * const jacobian,
                                    Epetra_FEVector * const residual ) override;
  
  /**
   * @brief set the sparsity pattern for the linear system
   * @param domain the domain partition
   * @param sparsity the sparsity pattern matrix
   */
   virtual void SetSparsityPattern( DomainPartition const * const domain,
                                    Epetra_FECrsGraph * const sparsity,
                                    globalIndex firstWellElemDofNumber,
                                    localIndex numDofPerResElement) override;

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
  virtual void SetNumRowsAndTrilinosIndices( DomainPartition const * const domain,
                                             localIndex & numLocalRows,
                                             globalIndex & numGlobalRows,
                                             localIndex offset ) override;
  
  /**@}*/

  struct viewKeyStruct : WellSolverBase::viewKeyStruct
  {
    // degrees of freedom numbers on the well elements
    static constexpr auto dofNumberString = "segmentLocalDofNumber";

    // primary solution field
    static constexpr auto pressureString      = "segmentPressure";
    static constexpr auto deltaPressureString = "deltaSegmentPressure";
    static constexpr auto connRateString      = "connectionRate";
    static constexpr auto deltaConnRateString = "deltaConnectionRate";

    // perforation rates
    static constexpr auto perforationRateString        = "perforationRate";
    static constexpr auto dPerforationRate_dPresString = "dPerforationRate_dPres";
    
    using ViewKey = dataRepository::ViewKey;

    // degrees of freedom numbers of the well elements
    ViewKey dofNumber = { dofNumberString };

    // primary solution field
    ViewKey pressure      = { pressureString };
    ViewKey deltaPressure = { deltaPressureString };
    ViewKey rate          = { connRateString };
    ViewKey deltaRate     = { deltaConnRateString };
    
    // perforation rates
    ViewKey perforationRate        = { perforationRateString };
    ViewKey dPerforationRate_dPres = { dPerforationRate_dPresString };
    
  } viewKeysSinglePhaseWell;

  struct groupKeyStruct : SolverBase::groupKeyStruct
  {
  } groupKeysSinglePhaseWell;

protected:

  virtual void InitializePreSubGroups( ManagedGroup * const rootGroup ) override;

private:

  /**
   * @brief Setup stored reservoir views into domain data for the current step
   * @param domain the domain containing the well manager to access individual wells
   */
  void ResetViews( DomainPartition * const domain ) override;

  /**
   * @brief Initialize all the primary and secondary variables in all the wells
   * @param domain the domain containing the well manager to access individual wells
   */
  void InitializeWells( DomainPartition * const domain ) override;

  /**
   * @brief Check if the controls are viable; if not, switch the controls
   * @param domain the domain containing the well manager to access individual wells
   */
  void CheckWellControlSwitch( DomainPartition * const domain ) override;

  /**
   * @brief Compute all the perforation rates for this well
   * @param well the well with its perforations
   */
  void ComputeAllPerforationRates( Well const * const well );

  /**
   * @brief Save all the rates and pressures in the well for reporting purposes
   * @param well the well with its perforations
   */
  void RecordWellData( Well const * const well );
  
  ElementRegionManager::ElementViewAccessor<arrayView1d<globalIndex>> m_resDofNumber; // TODO will move to DofManager
  
  /// views into reservoir primary variable fields

  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> m_resPressure;
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> m_deltaResPressure;

  /// views into reservoir material fields

  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> m_resDensity;
  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> m_dResDens_dPres;

  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> m_resViscosity;
  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> m_dResVisc_dPres;

};

} // namespace geosx


#endif //SRC_COMPONENTS_CORE_SRC_WELLS_SINGLEPHASEWELL_HPP_
