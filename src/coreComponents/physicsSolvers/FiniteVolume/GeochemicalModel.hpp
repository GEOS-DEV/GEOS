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
 * @file GeochemicalModel.hpp
 */

#ifndef SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_GEOCHEMICALMODEL_HPP_
#define SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_GEOCHEMICALMODEL_HPP_

#include "physicsSolvers/FiniteVolume/FlowSolverBase.hpp"
#include "constitutive/Fluid/ReactiveFluidBase.hpp"

class Epetra_FECrsGraph;

namespace geosx
{

namespace dataRepository
{
class ManagedGroup;
}
class FieldSpecificationBase;
class FiniteElementBase;
class DomainPartition;

/**
 * @class GeochemicalModel
 *
 * class to perform a proppant finite volume solve.
 */
class GeochemicalModel : public FlowSolverBase
{
public:
  /**
   * @brief main constructor for ManagedGroup Objects
   * @param name the name of this instantiation of ManagedGroup in the repository
   * @param parent the parent group of this instantiation of ManagedGroup
   */
  GeochemicalModel( const std::string& name,
                   ManagedGroup * const parent );


  /// deleted default constructor
  GeochemicalModel() = delete;

  /// deleted copy constructor
  GeochemicalModel( GeochemicalModel const & ) = delete;

  /// default move constructor
  GeochemicalModel( GeochemicalModel && ) = default;

  /// deleted assignment operator
  GeochemicalModel & operator=( GeochemicalModel const & ) = delete;

  /// deleted move operator
  GeochemicalModel & operator=( GeochemicalModel && ) = delete;

  /**
   * @brief default destructor
   */
  virtual ~GeochemicalModel() override = default;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string CatalogName() { return "GeochemicalModel"; }

  virtual void InitializePreSubGroups(ManagedGroup * const rootGroup) override;

  virtual void RegisterDataOnMesh(ManagedGroup * const MeshBodies) override;

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
                               real64 const time_n,
                               real64 const dt ) override;

  virtual void ApplyBoundaryConditions( DomainPartition * const domain,
                                        systemSolverInterface::EpetraBlockSystem * const blockSystem,
                                        real64 const time_n,
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

  /**
   * @brief assembles the accumulation terms for all cells
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
   * @brief assembles the flux terms for all cells
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


  void ResizeFields( MeshLevel * const meshLevel );

  void WriteSpeciesToFile(DomainPartition * const domain);


  /**@}*/


  struct viewKeyStruct : FlowSolverBase::viewKeyStruct
  {

    static constexpr auto blockLocalDofNumberString = "blockLocalDofNumber_GeochemicalModel" ;

    static constexpr auto reactiveFluidNameString      = "reactiveFluidName";
    static constexpr auto reactiveFluidIndexString      = "reactiveFluidIndex";
    

    // primary solution field
    static constexpr auto pressureString      = "pressure";
    static constexpr auto deltaPressureString = "deltaPressure";

    static constexpr auto temperatureString      = "temperature";
    static constexpr auto deltaTemperatureString = "deltaTemperature";

    static constexpr auto concentrationString      = "concentration";
    static constexpr auto deltaConcentrationString      = "deltaConcentration";

    static constexpr auto totalConcentrationString      = "totalConcentration";
    static constexpr auto concentrationNewString      = "concentrationNew";
    static constexpr auto outputSpeciesFileNameString      = "outputSpeciesFileName";        

    using ViewKey = dataRepository::ViewKey;

    // dof numbering
    ViewKey blockLocalDofNumber = { blockLocalDofNumberString };

    // primary solution field
    ViewKey pressure      = { pressureString };
    ViewKey deltaPressure = { deltaPressureString };

    ViewKey temperature      = { temperatureString };
    ViewKey deltaTemperature = { deltaTemperatureString };    

    ViewKey concentration      = { concentrationString };
    ViewKey deltaConcentration = { deltaConcentrationString };

    ViewKey concentrationNew      = { concentrationNewString };

    ViewKey outputSpeciesFileName  = { outputSpeciesFileNameString };        
    
  } viewKeysGeochemicalModel;

  viewKeyStruct & viewKeys() { return viewKeysGeochemicalModel; }
  viewKeyStruct const & viewKeys() const { return viewKeysGeochemicalModel; }

  struct groupKeyStruct : SolverBase::groupKeyStruct
  {
  } groupKeysGeochemicalModel;

  groupKeyStruct & groupKeys() { return groupKeysGeochemicalModel; }
  groupKeyStruct const & groupKeys() const { return groupKeysGeochemicalModel; }

protected:

  virtual void InitializePostInitialConditions_PreSubGroups( dataRepository::ManagedGroup * const rootGroup ) override;

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
                                     localIndex offset );

  /**
   * @brief Setup stored views into domain data for the current step
   */
  void ResetViews( DomainPartition * const domain ) override;

  /**
   * @brief Function to update all constitutive models
   * @param domain the domain
   */
  void UpdateReactiveFluidModel( ManagedGroup * const dataGroup );

  void UpdateState( ManagedGroup * dataGroup );

  /// views into primary variable fields

  ElementRegionManager::ElementViewAccessor<arrayView1d<globalIndex>> m_dofNumber; // TODO will move to DofManager

  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> m_pressure;
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> m_deltaPressure;

  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> m_temperature;
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> m_deltaTemperature;  

  ElementRegionManager::ElementViewAccessor<arrayView2d<real64>> m_concentration;
  ElementRegionManager::ElementViewAccessor<arrayView2d<real64>> m_deltaConcentration;  

  ElementRegionManager::ElementViewAccessor<arrayView2d<real64>> m_totalConcentration;

  ElementRegionManager::ElementViewAccessor<arrayView2d<real64>> m_concentrationNew;  
  
  /// views into material fields

  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> m_dependentConc;
  ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> m_dDependentConc_dConc;

  string m_reactiveFluidName;
  localIndex m_reactiveFluidIndex;  
  localIndex m_numBasisSpecies;
  localIndex m_numDependentSpecies;  

  string m_outputSpeciesFileName;
  
  //Below is not used in GeochemicalModel model
  
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> m_porosity;
  
  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> m_pvMult;
  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> m_dPvMult_dPres;
  
};


} /* namespace geosx */

#endif //SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_GEOCHEMICALMODEL_HPP_
