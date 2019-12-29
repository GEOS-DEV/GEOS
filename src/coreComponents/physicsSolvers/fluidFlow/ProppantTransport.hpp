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
 * @file ProppantTransport.hpp
 */

#ifndef SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_PROPPANTTRANSPORT_HPP_
#define SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_PROPPANTTRANSPORT_HPP_

#include "physicsSolvers/fluidFlow/FlowSolverBase.hpp"
#include "constitutive/fluid/SlurryFluidBase.hpp"

namespace geosx
{

namespace dataRepository
{
class Group;
}
class FieldSpecificationBase;
class FiniteElementBase;
class DomainPartition;

/**
 * @class ProppantTransport
 *
 * class to perform a proppant finite volume solve.
 */
class ProppantTransport : public FlowSolverBase
{
public:
  /**
   * @brief main constructor for Group Objects
   * @param name the name of this instantiation of Group in the repository
   * @param parent the parent group of this instantiation of Group
   */
  ProppantTransport( const std::string& name,
                   Group * const parent );


  /// deleted default constructor
  ProppantTransport() = delete;

  /// deleted copy constructor
  ProppantTransport( ProppantTransport const & ) = delete;

  /// default move constructor
  ProppantTransport( ProppantTransport && ) = default;

  /// deleted assignment operator
  ProppantTransport & operator=( ProppantTransport const & ) = delete;

  /// deleted move operator
  ProppantTransport & operator=( ProppantTransport && ) = delete;

  /**
   * @brief default destructor
   */
  virtual ~ProppantTransport() override = default;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string CatalogName() { return "ProppantTransport"; }

  virtual void InitializePreSubGroups(Group * const rootGroup) override;
  
  virtual void RegisterDataOnMesh(Group * const MeshBodies) override;

  virtual real64 SolverStep( real64 const& time_n,
                             real64 const& dt,
                             integer const cycleNumber,
                             DomainPartition * domain ) override;

  void PreStepUpdate( real64 const& time_n,
                      real64 const& dt,
                      integer const cycleNumber,
                      DomainPartition * domain );  

  void PostStepUpdate( real64 const& time_n,
                       real64 const& dt,
                       integer const cycleNumber,
                       DomainPartition * domain );  
  
  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/

  virtual void ImplicitStepSetup( real64 const & time_n,
                                  real64 const & dt,
                                  DomainPartition * const domain,
                                  DofManager & dofManager,
                                  ParallelMatrix & matrix,
                                  ParallelVector & rhs,
                                  ParallelVector & solution ) override;

  virtual void
  SetupDofs( DomainPartition const * const domain,
             DofManager & dofManager ) const override;

  virtual void AssembleSystem( real64 const time,
                               real64 const dt,
                               DomainPartition * const domain,
                               DofManager const & dofManager,
                               ParallelMatrix & matrix,
                               ParallelVector & rhs ) override;

  virtual void ApplyBoundaryConditions( real64 const time,
                                        real64 const dt,
                                        DomainPartition * const domain,
                                        DofManager const & dofManager,
                                        ParallelMatrix & matrix,
                                        ParallelVector & rhs ) override;

  virtual real64
  CalculateResidualNorm( DomainPartition const * const domain,
                         DofManager const & dofManager,
                         ParallelVector const & rhs ) override;

  virtual void SolveSystem( DofManager const & dofManager,
                            ParallelMatrix & matrix,
                            ParallelVector & rhs,
                            ParallelVector & solution ) override;

  virtual void
  ApplySystemSolution( DofManager const & dofManager,
                       ParallelVector const & solution,
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

  void AssembleAccumulationTerms( DomainPartition const * const domain,
                                  DofManager const * const dofManager,
                                  ParallelMatrix * const matrix,
                                  ParallelVector * const rhs );

  /**
   * @brief assembles the flux terms for all cells
   * @param domain the physical domain object
   * @param blockSystem the entire block system
   * @param time_n previous time value
   * @param dt time step
   */
  void AssembleFluxTerms( real64 const time_n,
                          real64 const dt,
                          DomainPartition const * const domain,
                          DofManager const * const dofManager,
                          ParallelMatrix * const matrix,
                          ParallelVector * const rhs );

  /**@}*/

  void ResizeFractureFields(  real64 const & time_n,
                              real64 const & dt,
                              DomainPartition * const domain);  
  
  
  struct viewKeyStruct : FlowSolverBase::viewKeyStruct
  {

    static constexpr auto proppantNameString      = "proppantName";
    static constexpr auto proppantIndexString      = "proppantIndex";
    
    // primary solution field
    static constexpr auto proppantConcentrationString      = "proppantConcentration";
    static constexpr auto deltaProppantConcentrationString      = "deltaProppantConcentration";    

    static constexpr auto componentConcentrationString      = "componentConcentration";
    
    static constexpr auto deltaComponentConcentrationString      = "deltaComponentConcentration";    

    static constexpr auto bcComponentConcentrationString      = "bcComponentConcentration";
    
    static constexpr auto updatedComponentConcentrationString      = "updatedComponentConcentration";
    
    // these are used to store last converged time step values
    static constexpr auto oldProppantConcentrationString  = "oldProppantConcentration";
    static constexpr auto oldComponentDensityString  = "oldComponentDensity";    
    
    static constexpr auto updateProppantPackingString  = "updateProppantPacking";

    static constexpr auto cellBasedFluxString  = "cellBasedFlux";    

    static constexpr auto isInterfaceElementString   = "isInterfaceElement";            

    static constexpr auto isProppantBoundaryString   = "isProppantBoundary";

    static constexpr auto isProppantMobileString   = "isProppantMobile";                    
    
    static constexpr auto proppantPackVolumeFractionString  = "proppantPackVolumeFraction";

    static constexpr auto proppantExcessPackVolumeString  = "proppantExcessPackVolume";

    static constexpr auto proppantLiftFluxString  = "proppantLiftFlux";

    static constexpr auto poroMultiplierString  = "poroMultiplier";

    static constexpr auto transTMultiplierString  = "transTMultiplier";            

    static constexpr auto bridgingFactorString  = "bridgingFactor";

    static constexpr auto maxProppantConcentrationString  = "maxProppantConcentration";

    static constexpr auto proppantDiameterString  = "proppantDiameter";

    static constexpr auto proppantDensityString  = "proppantDensity";

    static constexpr auto criticalShieldsNumberString  = "criticalShieldsNumber";

    static constexpr auto frictionCoefficientString  = "frictionCoefficient";                                    

    /*
    using ViewKey = dataRepository::ViewKey;

    // primary solution field
    ViewKey proppantConcentration      = { proppantConcentrationString };
    ViewKey deltaProppantConcentration = { deltaProppantConcentrationString };

    ViewKey componentConcentration      = { componentConcentrationString };
    ViewKey deltaComponentConcentration = { deltaComponentConcentrationString };    
    
    ViewKey updatedComponentConcentration      = { updatedComponentConcentrationString };

    
    ViewKey oldProppantConcentration    = { oldProppantConcentrationString };
    ViewKey oldComponentDensity    = { oldComponentDensityString };
    

    ViewKey proppantName      = { proppantNameString };
    ViewKey proppantIndex      = { proppantIndexString };
    */
    
  } viewKeysProppantTransport;

  viewKeyStruct & viewKeys() { return viewKeysProppantTransport; }
  viewKeyStruct const & viewKeys() const { return viewKeysProppantTransport; }

  struct groupKeyStruct : SolverBase::groupKeyStruct
  {
  } groupKeysProppantTransport;

  groupKeyStruct & groupKeys() { return groupKeysProppantTransport; }
  groupKeyStruct const & groupKeys() const { return groupKeysProppantTransport; }

  static constexpr localIndex MAX_NUM_COMPONENTS = 3;

protected:

  virtual void InitializePostInitialConditions_PreSubGroups( dataRepository::Group * const rootGroup ) override;

private:

  /**
   * @brief Setup stored views into domain data for the current step
   */
  void ResetViews( DomainPartition * const domain ) override;

  /**
   * @brief Function to update all constitutive models
   * @param domain the domain
   */
  void UpdateFluidModel( Group * const dataGroup );

  void UpdateComponentDensity( Group * const dataGroup );  

  void UpdateProppantModel( Group * const dataGroup );

  void UpdateProppantMobility( Group * const dataGroup );

  void UpdateProppantPackVolume( real64 const time_n,
                                 real64 const dt,
                                 DomainPartition * const domain);

  void UpdateCellBasedFlux( real64 const time_n,
                            DomainPartition * const domain);  
  
  void UpdateState( Group * dataGroup );

  /// views into primary variable fields

  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> m_pressure;
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> m_deltaPressure;

  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> m_proppantConcentration;
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> m_deltaProppantConcentration;

  ElementRegionManager::ElementViewAccessor<arrayView2d<real64>> m_componentConcentration;
  ElementRegionManager::ElementViewAccessor<arrayView2d<real64>> m_deltaComponentConcentration;    

  ElementRegionManager::ElementViewAccessor<arrayView2d<real64>> m_updatedComponentConcentration;

  ElementRegionManager::ElementViewAccessor<arrayView1d<R1Tensor>> m_cellBasedFlux;

  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> m_proppantLiftFlux;
  
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> m_proppantPackVolumeFraction;

  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> m_proppantExcessPackVolume;

  ElementRegionManager::ElementViewAccessor<arrayView1d<integer>> m_isProppantBoundaryElement;

  ElementRegionManager::ElementViewAccessor<arrayView1d<integer>> m_isInterfaceElement;      
  
  ElementRegionManager::ElementViewAccessor<arrayView1d<integer>> m_isProppantMobile;

  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> m_poroMultiplier;

  ElementRegionManager::ElementViewAccessor<arrayView1d<R1Tensor>> m_transTMultiplier;    

  /// views into backup fields

  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> m_proppantConcentrationOld;
  ElementRegionManager::ElementViewAccessor<arrayView2d<real64>> m_componentDensityOld;  

  /// views into material fields

  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> m_density;
  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> m_dDensity_dPressure;
  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> m_dDensity_dProppantConcentration;
  ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> m_dDensity_dComponentConcentration;

  ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> m_componentDensity;
  ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> m_dComponentDensity_dPressure;
  ElementRegionManager::MaterialViewAccessor<arrayView4d<real64>> m_dComponentDensity_dComponentConcentration;
  
  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> m_fluidDensity;  
  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> m_dFluidDensity_dPressure;
  ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> m_dFluidDensity_dComponentConcentration;      

  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> m_fluidViscosity;  
  
  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> m_viscosity;
  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> m_dViscosity_dPressure;
  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> m_dViscosity_dProppantConcentration;
  ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> m_dViscosity_dComponentConcentration;  

  ElementRegionManager::MaterialViewAccessor<arrayView1d<real64>> m_settlingFactor;
  ElementRegionManager::MaterialViewAccessor<arrayView1d<real64>> m_dSettlingFactor_dPressure;
  ElementRegionManager::MaterialViewAccessor<arrayView1d<real64>> m_dSettlingFactor_dProppantConcentration;
  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> m_dSettlingFactor_dComponentConcentration;  

  ElementRegionManager::MaterialViewAccessor<arrayView1d<real64>> m_collisionFactor;
  ElementRegionManager::MaterialViewAccessor<arrayView1d<real64>> m_dCollisionFactor_dProppantConcentration;

  string m_proppantName;
  localIndex m_proppantIndex;  
  integer m_numComponents;

  R1Tensor m_downVector;

  integer m_updateProppantPacking;
  real64 m_proppantPackPermeability;
  real64 m_bridgingFactor;
  real64 m_minAperture;  
  real64 m_maxProppantConcentration;
  real64 m_proppantDiameter;
  real64 m_proppantDensity;
  real64 m_criticalShieldsNumber;        
  real64 m_frictionCoefficient;          
};


} /* namespace geosx */

#endif //SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_PROPPANTTRANSPORT_HPP_
