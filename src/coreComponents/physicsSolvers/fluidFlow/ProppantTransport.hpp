/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
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
  ProppantTransport( const std::string & name,
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

  virtual void InitializePreSubGroups( Group * const rootGroup ) override;

  virtual void RegisterDataOnMesh( Group * const MeshBodies ) override;

  virtual real64 SolverStep( real64 const & time_n,
                             real64 const & dt,
                             integer const cycleNumber,
                             DomainPartition * domain ) override;

  void PreStepUpdate( real64 const & time_n,
                      real64 const & dt,
                      integer const cycleNumber,
                      DomainPartition * domain );

  void PostStepUpdate( real64 const & time_n,
                       real64 const & dt,
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

  virtual void ImplicitStepComplete( real64 const & time,
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

  void ResizeFractureFields( real64 const & time_n,
                             real64 const & dt,
                             DomainPartition * const domain );

  arrayView1d< string const > const & proppantModelNames() const { return m_proppantModelNames; }

  struct viewKeyStruct : FlowSolverBase::viewKeyStruct
  {

    static constexpr auto proppantNamesString      = "proppantNames";

    // primary solution field
    static constexpr auto proppantConcentrationString      = "proppantConcentration";
    static constexpr auto deltaProppantConcentrationString      = "deltaProppantConcentration";
    static constexpr auto componentConcentrationString      = "componentConcentration";
    static constexpr auto deltaComponentConcentrationString      = "deltaComponentConcentration";
    static constexpr auto bcComponentConcentrationString      = "bcComponentConcentration";

    // these are used to store last converged time step values

    static constexpr auto oldComponentDensityString  = "oldComponentDensity";

    static constexpr auto updateProppantPackingString  = "updateProppantPacking";
    static constexpr auto cellBasedFluxString  = "cellBasedFlux";

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

  } viewKeysProppantTransport;

  viewKeyStruct & viewKeys() { return viewKeysProppantTransport; }
  viewKeyStruct const & viewKeys() const { return viewKeysProppantTransport; }

  struct groupKeyStruct : SolverBase::groupKeyStruct
  {} groupKeysProppantTransport;

  groupKeyStruct & groupKeys() { return groupKeysProppantTransport; }
  groupKeyStruct const & groupKeys() const { return groupKeysProppantTransport; }

  static constexpr localIndex MAX_NUM_COMPONENTS = 3;

protected:

  virtual void PostProcessInput() override;

  virtual void InitializePostInitialConditions_PreSubGroups( dataRepository::Group * const rootGroup ) override;

private:

  /**
   * @brief Setup stored views into domain data for the current step
   */
  void ResetViews( DomainPartition * const domain ) override;

  /**
   * @brief Function to update fluid properties
   * @param domain the domain
   */
  void UpdateFluidModel( Group & dataGroup, localIndex const targetIndex );

  void UpdateComponentDensity( Group & dataGroup, localIndex const targetIndex );

  void UpdateProppantModel( Group & dataGroup, localIndex const targetIndex );

  void UpdateProppantMobility( Group & dataGroup );

  /**
   * @brief Function to update proppant pack volume fraction
   */
  void UpdateProppantPackVolume( real64 const time_n,
                                 real64 const dt,
                                 DomainPartition * const domain );

  /**
   * @brief Function to update cell-based fluid flux
   */
  void UpdateCellBasedFlux( real64 const time_n,
                            DomainPartition * const domain );

  /**
   * @brief Function to update fluid and proppant properties
   * @param domain the domain
   */
  void UpdateState( Group & dataGroup, localIndex const targetIndex );

  /// views into primary variable fields

  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > m_pressure;
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > m_deltaPressure;

  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > m_proppantConcentration;
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > m_deltaProppantConcentration;

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > m_componentConcentration;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > m_deltaComponentConcentration;

  ElementRegionManager::ElementViewAccessor< arrayView1d< R1Tensor > > m_cellBasedFlux;

  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > m_proppantLiftFlux;

  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > m_proppantPackVolumeFraction;

  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > m_proppantExcessPackVolume;

  ElementRegionManager::ElementViewAccessor< arrayView1d< integer > > m_isProppantBoundaryElement;

  ElementRegionManager::ElementViewAccessor< arrayView1d< integer > > m_isProppantMobile;

  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > m_poroMultiplier;

  ElementRegionManager::ElementViewAccessor< arrayView1d< R1Tensor > > m_transTMultiplier;

  /// views into backup fields

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > m_componentDensityOld;

  /// views into material fields

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > m_density;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > m_dDensity_dPressure;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > m_dDensity_dProppantConcentration;
  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 > > m_dDensity_dComponentConcentration;

  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 > > m_componentDensity;
  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 > > m_dComponentDensity_dPressure;
  ElementRegionManager::ElementViewAccessor< arrayView4d< real64 > > m_dComponentDensity_dComponentConcentration;

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > m_fluidDensity;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > m_dFluidDensity_dPressure;
  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 > > m_dFluidDensity_dComponentConcentration;

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > m_fluidViscosity;

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > m_viscosity;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > m_dViscosity_dPressure;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > m_dViscosity_dProppantConcentration;
  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 > > m_dViscosity_dComponentConcentration;

  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > m_settlingFactor;
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > m_dSettlingFactor_dPressure;
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > m_dSettlingFactor_dProppantConcentration;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > m_dSettlingFactor_dComponentConcentration;

  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > m_collisionFactor;
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > m_dCollisionFactor_dProppantConcentration;

  array1d< string > m_proppantModelNames;

  integer m_numComponents;

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
