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
 * @file CompositionalMultiphaseReservoir.hpp
 *
 */

#ifndef SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_COUPLEDSOLVERS_COMPOSITIONALMULTIPHASERESERVOIR_HPP_
#define SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_COUPLEDSOLVERS_COMPOSITIONALMULTIPHASERESERVOIR_HPP_

#include "physicsSolvers/CoupledSolvers/ReservoirSolverBase.hpp"

namespace geosx
{

class CompositionalMultiphaseReservoir : public ReservoirSolverBase
{
public:

  /**
   * @brief main constructor for ManagedGroup Objects
   * @param name the name of this instantiation of ManagedGroup in the repository
   * @param parent the parent group of this instantiation of ManagedGroup
   */
  CompositionalMultiphaseReservoir( const std::string& name,
                                    Group * const parent );

  /**
   * @brief default destructor
   */
  virtual ~CompositionalMultiphaseReservoir() override;

  /// deleted copy constructor
  CompositionalMultiphaseReservoir( CompositionalMultiphaseReservoir const & ) = delete;

  /// default move constructor
  CompositionalMultiphaseReservoir( CompositionalMultiphaseReservoir && ) = default;

  /// deleted assignment operator
  CompositionalMultiphaseReservoir & operator=( CompositionalMultiphaseReservoir const & ) = delete;

  /// deleted move operator
  CompositionalMultiphaseReservoir & operator=( CompositionalMultiphaseReservoir && ) = delete;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string CatalogName() { return "CompositionalMultiphaseReservoir"; }

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/
  
  virtual void SetupSystem( DomainPartition * const domain,
                            DofManager & dofManager,
                            ParallelMatrix & matrix,
                            ParallelVector & rhs,
                            ParallelVector & solution ) override;

  virtual void SetupDofs( DomainPartition const * const domain,
                          DofManager & dofManager ) const override;
  
protected:

  virtual void AssembleCouplingTerms( real64 const time_n,
                                      real64 const dt,
                                      DomainPartition const * const domain,
                                      DofManager const * const dofManager,
                                      ParallelMatrix * const matrix,
                                      ParallelVector * const rhs ) override;

  virtual void ResetViews( DomainPartition * const domain ) override;

  /**
   * @brief Compute all the perforation rates for this well
   * @param well the well with its perforations
   */
  void ComputeAllPerforationRates( WellElementSubRegion const * const subRegion );

  /// views into reservoir primary variable fields

  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> m_resPressure;
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> m_deltaResPressure;

  /// views into other reservoir variable fields

  ElementRegionManager::ElementViewAccessor<arrayView2d<real64>> m_dResPhaseVolFrac_dPres;
  ElementRegionManager::ElementViewAccessor<arrayView3d<real64>> m_dResPhaseVolFrac_dCompDens;

  ElementRegionManager::ElementViewAccessor<arrayView3d<real64>> m_dResCompFrac_dCompDens;

  ElementRegionManager::ElementViewAccessor<arrayView2d<real64>> m_resPhaseMob;
  ElementRegionManager::ElementViewAccessor<arrayView2d<real64>> m_dResPhaseMob_dPres;
  ElementRegionManager::ElementViewAccessor<arrayView3d<real64>> m_dResPhaseMob_dCompDens;

  /// views into reservoir material fields

  ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> m_resPhaseVisc;
  ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> m_dResPhaseVisc_dPres;
  ElementRegionManager::MaterialViewAccessor<arrayView4d<real64>> m_dResPhaseVisc_dComp;

  ElementRegionManager::MaterialViewAccessor<arrayView4d<real64>> m_resPhaseCompFrac;
  ElementRegionManager::MaterialViewAccessor<arrayView4d<real64>> m_dResPhaseCompFrac_dPres;
  ElementRegionManager::MaterialViewAccessor<arrayView5d<real64>> m_dResPhaseCompFrac_dComp;

  ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> m_resPhaseRelPerm;
  ElementRegionManager::MaterialViewAccessor<arrayView4d<real64>> m_dResPhaseRelPerm_dPhaseVolFrac;


};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_COUPLEDSOLVERS_COMPOSITIONALMULTIPHASERESERVOIR_HPP_ */
