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
 * @file SinglePhaseReservoir.hpp
 *
 */

#ifndef SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_COUPLEDSOLVERS_SINGLEPHASERESERVOIR_HPP_
#define SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_COUPLEDSOLVERS_SINGLEPHASERESERVOIR_HPP_

#include "physicsSolvers/CoupledSolvers/ReservoirSolverBase.hpp"

namespace geosx
{

class SinglePhaseReservoir : public ReservoirSolverBase
{
public:

  /**
   * @brief main constructor for ManagedGroup Objects
   * @param name the name of this instantiation of ManagedGroup in the repository
   * @param parent the parent group of this instantiation of ManagedGroup
   */
  SinglePhaseReservoir( const std::string& name,
                        Group * const parent );

  /**
   * @brief default destructor
   */
  virtual ~SinglePhaseReservoir() override;

  /// deleted copy constructor
  SinglePhaseReservoir( SinglePhaseReservoir const & ) = delete;

  /// default move constructor
  SinglePhaseReservoir( SinglePhaseReservoir && ) = default;

  /// deleted assignment operator
  SinglePhaseReservoir & operator=( SinglePhaseReservoir const & ) = delete;

  /// deleted move operator
  SinglePhaseReservoir & operator=( SinglePhaseReservoir && ) = delete;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string CatalogName() { return "SinglePhaseReservoir"; }

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

  /**@}*/

  virtual void AssembleCouplingTerms( real64 const time_n,
                                      real64 const dt,
                                      DomainPartition const * const domain,
                                      DofManager const * const dofManager,
                                      ParallelMatrix * const matrix,
                                      ParallelVector * const rhs ) override;

protected:

  virtual void ResetViews( DomainPartition * const domain ) override;

  /**
   * @brief Compute all the perforation rates for this well
   * @param well the well with its perforations
   */
  void ComputeAllPerforationRates( WellElementSubRegion const * const subRegion );

private:

  /// views into reservoir primary variable fields

  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> m_resPressure;
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> m_deltaResPressure;

  /// views into reservoir material fields

  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> m_resDensity;
  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> m_dResDens_dPres;

  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> m_resViscosity;
  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> m_dResVisc_dPres;


};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_COUPLEDSOLVERS_SINGLEPHASERESERVOIR_HPP_ */
