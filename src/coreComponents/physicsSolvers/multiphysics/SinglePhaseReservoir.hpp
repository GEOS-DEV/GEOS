/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SinglePhaseReservoir.hpp
 *
 */

#ifndef GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASERESERVOIR_HPP_
#define GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASERESERVOIR_HPP_

#include "physicsSolvers/multiphysics/ReservoirSolverBase.hpp"

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
  SinglePhaseReservoir( const string & name,
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
  static string catalogName() { return "SinglePhaseReservoir"; }

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/

  virtual void
  setupSystem( DomainPartition & domain,
               DofManager & dofManager,
               CRSMatrix< real64, globalIndex > & localMatrix,
               array1d< real64 > & localRhs,
               array1d< real64 > & localSolution,
               bool const setSparsity = true ) override;

  /**@}*/

  virtual void addCouplingSparsityPattern( DomainPartition const & domain,
                                           DofManager const & dofManager,
                                           SparsityPatternView< globalIndex > const & pattern ) const override;

  virtual void assembleCouplingTerms( real64 const time_n,
                                      real64 const dt,
                                      DomainPartition const & domain,
                                      DofManager const & dofManager,
                                      CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                      arrayView1d< real64 > const & localRhs ) override;

protected:

  virtual void initializePostInitialConditionsPreSubGroups() override;

};

} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASERESERVOIR_HPP_ */
