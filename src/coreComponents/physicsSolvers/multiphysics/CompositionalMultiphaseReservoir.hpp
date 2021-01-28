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

#ifndef GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_COMPOSITIONALMULTIPHASERESERVOIR_HPP_
#define GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_COMPOSITIONALMULTIPHASERESERVOIR_HPP_

#include "physicsSolvers/multiphysics/ReservoirSolverBase.hpp"

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
  CompositionalMultiphaseReservoir( const string & name,
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
  static string catalogName() { return "CompositionalMultiphaseReservoir"; }

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/

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

};

} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_COMPOSITIONALMULTIPHASERESERVOIR_HPP_ */
