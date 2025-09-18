/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */


/**
 * @file GravityFE.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_GRAVITY_GRAVITYFE_HPP_
#define GEOS_PHYSICSSOLVERS_GRAVITY_GRAVITYFE_HPP_

#include "GravitySolverBase.hpp"
#include "mesh/MeshFields.hpp"

namespace geos
{


class GravityFE : public GravitySolverBase
{
public:

  using EXEC_POLICY = parallelDevicePolicy< 32 >;
  using ATOMIC_POLICY = parallelDeviceAtomic;


  GravityFE() = delete;

  GravityFE( const std::string & name,
             Group * const parent );

  virtual ~GravityFE() override;


  GravityFE( GravityFE const & ) = delete;
  GravityFE( GravityFE && ) = delete;

  GravityFE & operator=( GravityFE const & ) = delete;
  GravityFE & operator=( GravityFE && ) = delete;


  static string catalogName() { return "GravityFE"; }
  string getCatalogName() const override { return catalogName(); }


  virtual void initializePreSubGroups() override;

  virtual void registerDataOnMesh( Group & meshBodies ) override final;


  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/
  virtual
  real64 explicitStepModeling( real64 const & time_n,
                               real64 const & dt,
                               integer const cycleNumber,
                               DomainPartition & domain ) override;
  virtual
  real64 explicitStepAdjoint( real64 const & time_n,
                              real64 const & dt,
                              integer const cycleNumber,
                              DomainPartition & domain ) override;
  /**@}*/

};

} // namespace geos

#endif // GEOS_PHYSICSSOLVERS_GRAVITY_GRAVITYFE_HPP_
