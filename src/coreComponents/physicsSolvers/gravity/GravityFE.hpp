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

#include "physicsSolvers/PhysicsSolverBase.hpp"
#include "GravitySolverBase.hpp"



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
  GravityFE( GravityFE && ) = default;

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


protected:

  virtual void postInputInitialization() override final;

  virtual void initializePostInitialConditionsPreSubGroups() override final;


private:


};



namespace fields
{
DECLARE_FIELD( MediumDensity,
               "mediumDensity",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Medium density of the cell" );

DECLARE_FIELD( VolumeIntegral,
               "volumeIntegral",
               array1d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "VolumeIntegral of the cell." );

DECLARE_FIELD( VolumeIntegral2d,
               "volumeIntegral2d",
               array2d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "VolumeIntegral for adjoint computation." );

DECLARE_FIELD( Adjoint,
               "adjoint",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Adjoint field." );
}



} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_GRAVITY_GRAVITYFE_HPP_ */
