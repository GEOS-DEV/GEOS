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
 * @file GravityFE_CompositionalMultiphaseFVM.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_GRAVITY_GRAVITYFE_COMPOSITIONALMULTUPHASEFVM_HPP_
#define GEOS_PHYSICSSOLVERS_GRAVITY_GRAVITYFE_COMPOSITIONALMULTUPHASEFVM_HPP_

#include "GravitySolverBase.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseFVM.hpp"

namespace geos
{

class GravityFE_CompositionalMultiphaseFVM : public GravitySolverBase
{

public:
  using EXEC_POLICY = parallelDevicePolicy< 32 >;
  using ATOMIC_POLICY = parallelDeviceAtomic;

  GravityFE_CompositionalMultiphaseFVM() = delete;
  GravityFE_CompositionalMultiphaseFVM( const std::string & name,
                                        Group * const parent );
  virtual ~GravityFE_CompositionalMultiphaseFVM() override;

  GravityFE_CompositionalMultiphaseFVM( GravityFE_CompositionalMultiphaseFVM const & ) = delete;
  GravityFE_CompositionalMultiphaseFVM( GravityFE_CompositionalMultiphaseFVM && ) = delete;
  GravityFE_CompositionalMultiphaseFVM & operator=( GravityFE_CompositionalMultiphaseFVM const & ) = delete;
  GravityFE_CompositionalMultiphaseFVM & operator=( GravityFE_CompositionalMultiphaseFVM && ) = delete;

  static string catalogName() { return "GravityFE_CompositionalMultiphaseFVM"; }
  string getCatalogName() const override { return catalogName(); }

  virtual void initializePreSubGroups() override;
  virtual void registerDataOnMesh( Group & meshBodies ) override final;
  virtual void postInputInitialization() override final;


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

  struct viewKeyStruct : GravitySolverBase::viewKeyStruct
  {
    static constexpr char const * useRockDensityString() { return "useRockDensity"; }
    static constexpr char const * usePorosityString() { return "usePorosity"; }
    static constexpr char const * useReferencePorosityString() { return "useReferencePorosity"; }
  } gravityViewKeys;


protected:



  virtual void initializePostInitialConditionsPreSubGroups() override final;

  /// Use rock density in addition to fluid density
  localIndex m_useRockDensity;
  localIndex m_useReferencePorosity;
  localIndex m_usePorosity;

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
DECLARE_FIELD( FluidDensity,
               "fluidDensity",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Fluid density of the cell" );
DECLARE_FIELD( RockDensity,
               "rockDensity",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Rock density of the cell" );
DECLARE_FIELD( Porosity,
               "porosity",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Porosity of the cell" );
DECLARE_FIELD( VolumeIntegral,
               "volumeIntegral",
               array1d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "VolumeIntegral of the cell." );

}

} // namespace geos

#endif // GEOS_PHYSICSSOLVERS_GRAVITY_GRAVITYFE_COMPOSITIONALMULTUPHASEFVM_HPP_
