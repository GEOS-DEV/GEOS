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

  GravityFE_CompositionalMultiphaseFVM( GravityFE_CompositionalMultiphaseFVM const & ) = delete;
  GravityFE_CompositionalMultiphaseFVM( GravityFE_CompositionalMultiphaseFVM && ) = delete;
  GravityFE_CompositionalMultiphaseFVM & operator=( GravityFE_CompositionalMultiphaseFVM const & ) = delete;
  GravityFE_CompositionalMultiphaseFVM & operator=( GravityFE_CompositionalMultiphaseFVM && ) = delete;

  static string catalogName() { return "GravityFE_CompositionalMultiphaseFVM"; }
  string getCatalogName() const override { return catalogName(); }

  virtual void registerDataOnMesh( Group & meshBodies ) override final;
  virtual void initializePreSubGroups() override;


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

  /// Use rock density in addition to fluid density
  localIndex m_useRockDensity;
  localIndex m_useReferencePorosity;
  localIndex m_usePorosity;

private:



};

} // namespace geos

#endif // GEOS_PHYSICSSOLVERS_GRAVITY_GRAVITYFE_COMPOSITIONALMULTUPHASEFVM_HPP_
