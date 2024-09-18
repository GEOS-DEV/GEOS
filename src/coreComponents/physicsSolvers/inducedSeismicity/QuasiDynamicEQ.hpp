/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#ifndef GEOS_PHYSICSSOLVERS_INDUCED_QUASIDYNAMICEQ_HPP
#define GEOS_PHYSICSSOLVERS_INDUCED_QUASIDYNAMICEQ_HPP

#include "physicsSolvers/SolverBase.hpp"

namespace geos
{

class QuasiDynamicEQ : public SolverBase
{
public:
  /// The default nullary constructor is disabled to avoid compiler auto-generation:
  QuasiDynamicEQ() = delete;

  /// The constructor needs a user-defined "name" and a parent Group (to place this instance in the tree structure of classes)
  QuasiDynamicEQ( const string & name,
                  Group * const parent );

  /// Destructor
  virtual ~QuasiDynamicEQ() override;

  static string catalogName() { return "QuasiDynamicEQ"; }

  /**
   * @return Get the final class Catalog name
   */
  virtual string getCatalogName() const override { return catalogName(); }

  /// This method ties properties with their supporting mesh
  virtual void registerDataOnMesh( Group & meshBodies ) override;

  struct viewKeyStruct : public SolverBase::viewKeyStruct
  {
    /// stress solver name
    static constexpr char const * stressSolverNameString() { return "stressSolverName"; }
    /// Friction law name string
    constexpr static char const * frictionLawNameString() { return "frictionLawName"; }
    /// Friction law name string
    constexpr static char const * shearImpedanceString() { return "shearImpedance"; }
    /// max number of Newton iterations string
    constexpr static char const * maxNumberOfNewtonIterationsString() { return "maxNumberOfNewtonIterations"; }
  };

  virtual real64 solverStep( real64 const & time_n,
                             real64 const & dt,
                             integer const cycleNumber,
                             DomainPartition & domain ) override final;

private:

  virtual void postInputInitialization() override;

  real64 updateStresses( real64 const & time_n,
                         real64 const & dt,
                         const int cycleNumber,
                         DomainPartition & domain ) const;

  /**
   * @brief save the old state
   * @param subRegion
   */
  void saveOldState( ElementSubRegionBase & subRegion ) const;

  /// pointer to stress solver
  SolverBase * m_stressSolver;

  /// stress solver name
  string m_stressSolverName;

  /// max number of newton iterations for rate and state solver
  integer m_maxNewtonIterations;

  /// Shear impedance
  real64 m_shearImpedance;
};

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_INDUCED_QUASIDYNAMICEQ_HPP */
