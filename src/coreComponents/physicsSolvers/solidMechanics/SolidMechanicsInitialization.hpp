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
 * @file SolidMechanicsInitialization.hpp
 */

#ifndef SRC_CORECOMPONENTS_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSINITIALIZATION_HPP_
#define SRC_CORECOMPONENTS_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSINITIALIZATION_HPP_

#include "events/tasks/TaskBase.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsStateReset.hpp"

namespace geos
{
class SolidMechanicsStatistics;

/**
 * @class SolidMechanicsInitialization
 *
 * Task class used to perform stress initialization in solid mechanics-only simulations.
 */
template< typename SOLID_SOLVER >
class SolidMechanicsInitialization : public TaskBase
{
public:

  SolidMechanicsInitialization( const string & name,
                                Group * const parent );

  ~SolidMechanicsInitialization() override;

  static string catalogName()
  {
    return SOLID_SOLVER::catalogName() + "Initialization";
  }

  virtual bool execute( real64 const time_n,
                        real64 const dt,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        DomainPartition & domain ) override;

private:

  struct viewKeyStruct
  {
    constexpr static char const * solidSolverNameString() { return "solidSolverName"; }
    constexpr static char const * solidMechanicsStatisticsNameString() { return "solidMechanicsStatisticsName"; }
  };

  void postInputInitialization() override;

  string m_solidSolverName;
  string m_solidMechanicsStatisticsName;

  SOLID_SOLVER * m_solidSolver;
  SolidMechanicsStatistics * m_solidMechanicsStatistics;

  SolidMechanicsStateReset m_solidMechanicsStateResetTask;
};

}

#endif /* SRC_CORECOMPONENTS_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSINITIALIZATION_HPP_ */
