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
 * @file FieldStatisticsBase.hpp
 */

#ifndef SRC_CORECOMPONENTS_PHYSICSSOLVERS_FIELDSTATISTICSBASE_HPP_
#define SRC_CORECOMPONENTS_PHYSICSSOLVERS_FIELDSTATISTICSBASE_HPP_

#include "events/tasks/TaskBase.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mesh/MeshLevel.hpp"

namespace geosx
{

/**
 * @class FieldStatisticsBase
 *
 * Base task class allowing for the computation of aggregate statistics
 */
template< typename SOLVER >
class FieldStatisticsBase : public TaskBase
{
public:

  /**
   * @brief Constructor for the state reset class
   * @param[in] name the name of the task coming from the xml
   * @param[in] parent the parent group of the task
   */
  FieldStatisticsBase( const string & name,
                       Group * const parent )
    : TaskBase( name, parent )
  {
    enableLogLevelInput();

    string const key = SOLVER::coupledSolverAttributePrefix() + "SolverName";
    registerWrapper( key, &m_solverName ).
      setInputFlag( dataRepository::InputFlags::REQUIRED ).
      setDescription( "Name of the " + SOLVER::coupledSolverAttributePrefix() + " solver" );
  }

  /**
   * @defgroup Tasks Interface Functions
   *
   * This function implements the interface defined by the abstract TaskBase class
   */
  /**@{*/

  virtual bool execute( real64 const time_n,
                        real64 const dt,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        DomainPartition & domain ) override = 0;

  /**@}*/

protected:

  void postProcessInput() override
  {
    ProblemManager & problemManager = this->getGroupByPath< ProblemManager >( "/Problem" );
    PhysicsSolverManager & physicsSolverManager = problemManager.getPhysicsSolverManager();

    m_solver = physicsSolverManager.getGroupPointer< SOLVER >( m_solverName );
    GEOSX_THROW_IF( m_solver == nullptr,
                    GEOSX_FMT( "Could not find solver '{}' of type {}",
                               m_solverName, LvArray::system::demangleType< SOLVER >() ),
                    InputError );
  }

  /// Name of the solver
  string m_solverName;

  /// Pointer to the physics solver
  SOLVER * m_solver;

};


} /* namespace geosx */

#endif /* SRC_CORECOMPONENTS_PHYSICSSOLVERS_FIELDSTATISTICSBASE_HPP_ */
