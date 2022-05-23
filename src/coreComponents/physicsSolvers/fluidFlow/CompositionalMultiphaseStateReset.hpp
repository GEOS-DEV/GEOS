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
 * @file CompositionalMultiphaseStateReset.hpp
 */

#ifndef SRC_CORECOMPONENTS_PHYSICSSOLVERS_MULTIPHYSICS_COMPOSITIONALMULTIPHASESTATERESET_HPP_
#define SRC_CORECOMPONENTS_PHYSICSSOLVERS_MULTIPHYSICS_COMPOSITIONALMULTIPHASESTATERESET_HPP_

#include "events/tasks/TaskBase.hpp"

namespace geosx
{

class CompositionalMultiphaseBase;

/**
 * @class CompositionalMultiphaseStateReset
 *
 * Task class allowing for various "reset" actions on solver configuration and state variables,
 * useful for pre- and post-equilibration actions. Currently, this class gives the user the option
 * to reset change the MGR recipe during the initialization step. A common
 * workflow might be:
 * (1) freezeFlowVariablesDuringStep = true
 * (2) solve an equilibration problem for initial conditions,
 * (3) freezeFlowVariablesDuringStep = false
 * (4) run a normal simulation
 *
 */
class CompositionalMultiphaseStateReset : public TaskBase
{
public:

  /**
   * @brief Constructor for the state reset class
   * @param[in] name the name of the task coming from the xml
   * @param[in] parent the parent group of the task
   */
  CompositionalMultiphaseStateReset( const string & name,
                                     Group * const parent );

  /// Destructor for the class
  ~CompositionalMultiphaseStateReset() override;

  /// Accessor for the catalog name
  static string catalogName() { return "CompositionalMultiphaseStateReset"; }

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
                        DomainPartition & domain ) override;

  /**@}*/

private:

  /**
   * @struct viewKeyStruct holds char strings and viewKeys for fast lookup
   */
  struct viewKeyStruct
  {
    /// String for the flow solver name
    constexpr static char const * flowSolverNameString() { return "flowSolverName"; }
    /// String for the solver configuration
    constexpr static char const * freezeFlowVariablesDuringStepString() { return "freezeFlowVariablesDuringStep"; }
  };


  void postProcessInput() override;

  /// Name of the flow solver
  string m_flowSolverName;

  /// Pointer to the flow solver
  CompositionalMultiphaseBase * m_flowSolver;

  /// Flag for freeze the flow variables during step
  integer m_freezeFlowVariablesDuringStep;

};


} /* namespace geosx */

#endif /* SRC_CORECOMPONENTS_PHYSICSSOLVERS_MULTIPHYSICS_COMPOSITIONALMULTIPHASESTATERESET_HPP_ */
