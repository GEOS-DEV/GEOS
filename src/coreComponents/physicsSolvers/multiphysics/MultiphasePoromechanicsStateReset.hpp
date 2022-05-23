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
 * @file MultiphasePoromechanicsStateReset.hpp
 */

#ifndef SRC_CORECOMPONENTS_PHYSICSSOLVERS_MULTIPHYSICS_MULTIPHASEPOROMECHANICSSTATERESET_HPP_
#define SRC_CORECOMPONENTS_PHYSICSSOLVERS_MULTIPHYSICS_MULTIPHASEPOROMECHANICSSTATERESET_HPP_

#include "events/tasks/TaskBase.hpp"

namespace geosx
{

class MultiphasePoromechanicsSolver;

/**
 * @class MultiphasePoromechanicsStateReset
 *
 * Task class allowing for various "reset" actions on solver configuration and state variables,
 * useful for pre- and post-equilibration actions. Currently, this class gives the user the option
 * to reset change the MGR recipe during the initialization step. A common
 * workflow might be:
 * (1) useInitializationSolverConfig = true
 * (2) solve an equilibration problem for initial conditions,
 * (3) useInitializationSolverConfig = false
 * (4) run a normal simulation
 *
 */
class MultiphasePoromechanicsStateReset : public TaskBase
{
public:

  /**
   * @brief Constructor for the state reset class
   * @param[in] name the name of the task coming from the xml
   * @param[in] parent the parent group of the task
   */
  MultiphasePoromechanicsStateReset( const string & name,
                                     Group * const parent );

  /// Destructor for the class
  ~MultiphasePoromechanicsStateReset() override;

  /// Accessor for the catalog name
  static string catalogName() { return "MultiphasePoromechanicsStateReset"; }

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
    /// String for the multiphase poromechanics solver name
    constexpr static char const * multiphasePoromechanicsSolverNameString() { return "multiphasePoromechanicsSolverName"; }
    /// String for the solver configuration
    constexpr static char const * useInitializationSolverConfigurationString() { return "useInitializationSolverConfiguration"; }
  };


  void postProcessInput() override;

  /// Name of the multiphase poromechanics solver
  string m_multiphasePoromechanicsSolverName;

  /// Pointer to the multiphase poromechanics solver
  MultiphasePoromechanicsSolver * m_multiphasePoromechanicsSolver;

  /// Flag for reset the solver configuration
  integer m_useInitializationSolverConfiguration;

};


} /* namespace geosx */

#endif /* SRC_CORECOMPONENTS_PHYSICSSOLVERS_MULTIPHYSICS_MULTIPHASEPOROMECHANICSSTATERESET_HPP_ */
