/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file ExecutableGroup.hpp
 */

#ifndef EXECUTABLEGROUP_HPP_
#define EXECUTABLEGROUP_HPP_

#include "common/DataTypes.hpp"
#include "Group.hpp"


/** Objects that are executable and/or able to request changes dt, should
    be derived from this type. */

namespace geosx
{

/**
 * @class ExecutableGroup
 *
 * Objects that are executable and/or able to request changes dt
 * should be derived from this type.
 */
class ExecutableGroup : public dataRepository::Group
{
public:

  using dataRepository::Group::Group;

  /*
   * If the start criteria are satisfied, then the event manager
   * will call this method
   */
  virtual void Execute( real64 const time_n,
                        real64 const dt,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        dataRepository::Group * domain ) = 0;

  /*
   * This method will inform the object that it expects to execute
   * during the next timestep
   */
  virtual void SignalToPrepareForExecution( real64 const time,
                                            real64 const dt,
                                            integer const cycle,
                                            dataRepository::Group * domain ) {}
  /*
   * This method is called as the code exits the main run loop
   */
  virtual void Cleanup( real64 const time_n,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        dataRepository::Group * domain ) {}

  /*
   * This supplies the timestep request for each executable
   * target to the event manager
   */
  virtual real64 GetTimestepRequest( real64 const time ) {return std::numeric_limits< integer >::max();}


  /// These set and supply the timestep behavior for a target
  void SetTimestepBehavior( integer behavior ){ m_timestepType = behavior; }

  integer GetTimestepBehavior(){ return m_timestepType; }


private:
  integer m_timestepType = 0;
};


} /* namespace ANST */


#endif /* EXECUTABLEGROUP_HPP_ */
