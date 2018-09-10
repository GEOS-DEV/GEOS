/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistrubute it and/or modify it under
 * the terms of the GNU Lesser General Public Liscense (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file ExecutableGroup.hpp
 */

#ifndef EXECUTABLEGROUP_HPP_
#define EXECUTABLEGROUP_HPP_

#include "ManagedGroup.hpp"
#include "common/DataTypes.hpp"


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
class ExecutableGroup : public dataRepository::ManagedGroup
{
public:

  using dataRepository::ManagedGroup::ManagedGroup;

  /*
   * If the start criteria are satisfied, then the event manager
   * will call this method
   */
  virtual void Execute( real64 const & time_n,
                        real64 const & dt,
                        int const cycleNumber,
                        real64 const & eventPosition,
                        dataRepository::ManagedGroup * domain ) = 0;

  /*
   * This method will inform the object that it expects to execute
   * during the next timestep
   */
  virtual void SignalToPrepareForExecution(real64 const time,
                                           real64 const dt,  
                                           integer const cycle,
                                           dataRepository::ManagedGroup * domain) {}
  /*
   * This method is called as the code exits the main run loop
   */
  virtual void Cleanup( real64 const & time_n,
                        int const cycleNumber,
                        real64 const & eventPosition,
                        dataRepository::ManagedGroup * domain ) {}

  /*
   * This supplies the timestep request for each executable
   * target to the event manager
   */
  virtual real64 GetTimestepRequest(real64 const time) {return std::numeric_limits<integer>::max();}
};


} /* namespace ANST */


#endif /* EXECUTABLEGROUP_HPP_ */
