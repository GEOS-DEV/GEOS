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

#ifndef EXECUTABLEGROUP_HPP_
#define EXECUTABLEGROUP_HPP_

#include "ManagedGroup.hpp"
#include "common/DataTypes.hpp"


/** Objects that are executable and/or able to request changes dt, should
    be derived from this type. */

namespace geosx
{

class ExecutableGroup : public dataRepository::ManagedGroup
{
public:

  using dataRepository::ManagedGroup::ManagedGroup;

  /*
   * If the start criteria are satisfied, then the event manager will call this method
   */
  virtual void Execute( real64 const & time_n,
                        real64 const & dt,
                        int const cycleNumber,
                        dataRepository::ManagedGroup * domain ) = 0;

  /*
   * This supplies the timestep request for each executable target to the event manager
   */
  virtual real64 GetTimestepRequest(real64 const time) {return std::numeric_limits<integer>::max();}
};


} /* namespace ANST */


#endif /* EXECUTABLEGROUP_HPP_ */
