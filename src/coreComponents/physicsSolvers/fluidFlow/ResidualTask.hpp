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

#ifndef GEOSX_RESIDUALTASK_HPP
#define GEOSX_RESIDUALTASK_HPP

#include "events/tasks/TaskBase.hpp"
#include "FlowSolverBase.hpp"
#include "mesh/DomainPartition.hpp"

namespace geos
{

using namespace dataRepository;


//template< typename SOLVER >
class ResidualTask : public TaskBase
{

public:
  ResidualTask( string const & name, Group * parent );
  virtual ~ResidualTask() override {};

  static string catalogName() { return "Residual"; }

  void initializePostSubGroups() override;

  bool execute( const geos::real64 time_n,
                const geos::real64 dt,
                const geos::integer cycleNumber,
                const geos::integer eventCounter,
                const geos::real64 eventProgress,
                geos::DomainPartition & domain ) override;

private:

  struct viewKeyStruct
  {
    constexpr static char const * solverTargetString()
    {return "solverTarget";}

    constexpr static char const * outputTargetString()
    { return "outputTarget";}
  };

  string m_solverTargetString, m_outputTargetString;
  ExecutableGroup * m_solverTarget;
  ExecutableGroup * m_outputTarget;

  //TODO
//  void registerDataOnMesh(geos::dataRepository::Group &meshBodies) override;

};

}//geos

#endif //GEOSX_RESIDUALTASK_HPP
