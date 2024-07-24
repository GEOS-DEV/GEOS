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

#include "ResidualTask.hpp"

namespace geos
{

ResidualTask::ResidualTask( const geos::string & name, geos::dataRepository::Group *parent )
  : TaskBase( name, parent ),
  m_solverTarget( nullptr ),
  m_outputTarget( nullptr )
{

  registerWrapper( viewKeyStruct::outputTargetString(), &m_outputTargetString ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Name of the output target." );

  registerWrapper( viewKeyStruct::solverTargetString(), &m_solverTargetString ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Name of the solver target." );

}

void ResidualTask::initializePostSubGroups()
{

  TaskBase::initializePostSubGroups();

  if( !m_solverTargetString.empty() )
  {
    try
    {
      m_solverTarget = &this->getGroupByPath< ExecutableGroup >( m_solverTargetString );
    }
    catch( std::exception const & e )
    {
      throw InputError( e, GEOS_FMT( "Error while reading {}:\n",
                                     getWrapperDataContext( viewKeyStruct::solverTargetString() ) ) );
    }
  }
  if( !m_outputTargetString.empty() )
  {
    try
    {
      m_outputTarget = &this->getGroupByPath< ExecutableGroup >( m_outputTargetString );
    }
    catch( std::exception const & e )
    {
      throw InputError( e, GEOS_FMT( "Error while reading {}:\n",
                                     getWrapperDataContext( viewKeyStruct::outputTargetString() ) ) );
    }
  }
}

bool ResidualTask::execute( const geos::real64 time_n, const geos::real64 dt, const geos::integer cycleNumber,
                            const geos::integer eventCounter, const geos::real64 eventProgress,
                            geos::DomainPartition & domain )
{

  if( m_solverTarget!= nullptr && m_outputTarget!= nullptr )
  {
    m_solverTarget->execute( time_n, dt, cycleNumber, eventCounter, eventProgress, domain );
    if( dynamicCast< SolverBase * >( m_solverTarget )->hasNonlinearIssues() )
      m_outputTarget->execute( time_n, dt, cycleNumber, eventCounter, eventProgress, domain );

  }

  return false;
}

REGISTER_CATALOG_ENTRY( TaskBase,
                        ResidualTask,
                        string const &, dataRepository::Group * const )

}//geos
