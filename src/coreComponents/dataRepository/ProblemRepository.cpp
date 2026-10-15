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
 * @file MeshLevel.cpp
 */

#include "ProblemRepository.hpp"
#include "codingUtilities/RTTypes.hpp"

namespace geos
{
namespace dataRepository
{

/**
 * @name Inline functions implementation
 */
///@{

string ProblemRepositoryABC::getNoProblemPath( string const & originalPath )
{
  // In the Conduit node hierarchy everything begins with 'Problem', we should change it so that
  // the ProblemManager actually uses the root Conduit Node but that will require a full rebaseline.
  std::string_view const path( originalPath );
  std::string_view const problemPrefix( ProblemGroupKeys::problemString() );

  if( !stringutilities::startsWith( path, problemPrefix ) )
    return string( path );

  string const noProblem = path.substr( problemPrefix.size() );
  return noProblem.empty() ? "/" : noProblem;
}

ProblemRepository & ProblemRepository::get( Group & group )
{
  // the ProblemRepository is expected to always be the root Group instance.
  Group * current = &group;
  while( current->hasParent() )
  {
    current = &current->getParent();
  }

  // if wrong/unknown type requested, we will have a GEOS_ERROR here
  ProblemRepository & root = dynamicCast< ProblemRepository & >( *current );
  return root;
}

ProblemRepository const & ProblemRepository::get( Group const & group )
{
  // the ProblemRepository is expected to always be the root Group instance.
  Group const * current = &group;
  while( current->hasParent() )
  {
    current = &current->getParent();
  }

  // if wrong/unknown type requested, we will have a GEOS_ERROR here
  ProblemRepository const & root = dynamicCast< ProblemRepository const & >( *current );
  return root;
}

} /* namespace dataRepository */
} /* namespace geos */
