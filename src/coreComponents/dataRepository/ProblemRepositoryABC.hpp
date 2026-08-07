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
 * @file ProblemRepository.hpp
 */

#ifndef GEOS_DATAREPOSITORY_PROBLEMREPOSITORYABC_HPP_
#define GEOS_DATAREPOSITORY_PROBLEMREPOSITORYABC_HPP_

#include "common/DataTypes.hpp"

namespace geos
{
namespace dataRepository
{

/**
 * @brief Interface for an object which gives access all the problem data-repository.
 *        No dependancy on Group since Group needs some features (getNoProblemPath()).
 *        - No exposure of the problem GroupKey string, to avoid "invisible dependancies":
 *       being dependent of the existence
 *        of a Group at a specific data-repository path without mentionning its type (making it
 *        hard to see/find the dependency).
 */
class ProblemRepositoryABC
{
public:

  /**
   * @brief Utility function to output a Group/Wrapper path without "Problem": "Problem/commandLine" -> "/commandLine"
   * @param originalPath Group/Wrapper conduit node path
   * @return string without "Problem" at start
   */
  static string getNoProblemPath( string const & originalPath );

protected:

  /**
   * @brief Group keys for faster lookup.
   *        Remain internal data to avoid "invisible dependancies": being dependent of the existence
   *        of a Group at a specific data-repository path without mentionning its type (making it
   *        hard to see/find the dependency).
   */
  struct ProblemGroupKeys
  {
    static constexpr char const * problemString() { return "Problem"; }
  };

  ProblemRepositoryABC()
  {}

};

} /* namespace dataRepository */
} /* namespace geos */

#endif /* GEOS_DATAREPOSITORY_PROBLEMREPOSITORYABC_HPP_ */
