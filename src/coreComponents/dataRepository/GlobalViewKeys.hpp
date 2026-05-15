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
 * @file GlobalViewKeys.hpp
 */

#ifndef GEOS_DATAREPOSITORY_GLOBALVIEWKEYS_HPP_
#define GEOS_DATAREPOSITORY_GLOBALVIEWKEYS_HPP_


namespace geos
{
namespace dataRepository
{


/**
 * @struct GlobalViewKeys
 */
struct GlobalViewKeys
{
  static constexpr char const * problem()                   { return "Problem"; }
  static constexpr char const * commandLine()               { return "commandLine"; }
  static constexpr char const * domain()                    { return "domain"; }
  static constexpr char const * constitutiveManager()       { return "Constitutive"; }
  static constexpr char const * eventManager()              { return "Events"; }
  static constexpr char const * externalDataSourceManager() { return "ExternalDataSource"; }
  static constexpr char const * fieldSpecificationManager() { return "FieldSpecifications"; }
  static constexpr char const * functionManager()           { return "Functions"; }
  static constexpr char const * geometricObjectManager()    { return "Geometry"; }
  static constexpr char const * meshManager()               { return "Mesh"; }
  static constexpr char const * numericalMethodsManager()   { return "NumericalMethods"; }
  static constexpr char const * outputManager()             { return "Outputs"; }
  static constexpr char const * physicSolverManager()       { return "Solvers"; }
  static constexpr char const * taskManager()               { return "Tasks"; }
  static constexpr char const * cellManager()               { return "cellManager"; }
  static constexpr char const * particleManager()           { return "particleManager"; }
  static constexpr char const * partitionManager()          { return "partitionManager"; }
};

} /* namespace dataRepository */
} /* namespace geos */


#endif /* GEOS_DATAREPOSITORY_GLOBALVIEWKEYS_HPP_ */
