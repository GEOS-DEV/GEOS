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
  /// @return Root problem group name
  static constexpr char const * problem()                   { return "Problem"; }
  /// @return Command-line group name
  static constexpr char const * commandLine()               { return "commandLine"; }
  /// @return Domain partition group name
  static constexpr char const * domain()                    { return "domain"; }
  /// @return Constitutive manager group name
  static constexpr char const * constitutiveManager()       { return "Constitutive"; }
  /// @return Event manager group name
  static constexpr char const * eventManager()              { return "Events"; }
  /// @return External data source manager group name
  static constexpr char const * externalDataSourceManager() { return "ExternalDataSource"; }
  /// @return FieldSpecification manager group name
  static constexpr char const * fieldSpecificationManager() { return "FieldSpecifications"; }
  /// @return Function manager group name
  static constexpr char const * functionManager()           { return "Functions"; }
  /// @return Geometric object manager group name
  static constexpr char const * geometricObjectManager()    { return "Geometry"; }
  /// @return Mesh manager group name
  static constexpr char const * meshManager()               { return "Mesh"; }
  /// @return Numerical methods manager group name
  static constexpr char const * numericalMethodsManager()   { return "NumericalMethods"; }
  /// @return Outputs manager group name
  static constexpr char const * outputManager()             { return "Outputs"; }
  /// @return Physics solvers manager group name
  static constexpr char const * physicsSolverManager()      { return "Solvers"; }
  /// @return Tasks manager group name
  static constexpr char const * tasksManager()              { return "Tasks"; }
};

} /* namespace dataRepository */
} /* namespace geos */


#endif /* GEOS_DATAREPOSITORY_GLOBALVIEWKEYS_HPP_ */
