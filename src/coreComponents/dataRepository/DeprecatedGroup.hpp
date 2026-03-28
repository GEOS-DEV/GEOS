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
 * @file DeprecatedGroup.hpp
 *
 * This file provides utilities for marking Group classes as deprecated while maintaining
 * full functionality.
 *
 * @section deprecated_lifecycle Deprecated Stage in Group Lifecycle
 *
 * Classes in GEOS follow a four-stage lifecycle:
 *
 * 1. **Active**: Normal implementation with full functionality
 * 2. **Deprecated**: Full functionality maintained, warnings logged on use (this file)
 * 3. **Obsolete**: Type exists for catalog/schema compatibility, no functionality, errors on use
 * 4. **Removed**: Class deleted entirely when version support window expires
 *
 * Use DeprecatedGroup<T> to mark a class as deprecated while maintaining full functionality.
 * Warnings are logged to inform users to migrate to the replacement. The class continues
 * to function normally but alerts users to update their workflows.
 *
 * @see ObsoleteGroup.hpp for the next stage in the lifecycle
 */

#ifndef GEOS_DATAREPOSITORY_DEPRECATEDGROUP_HPP_
#define GEOS_DATAREPOSITORY_DEPRECATEDGROUP_HPP_

#include "Group.hpp"
#include "common/logger/Logger.hpp"
#include "common/format/Format.hpp"
#include "LvArray/src/system.hpp"
#include <string>
#include <type_traits>

namespace geos
{
namespace dataRepository
{

/**
 * @class DeprecatedGroup
 * @brief A template wrapper that injects deprecation warnings into the inheritance hierarchy.
 *
 * This class template injects itself between a derived class and its base class T.
 * When constructed, it logs a deprecation warning but maintains full functionality.
 *
 * @tparam T The base Group class to inherit from (e.g., Group, SolverBase, etc.)
 *
 * Usage pattern:
 * @code
 * // Original:
 * class MyGroup : public Group { ... };
 *
 * // After deprecation:
 * class MyGroup : public DeprecatedGroup<Group>
 * {
 * public:
 *   static string deprecationMessage() { return "MyGroup is deprecated. Use NewGroup instead."; }
 *   // ... rest of implementation unchanged
 * };
 * @endcode
 */
template< typename T >
class DeprecatedGroup : public T
{
public:

  /**
   * @brief Constructor
   * @param name The name of this Group.
   * @param parent The parent Group.
   */
  explicit DeprecatedGroup( string const & name, Group * const parent )
    : T( name, parent )
  {
    logDeprecationWarning();
  }

  /**
   * @brief Constructor
   * @param name The name of this Group.
   * @param rootNode The root node of the data repository.
   */
  explicit DeprecatedGroup( string const & name, conduit::Node & rootNode )
    : T( name, rootNode )
  {
    logDeprecationWarning();
  }

  /**
   * @brief Destructor
   */
  virtual ~DeprecatedGroup() = default;

  /// Deleted default constructor
  DeprecatedGroup() = delete;

  /// Deleted copy constructor
  DeprecatedGroup( DeprecatedGroup const & ) = delete;

  /// Move constructor
  DeprecatedGroup( DeprecatedGroup && ) = default;

  /// Deleted copy assignment
  DeprecatedGroup & operator=( DeprecatedGroup const & ) = delete;

  /// Deleted move assignment
  DeprecatedGroup & operator=( DeprecatedGroup && ) = delete;

protected:
  /**
   * @brief Override initialization hook to keep a single, consistent touch point
   *
   * Deprecated groups currently do not need extra behavior at this stage, but
   * overriding this hook keeps the lifecycle symmetric with ObsoleteGroup and
   * provides a natural place for future deprecation-time checks if needed.
   */
  virtual void initializePreSubGroups() override
  {
    T::initializePreSubGroups();
  }

private:
  /**
   * @brief Log the deprecation warning message
   */
  void logDeprecationWarning()
  {
    // Get the actual derived class type name from runtime type information
    string const className = LvArray::system::demangleType< DeprecatedGroup< T > >( );
    string const message = getDeprecationMessage();

    GEOS_LOG_RANK_0( GEOS_FMT(
                       "\n"
                       "********************************************************************************\n"
                       "* DEPRECATION WARNING\n"
                       "* Group: {}\n"
                       "* {}\n"
                       "********************************************************************************",
                       className, message ) );
  }

  /**
   * @brief Get the deprecation message
   *
   * This uses runtime polymorphism to call the derived class's deprecationMessage()
   * method if it exists, otherwise returns a default message.
   *
   * @return The deprecation message string
   */
  virtual string getDeprecationMessage() const
  {
    return "This class is deprecated and may be removed in a future version.";
  }
};

} // namespace dataRepository
} // namespace geos

#endif // GEOS_DATAREPOSITORY_DEPRECATEDGROUP_HPP_
