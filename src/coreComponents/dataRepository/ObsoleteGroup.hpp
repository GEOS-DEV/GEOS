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
 * @file ObsoleteGroup.hpp
 *
 * This file provides utilities for marking Group classes as obsolete while maintaining
 * catalog and schema compatibility for input parsing.
 *
 * @section obsolete_lifecycle Obsolete Stage in Group Lifecycle
 *
 * Use ObsoleteGroup<T> when removing functionality but maintaining type for input compatibility:
 * - Construction succeeds (for XML parsing and catalog registration)
 * - All wrappers use NoopWrapper (zero memory footprint)
 * - Initialization/execution throws clear error messages and is considered a programming error
 * - No implementation beyond constructors required
 * - Catalog entry preserved for schema validation
 *
 * This allows old input decks to parse without immediate failure while providing clear
 * error messages if the obsolete functionality is actually invoked. Any attempt to
 * "use" an obsolete group (beyond catalog/schema queries) is treated as a logic error
 * in the calling code and is expected to terminate execution.
 *
 * @see DeprecatedGroup.hpp for the full lifecycle workflow documentation
 */

#ifndef GEOS_DATAREPOSITORY_OBSOLETEGROUP_HPP_
#define GEOS_DATAREPOSITORY_OBSOLETEGROUP_HPP_

#include "Group.hpp"
#include "Wrapper.hpp"
#include "common/logger/Logger.hpp"
#include "common/format/Format.hpp"
#include "codingUtilities/RTTypes.hpp"
#include "LvArray/src/system.hpp"
#include <cstdlib>
#include <string>
#include <type_traits>

namespace geos
{
namespace dataRepository
{

/**
 * @class NoopWrapper
 * @brief A wrapper shim that has type information but no memory footprint
 *
 * This class mimics the WrapperBase / Wrapper<T> surface API but does not actually
 * store any data. It is used by ObsoleteGroup to maintain type-ness and catalog
 * compatibility without memory overhead.
 *
 * Semantics:
 * - NoopWrapper is intended for schema/catalog use only (e.g. type discovery).
 * - It must not be used for reading or writing simulation data.
 * - All operations are implemented as no-ops and size/capacity/bytesAllocated are 0.
 *
 * Any code path that attempts to treat an obsolete group's wrappers as real
 * storage (e.g. calling Wrapper<T>::reference(), resize to hold state, etc.)
 * is a programming error in the caller. Tests intentionally restrict themselves
 * to WrapperBase-level queries when interacting with NoopWrapper.
 *
 * @tparam T The type that would normally be wrapped
 */
template< typename T >
class NoopWrapper final : public WrapperBase
{
public:
  /// Alias for the wrapped type
  using TYPE = T;

  /**
   * @brief Constructor
   * @param name name of the object
   * @param parent parent group which owns the Wrapper
   */
  explicit NoopWrapper( string const & name, Group & parent )
    : WrapperBase( name, parent, rtTypes::getTypeName( typeid( T ) ) )
  {
    // Mark as optional to avoid any required/default value checks
    this->setInputFlag( InputFlags::OPTIONAL );
  }

  /// @brief Destructor
  virtual ~NoopWrapper() noexcept override = default;

  /// @copydoc WrapperBase::size
  virtual localIndex size() const override { return 0; }

  /// @copydoc WrapperBase::voidPointer
  virtual void const * voidPointer() const override { return nullptr; }

  /// @copydoc WrapperBase::elementByteSize
  virtual localIndex elementByteSize() const override { return sizeof( T ); }

  /// @copydoc WrapperBase::bytesAllocated
  virtual size_t bytesAllocated() const override { return 0; }

  /// @copydoc WrapperBase::resize(int, localIndex const * const)
  virtual void resize( int, localIndex const * const ) override {}

  /// @copydoc WrapperBase::reserve
  virtual void reserve( localIndex const ) override {}

  /// @copydoc WrapperBase::capacity
  virtual localIndex capacity() const override { return 0; }

  /// @copydoc WrapperBase::resize(localIndex)
  virtual void resize( localIndex ) override {}

  /// @copydoc WrapperBase::copy
  virtual void copy( localIndex const, localIndex const ) override {}

  /// @copydoc WrapperBase::erase
  virtual void erase( std::set< localIndex > const & ) override {}

  /// @copydoc WrapperBase::move
  virtual void move( LvArray::MemorySpace const, bool const ) const override {}

  /// @copydoc WrapperBase::getTypeRegex
  virtual Regex const & getTypeRegex() const override
  {
    return rtTypes::getTypeRegex< T >();
  }

  /// @copydoc WrapperBase::getTypeId
  virtual const std::type_info & getTypeId() const noexcept override
  {
    return typeid( T );
  }

  /// @copydoc WrapperBase::hasDefaultValue
  virtual bool hasDefaultValue() const override { return false; }

  /// @copydoc WrapperBase::getDefaultValueString
  virtual string getDefaultValueString() const override { return ""; }

  /// @copydoc WrapperBase::clone
  virtual std::unique_ptr< WrapperBase > clone( string const & name, Group & parent ) override
  {
    return std::make_unique< NoopWrapper< T > >( name, parent );
  }

  /// @copydoc WrapperBase::copyWrapper
  virtual void copyWrapper( WrapperBase const & ) override {}

  /// @copydoc WrapperBase::copyWrapperAttributes
  virtual void copyWrapperAttributes( WrapperBase const & source ) override
  {
    WrapperBase::copyWrapperAttributes( source );
  }

  /// @copydoc WrapperBase::numArrayDims
  virtual int numArrayDims() const override { return 0; }

  /// @copydoc WrapperBase::numArrayComp
  virtual localIndex numArrayComp() const override { return 0; }

  /// @copydoc WrapperBase::setDimLabels
  virtual WrapperBase & setDimLabels( integer const, Span< string const > const ) override
  {
    return *this;
  }

  /// @copydoc WrapperBase::getDimLabels
  virtual Span< string const > getDimLabels( integer const ) const override
  {
    return {};
  }

  /// @copydoc WrapperBase::getHistoryMetadata
  virtual HistoryMetadata getHistoryMetadata( localIndex const = -1 ) const override final
  {
    return {};
  }

  /// @copydoc WrapperBase::isPackable
  virtual bool isPackable( bool ) const override { return false; }

  /// @copydoc WrapperBase::unpack
  virtual localIndex unpack( buffer_unit_type const * &, bool, bool, parallelDeviceEvents & ) override final
  {
    return 0;
  }

  /// @copydoc WrapperBase::unpackByIndex
  virtual localIndex unpackByIndex( buffer_unit_type const * &, arrayView1d< localIndex const > const &,
                                     bool, bool, parallelDeviceEvents &, MPI_Op = MPI_REPLACE ) override
  {
    return 0;
  }

  /// @copydoc WrapperBase::processInputFile
  virtual bool processInputFile( xmlWrapper::xmlNode const &, xmlWrapper::xmlNodePos const & ) override
  {
    return false;
  }

  /// @copydoc WrapperBase::addBlueprintField
  virtual void addBlueprintField( conduit::Node &, string const &, string const &,
                                  stdVector< string > const & = {} ) const override {}

  /// @copydoc WrapperBase::populateMCArray
  virtual void populateMCArray( conduit::Node &, stdVector< string > const & = {} ) const override {}

  /// @copydoc WrapperBase::averageOverSecondDim
  virtual std::unique_ptr< WrapperBase > averageOverSecondDim( string const & name, Group & parent ) const override
  {
    return std::make_unique< NoopWrapper< T > >( name, parent );
  }

  /// @copydoc WrapperBase::registerToWrite
  virtual void registerToWrite() const override {}

  /// @copydoc WrapperBase::finishWriting
  virtual void finishWriting() const override {}

  /// @copydoc WrapperBase::loadFromConduit
  virtual bool loadFromConduit() override { return false; }

  /// @copydoc WrapperBase::copyData
  virtual void copyData( WrapperBase const & ) override {}

#if defined(USE_TOTALVIEW_OUTPUT)
  /// @copydoc WrapperBase::totalviewTypeName
  virtual string totalviewTypeName() const override
  {
    return LvArray::system::demangleType< T >();
  }
#endif

  /// @copydoc WrapperBase::createPythonObject
  virtual PyObject * createPythonObject() override
  {
    return nullptr;
  }

private:
  /// @copydoc WrapperBase::packPrivate
  virtual localIndex packPrivate( buffer_unit_type * &, bool, bool, parallelDeviceEvents & ) const override
  {
    return 0;
  }

  /// @copydoc WrapperBase::packSizePrivate
  virtual localIndex packSizePrivate( bool, bool, parallelDeviceEvents & ) const override
  {
    return 0;
  }

  /// @copydoc WrapperBase::packByIndexPrivate
  virtual localIndex packByIndexPrivate( buffer_unit_type * &, arrayView1d< localIndex const > const &,
                                         bool, bool, parallelDeviceEvents & ) const override
  {
    return 0;
  }

  /// @copydoc WrapperBase::packByIndexSizePrivate
  virtual localIndex packByIndexSizePrivate( arrayView1d< localIndex const > const &, bool, bool,
                                             parallelDeviceEvents & ) const override
  {
    return 0;
  }
};

/**
 * @class ObsoleteGroup
 * @brief A template wrapper that injects obsolete behavior into the inheritance hierarchy.
 *
 * This class template injects itself between a derived class and its base class T.
 * The class type remains registered in the catalog for parsing compatibility, but
 * any attempt to actually use the functionality produces a fatal error.
 *
 * Key behaviors:
 * - Class remains in catalog for schema/parsing compatibility
 * - Construction succeeds (for parsing), but initialization/execution fails with clear error
 * - All wrapper registrations use NoopWrapper (no memory allocation)
 * - No functional implementation is provided or required
 * - Type information preserved for input validation
 *
 * Lifecycle workflow:
 * 1. Active class: Normal Group implementation
 * 2. Deprecated: Inherit from DeprecatedGroup<BaseClass>, full functionality, logs warnings
 * 3. Obsolete: Inherit from ObsoleteGroup<BaseClass>, no functionality, errors on use
 * 4. Removed: Delete class entirely when version support window expires
 *
 * @tparam T The base Group class to inherit from (e.g., Group, SolverBase, etc.)
 *
 * Usage pattern:
 * @code
 * // Step 1: Original active class
 * class OldSolver : public SolverBase
 * {
 * public:
 *   OldSolver( string const & name, Group * parent );
 *   virtual bool execute(...) override { ...implementation... }
 * };
 *
 * // Step 2: Mark as deprecated (transition period with warnings)
 * class OldSolver : public DeprecatedGroup<SolverBase>
 * {
 * public:
 *   static string deprecationMessage() { return "Use NewSolver instead"; }
 *   OldSolver( string const & name, Group * parent );
 *   virtual bool execute(...) override { ...implementation still works... }
 * };
 *
 * // Step 3: Mark as obsolete (type exists, no functionality)
 * class OldSolver : public ObsoleteGroup<SolverBase>
 * {
 * public:
 *   static string obsoleteMessage() { return "Removed in v2.0. Use NewSolver instead."; }
 *   // Only constructors needed, no other implementation
 * };
 *
 * // Step 4: Eventually remove the class definition entirely
 * @endcode
 */
template< typename T >
class ObsoleteGroup : public T
{
public:

  /**
   * @brief Constructor - allows object creation for catalog/parsing but logs warning
   * @param name The name of this Group.
   * @param parent The parent Group.
   */
  explicit ObsoleteGroup( string const & name, Group * const parent )
    : T( name, parent )
  {
    logObsoleteWarning();
  }

  /**
   * @brief Constructor - allows object creation for catalog/parsing but logs warning
   * @param name The name of this Group.
   * @param rootNode The root node of the data repository.
   */
  explicit ObsoleteGroup( string const & name, conduit::Node & rootNode )
    : T( name, rootNode )
  {
    logObsoleteWarning();
  }

  /**
   * @brief Destructor
   */
  virtual ~ObsoleteGroup() = default;

  /// Deleted default constructor
  ObsoleteGroup() = delete;

  /// Deleted copy constructor
  ObsoleteGroup( ObsoleteGroup const & ) = delete;

  /// Move constructor
  ObsoleteGroup( ObsoleteGroup && ) = default;

  /// Deleted copy assignment
  ObsoleteGroup & operator=( ObsoleteGroup const & ) = delete;

  /// Deleted move assignment
  ObsoleteGroup & operator=( ObsoleteGroup && ) = delete;

  /**
   * @brief Register a no-op wrapper instead of a real one
   * @tparam U The type to wrap
   * @param name The name of the wrapper
   * @return Reference to a NoopWrapper cast to Wrapper for API compatibility
   *
   * All wrapper registrations in obsolete groups use NoopWrapper to avoid memory
   * allocation for unused functionality while maintaining API compatibility.
   *
   * IMPORTANT:
   * - The returned reference is only valid to use through the WrapperBase API
   *   (e.g. size(), bytesAllocated(), getTypeId(), etc.).
   * - Calling Wrapper<T>-specific data accessors (such as reference(), setDefaultValue(),
   *   or methods that imply real storage) on this reference is undefined behavior and
   *   considered a programming error.
   * - Obsolete groups are for catalog/schema compatibility only; callers must not
   *   attempt to read or write simulation data through their wrappers.
   */
  template< typename U >
  Wrapper< U > & registerWrapper( string const & name )
  {
    // Create a NoopWrapper instead of a real Wrapper
    auto noopWrapper = std::make_unique< NoopWrapper< U > >( name, *this );
    WrapperBase & base = T::registerWrapper( std::move( noopWrapper ) );
    // Cast for API compatibility - the wrapper will never be used functionally
    return *reinterpret_cast< Wrapper< U > * >( &base );
  }

  /**
   * @brief Register a no-op wrapper with explicit default value (value is ignored)
   * @tparam U The type to wrap
   * @param name The name of the wrapper
   * @return Reference to a NoopWrapper cast to Wrapper for API compatibility
   */
  template< typename U >
  Wrapper< U > & registerWrapper( string const & name, U const & )
  {
    return registerWrapper< U >( name );
  }

protected:
  /**
   * @brief Override initialization to throw error - obsolete groups cannot be configured
   *
   * This hook is called on every group in the runtime tree after parsing.
   * For obsolete groups, any attempt to reach this point is treated as a
   * programming error (the group should never participate in configuration
   * or execution), so we unconditionally trigger a hard failure.
   */
  virtual void initializePreSubGroups() override
  {
    logObsoleteError();
  }

private:
  /**
   * @brief Log a warning when obsolete group is constructed (for parsing)
   */
  void logObsoleteWarning() const
  {
    // Get the actual derived class type name from runtime type information
    string const className = LvArray::system::demangleType< ObsoleteGroup< T > >(  );
    string const message = getObsoleteMessage();

    GEOS_LOG_RANK_0( GEOS_FMT(
                       "\n"
                       "********************************************************************************\n"
                       "* OBSOLETE GROUP WARNING\n"
                       "* Group: {}\n"
                       "* {}\n"
                       "* This group has been removed and will cause an error if used.\n"
                       "* (Group object created for input parsing compatibility only)\n"
                       "********************************************************************************",
                       className, message ) );
  }

  /**
   * @brief Log the obsolete error message and throw (called when actually trying to use)
   */
  [[noreturn]] void logObsoleteError() const
  {
    // Get the actual derived class type name from runtime type information
    string const className = LvArray::system::demangleType< ObsoleteGroup< T > >( );
    string const message = getObsoleteMessage();

    GEOS_ERROR( GEOS_FMT(
                  "\n"
                  "********************************************************************************\n"
                  "* OBSOLETE GROUP ERROR\n"
                  "* Group: {}\n"
                  "* {}\n"
                  "* This group has been removed and cannot be used.\n"
                  "* Please update your input deck to use the replacement.\n"
                  "********************************************************************************",
                  className, message ) );

    // This line should never be reached, but helps the compiler understand the noreturn attribute
    std::terminate();
  }

  /**
   * @brief Get the obsolete message
   *
   * This uses runtime polymorphism to call the derived class's obsoleteMessage()
   * method if it exists, otherwise returns a default message.
   *
   * @return The obsolete message string
   */
  virtual string getObsoleteMessage() const
  {
    return "This class is obsolete and has been removed.";
  }
};

} // namespace dataRepository
} // namespace geos

#endif // GEOS_DATAREPOSITORY_OBSOLETEGROUP_HPP_
