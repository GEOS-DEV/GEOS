/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */


/**
 * @file Group.hpp
 */

#ifndef GEOSX_DATAREPOSITORY_GROUP_HPP_
#define GEOSX_DATAREPOSITORY_GROUP_HPP_

#include "InputFlags.hpp"
#include "dataRepository/ObjectCatalog.hpp"
#include "MappedVector.hpp"
#include "RestartFlags.hpp"
#include "Wrapper.hpp"
#include "xmlWrapper.hpp"

#include <iostream>

#ifndef NOCHARTOSTRING_KEYLOOKUP
/// macro definition to enable/disable char * lookups
#define NOCHARTOSTRING_KEYLOOKUP 0
#endif

// Forward declaration of conduit::Node
namespace conduit
{
class Node;
}


/**
 * namespace to encapsulate GEOSX
 */
namespace geosx
{

/**
 * Encapsulates all dataRepository classes and functionality.
 */
namespace dataRepository
{

//START_SPHINX_INCLUDE_00
/// The default key type for entries in the hierarchy.
using keyType = string;

/// The default index type for entries the hierarchy.
using indexType = localIndex;

/**
 * @class Group
 *
 * The Group class serves as a "node" in a hierarchy of the dataRepository. The data structure is built as a
 * hierarchy of Group objects, or objects derived from group objects.
 */
class Group
{
public:
  //START_SPHINX_INCLUDE_01
  /// The template specialization of MappedVector to use for the collection of sub-Group objects.
  using subGroupMap = MappedVector< Group, Group *, keyType, indexType >;

  /// The template specialization of MappedVector to use for the collection wrappers objects.
  using wrapperMap = MappedVector< WrapperBase, WrapperBase *, keyType, indexType >;

  /**
   * @name Constructors/destructor
   */
  ///@{

  /**
   * @param[in] name the name of this object manager
   * @param[in] parent the parent Group
   */
  explicit Group( std::string const & name,
                  Group * const parent );


  /**
   * @brief Move constructor
   * @param[in] source source Group
   */
  Group( Group && source );

  /**
   * @brief Destructor, deletes all Groups and Wrappers owned by this Group
   */
  virtual ~Group();

  /**
   * @brief Deleted default constructor.
   */
  Group() = delete;

  /**
   * @brief Deleted copy constructor.
   */
  Group( Group const & ) = delete;

  /**
   * @brief Deleted move constructor.
   */
  Group( Group const && ) = delete;

  /**
   * @brief Deleted copy assignment operator.
   * @return
   */
  Group & operator=( Group const & ) = delete;

  /**
   * @brief Deleted move assignment operator.
   * @return
   */
  Group & operator=( Group && ) = delete;

  ///@}


  /**
   * @name Static Factory Catalog Functions
   */
  ///@{

  /**
   * @brief Type alias for catalog interface used by this class. See CatalogInterface.
   */
  using CatalogInterface = dataRepository::CatalogInterface< Group, std::string const &, Group * const >;

  /**
   * @brief Get the singleton catalog for this class.
   * @return reference to the catalog object
   */
  static CatalogInterface::CatalogType & GetCatalog();

  ///@}

  /**
   * @name Miscellaneous
   */
  ///@{

  /**
   * @brief Get typeid for current group
   * @return @p typeid(*this)
   */
  virtual const std::type_info & get_typeid() const
  {
    return typeid(*this);
  }

  /**
   * @brief Check a type_info against the type_info of this Group
   * @param typeToCheck value to check against
   * @return true of types are the same, false if not
   */
  bool CheckTypeID( std::type_info const & typeToCheck ) const
  {
    return typeToCheck == get_typeid();
  }

  /**
   * @brief Prints the data hierarchy recursively.
   * @param[in] indent The level of indentation to add to this level of output.
   */
  void PrintDataHierarchy( integer indent = 0 );

  /**
   * @brief Generates a table formatted string containing all input options.
   * @return a string containing a well formatted table containing input options.
   */
  string dumpInputOptions();

  ///@}

  //START_SPHINX_INCLUDE_REGISTER_GROUP
  /**
   * @name Sub-group registration interface
   */
  ///@{

  /**
   * @brief Register a new Group as a sub-group of current Group.
   *
   * @tparam T The type of the Group to add/register. This should be a type that derives from Group.
   * @param[in] name      The name of the group to use as a string key.
   * @param[in] newObject A unique_ptr to the object that is being registered.
   * @return              A pointer to the newly registered Group.
   *
   * Registers a Group or class derived from Group as a subgroup of this Group and takes ownership.
   */
  template< typename T = Group >
  T * RegisterGroup( std::string const & name, std::unique_ptr< Group > newObject );

  /**
   * @brief @copybrief RegisterGroup(std::string const &,std::unique_ptr<Group>)
   *
   * @tparam T The type of the Group to add/register. This should be a type that derives from Group.
   * @param[in] name          The name of the group to use as a string key.
   * @param[in] newObject     A unique_ptr to the object that is being registered.
   * @param[in] takeOwnership A flag to indicate whether or not the repository should
   *                          take ownership of the group.
   * @return                  A pointer to the newly registered Group.
   *
   * Registers a Group or class derived from Group as a subgroup of this Group and takes ownership
   * if @p takeOwnership is @p true.
   */
  template< typename T = Group >
  T * RegisterGroup( std::string const & name,
                     T * newObject,
                     bool const takeOwnership );

  /**
   * @brief @copybrief RegisterGroup(std::string const &,std::unique_ptr<Group>)
   *
   * @tparam T The type of the Group to add/register. This should be a type that derives from Group.
   * @param[in] name The name of the group to use as a string key.
   * @return         A pointer to the newly registered Group.
   *
   * Creates and registers a Group or class derived from Group as a subgroup of this Group.
   */
  template< typename T = Group >
  T * RegisterGroup( std::string const & name )
  {
    return RegisterGroup< T >( name, std::move( std::make_unique< T >( name, this )) );
  }

  /**
   * @brief @copybrief RegisterGroup(std::string const &,std::unique_ptr<Group>)
   *
   * @tparam T The type of the Group to add/register. This should be a type that derives from Group.
   * @param[in,out] keyIndex A KeyIndexT object that will be used to specify the name of
   *                         the new group. The index of the KeyIndex will also be set.
   * @return                 A pointer to the newly registered Group.
   *
   * Creates and registers a Group or class derived from Group as a subgroup of this Group.
   */
  template< typename T = Group >
  T * RegisterGroup( subGroupMap::KeyIndex & keyIndex )
  {
    T * rval = RegisterGroup< T >( keyIndex.Key(), std::move( std::make_unique< T >( keyIndex.Key(), this )) );
    keyIndex.setIndex( this->m_subGroups.getIndex( keyIndex.Key()) );
    return rval;
  }

  /**
   * @brief @copybrief RegisterGroup(std::string const &,std::unique_ptr<Group>)
   *
   * @tparam T The type of the Group to add/register. This should be a type that derives from Group.
   * @tparam TBASE The type whose type catalog will be used to look up the new sub-group type
   * @param[in] name        The name of the group to use as a string key.
   * @param[in] catalogName The catalog name of the new type.
   * @return                A pointer to the newly registered Group.
   *
   * Creates and registers a Group or class derived from Group as a subgroup of this Group.
   */
  template< typename T = Group, typename TBASE = Group >
  T * RegisterGroup( std::string const & name, std::string const & catalogName )
  {
    std::unique_ptr< TBASE > newGroup = TBASE::CatalogInterface::Factory( catalogName, name, this );
    return RegisterGroup< T >( name, std::move( newGroup ) );
  }

  /**
   * @brief Creates a new sub-Group using the ObjectCatalog functionality.
   * @param[in] childKey The name of the new object type's key in the
   *                     ObjectCatalog.
   * @param[in] childName The name of the new object in the collection of
   *                      sub-Groups.
   * @return A pointer to the new Group created by this function.
   */
  virtual Group * CreateChild( string const & childKey, string const & childName );

  ///@}

  //END_SPHINX_INCLUDE_REGISTER_GROUP

  /**
   * @name Group type-casting methods
   */
  ///@{

  /**
   * @brief       Downcast a Group *
   * @tparam T    Pointer type to downcast into
   * @param group A pointer the group to be casted
   * @return      A pointer to @p T that refers to the @p group
   */
  template< typename T >
  static T group_cast( Group * group )
  {
    return dynamicCast< T >( group );
  }

  /**
   * @brief       Downcast a Group const *
   * @tparam T    Pointer type to downcast into
   * @param group A pointer the group to be casted
   * @return      A pointer to @p T that refers to the @p group
   */
  template< typename T >
  static T group_cast( Group const * group )
  {
    return dynamicCast< T >( group );
  }

  /**
   * @brief       Downcast this Group
   * @tparam T    Pointer type to downcast into
   * @return      A pointer to \p T that refers to this object
   */
  template< typename T >
  T group_cast()
  {
    return dynamicCast< T >( this );
  }

  /**
   * @brief       Downcast this Group
   * @tparam T    Pointer type to downcast into
   * @return      A pointer to \p T that refers to this object
   */
  template< typename T >
  T group_cast() const
  {
    return dynamicCast< T >( this );
  }

  ///@}

  //START_SPHINX_INCLUDE_GET_GROUP
  /**
   * @name Sub-group retrieval methods.
   *
   * This collection of functions are used to get a sub-Group from the current group. Various methods
   * for performing the lookup are provided (localIndex, string, KeyIndex), and each have their
   * advantages and costs. The lowest cost lookup is the "localIndex" lookup. The KeyIndex lookup
   * will add a cost for checking to make sure the index stored in KeyIndex is valid (a string
   * compare, and a hash if it is incorrect). The string lookup is the full cost hash lookup every
   * time that it is called.
   *
   * The template parameter specifies the "type" that the caller expects to lookup, and thus attempts
   * to cast the pointer that is stored in m_subGroups to a pointer of the desired type. If this
   * cast fails, then a @p nullptr is returned. If no template parameter is specified then a default
   * type of Group is assumed.
   */
  ///@{

  /**
   * @brief Retrieve a sub-group from the current Group using an index.
   * @tparam T type of subgroup
   * @param[in] index the integral index of the group to retrieve.
   * @return A pointer to @p T that refers to the sub-group
   */
  template< typename T = Group >
  T * GetGroup( localIndex index )
  {
    return group_cast< T * >( m_subGroups[index] );
  }

  /**
   * @copydoc GetGroup(localIndex)
   */
  template< typename T = Group >
  T const * GetGroup( localIndex index ) const
  {
    return group_cast< T const * >( m_subGroups[index] );
  }

  /**
   * @brief Retrieve a sub-group from the current Group using a string.
   * @tparam T type of subgroup
   * @param[in] name the name/key of the group lookup and retrieve.
   * @return A pointer to @p T that refers to the sub-group
   */
  template< typename T = Group >
  T * GetGroup( string const & name )
  {
    return group_cast< T * >( m_subGroups[name] );
  }

  /**
   * @copydoc GetGroup(string const &)
   */
  template< typename T = Group >
  T const * GetGroup( string const & name ) const
  {
    return group_cast< T const * >( m_subGroups[name] );
  }

  /**
   * @brief Retrieve a sub-group from the current Group using a KeyIndexT.
   * @tparam T type of subgroup
   * @param[in,out] key the KeyIndex to use for the lookup
   * @return A pointer to @p T that refers to the sub-group
   */
  template< typename T = Group >
  T * GetGroup( subGroupMap::KeyIndex & key )
  {
    return group_cast< T * >( m_subGroups[key] );
  }

  /**
   * @copydoc GetGroup(subGroupMap::KeyIndex & key)
   */
  template< typename T = Group >
  T const * GetGroup( subGroupMap::KeyIndex & key ) const
  {
    return group_cast< T const * >( m_subGroups[key] );
  }

  /**
   * @copydoc GetGroup(subGroupMap::KeyIndex & key)
   * @note Const-correctness may be broken if the key is incorrect as
   *       @p key will be modified to contain the correct index.
   */
  template< typename T = Group >
  T * GetGroup( subGroupMap::KeyIndex const & key )
  {
    return group_cast< T * >( m_subGroups[key] );
  }

  /**
   * @copydoc GetGroup(subGroupMap::KeyIndex const & key)
   */
  template< typename T = Group >
  T const * GetGroup( subGroupMap::KeyIndex const & key ) const
  {
    return group_cast< T const * >( m_subGroups[key] );
  }

  /**
   * @brief Retrieve a group from the hierarchy using a path.
   * @tparam T type of subgroup
   * @param[in] path a unix-style string (absolute, relative paths valid)
   *                 to lookup the Group to return. Absolute paths search
   *                 from the tree root, while relative - from current group.
   * @return A pointer to @p T that refers to the sub-group
   */
  template< typename T = Group >
  T * GetGroupByPath( string const & path )
  {
    return const_cast< T * >(const_cast< Group const * >(this)->GetGroupByPath< T >( path ));
  }

  /**
   * @copydoc GetGroupByPath(string const &)
   */
  template< typename T = Group >
  T const * GetGroupByPath( string const & path ) const;

  //END_SPHINX_INCLUDE_GET_GROUP

  /**
   * @brief Get the subgroups object
   * @return a reference to the sub-group map.
   */
  subGroupMap & GetSubGroups()
  {
    return m_subGroups;
  }

  /**
   * @brief Get the subgroups object
   * @return a reference to const that points to the sub-group map.
   */
  subGroupMap const & GetSubGroups() const
  {
    return m_subGroups;
  }

  /**
   * @brief return the number of sub groups in this Group
   * @return number of sub groups in this Group
   */
  localIndex numSubGroups() const { return m_subGroups.size(); }

  /**
   * @brief Check whether a sub-group exists.
   * @param name the name of sub-group to search for
   * @return @p true if sub-group exists, @p false otherwise
   */
  bool hasGroup( std::string const & name ) const
  {
    return (m_subGroups[name] != nullptr);
  }

  ///@}

  /**
   * @name Functor application helpers
   *
   * This function is useful when trying to apply a functor that passes a pointer to an container,
   * but it is desired that the functor is only executed if the container can be casted to a certain
   * type. The variadic list consisting of CASTTYPE/S will be used recursively to check if the
   * container is able to be casted to the one of these types. The first type in the CASTTYPE/S list
   * will be used to execute the functor, and the function will return true.
   */
  ///@{

  /** @cond DO_NOT_DOCUMENT */
  template< typename CONTAINERTYPE, typename LAMBDA >
  static bool applyLambdaToContainer( CONTAINERTYPE const * const GEOSX_UNUSED_ARG( group ), LAMBDA && GEOSX_UNUSED_ARG( lambda ) )
  { return false; }

  template< typename CONTAINERTYPE, typename LAMBDA >
  static bool applyLambdaToContainer( CONTAINERTYPE * const GEOSX_UNUSED_ARG( group ), LAMBDA && GEOSX_UNUSED_ARG( lambda ) )
  { return false; }
  /** @endcond */

  /**
   * @brief Apply a given functor to a container if the container can be
   *        cast to one of the specified types.
   * @tparam CONTAINERTYPE the type of container
   * @tparam CASTTYPE      the first type that will be used in the attempted casting of container
   * @tparam CASTTYPES     a variadic list of types that will be used in the attempted casting of container
   * @tparam LAMBDA        the type of lambda function to call in the function
   * @param[in] container  a pointer to the container which will be passed to the lambda function
   * @param[in] lambda     the lambda function to call in the function
   * @return               a boolean to indicate whether the lambda was successfully applied to the container.
   */
  template< typename CONTAINERTYPE, typename CASTTYPE, typename ... CASTTYPES, typename LAMBDA >
  static bool applyLambdaToContainer( CONTAINERTYPE const * const container, LAMBDA && lambda )
  {
    bool rval = false;
    CASTTYPE const * const castedContainer = dynamic_cast< CASTTYPE const * >( container );
    if( castedContainer!= nullptr )
    {
      lambda( castedContainer );
      rval = true;
    }
    else
    {
      rval = applyLambdaToContainer< CONTAINERTYPE, CASTTYPES... >( container, std::forward< LAMBDA >( lambda ) );
    }
    return rval;
  }

  /**
   * @copydoc applyLambdaToContainer(CONTAINERTYPE const * const, LAMBDA &&)
   */
  template< typename CONTAINERTYPE, typename CASTTYPE, typename ... CASTTYPES, typename LAMBDA >
  static bool applyLambdaToContainer( CONTAINERTYPE * const container, LAMBDA && lambda )
  {
    bool rval = false;
    CASTTYPE * const castedContainer = dynamic_cast< CASTTYPE * >( container );
    if( castedContainer!= nullptr )
    {
      lambda( castedContainer );
      rval = true;
    }
    else
    {
      rval = applyLambdaToContainer< CONTAINERTYPE, CASTTYPES... >( container, std::forward< LAMBDA >( lambda ) );
    }
    return rval;
  }
  ///@}


  //START_SPHINX_INCLUDE_LOOP_INTERFACE
  /**
   * @name Functor-based subgroup iteration
   *
   * These functions loop over sub-groups and executes a functor that uses the sub-group as an
   * argument. The functor is only executed if the group can be cast to a certain type specified
   * by the @p ROUPTYPE/S pack. The variadic list consisting of @p GROUPTYPE/S will be used recursively
   * to check if the group is able to be cast to the one of these types. The first type in the
   * @p GROUPTYPE/S list will be used to execute the functor, and the next sub-group will be processed.
   */
  ///@{

  /**
   * @brief Apply the given functor to subgroups that can be casted to one of specified types.
   * @tparam GROUPTYPE  the first type that will be used in the attempted casting of group.
   * @tparam GROUPTYPES a variadic list of types that will be used in the attempted casting of group.
   * @tparam LAMBDA     the type of functor to call
   * @param[in] lambda  the functor to call on subgroups
   */
  template< typename GROUPTYPE = Group, typename ... GROUPTYPES, typename LAMBDA >
  void forSubGroups( LAMBDA lambda )
  {
    for( auto & subGroupIter : m_subGroups )
    {
      applyLambdaToContainer< Group, GROUPTYPE, GROUPTYPES... >( subGroupIter.second, [&]( auto * const castedSubGroup )
          {
            lambda( castedSubGroup );
          } );
    }
  }

  /**
   * @copydoc forSubGroups(LAMBDA)
   */
  template< typename GROUPTYPE = Group, typename ... GROUPTYPES, typename LAMBDA >
  void forSubGroups( LAMBDA lambda ) const
  {
    for( auto const & subGroupIter : m_subGroups )
    {
      applyLambdaToContainer< Group, GROUPTYPE, GROUPTYPES... >( subGroupIter.second, [&]( auto const * const castedSubGroup )
          {
            lambda( castedSubGroup );
          } );
    }
  }

  /**
   * @copybrief forSubGroups(LAMBDA)
   * @tparam GROUPTYPE  the first type that will be used in the attempted casting of group.
   * @tparam GROUPTYPES a variadic list of types that will be used in the attempted casting of group.
   * @tparam LAMBDA     the type of functor to call
   * @param[in] lambda        the functor to call
   * @param[in] subgroupNames list of subgroup names to apply the functor to
   */
  template< typename GROUPTYPE = Group, typename ... GROUPTYPES, typename LAMBDA >
  void forSubGroups( string_array const & subgroupNames, LAMBDA lambda )
  {
    for( string const & subgroupName : subgroupNames )
    {
      applyLambdaToContainer< Group, GROUPTYPE, GROUPTYPES... >( GetGroup( subgroupName ), [&]( auto * const castedSubGroup )
          {
            lambda( castedSubGroup );
          } );
    }
  }

  /**
   * @copydoc forSubGroups(string_array const &, LAMBDA)
   */
  template< typename GROUPTYPE = Group, typename ... GROUPTYPES, typename LAMBDA >
  void forSubGroups( string_array const & subgroupNames, LAMBDA lambda ) const
  {
    for( string const & subgroupName : subgroupNames )
    {
      applyLambdaToContainer< Group, GROUPTYPE, GROUPTYPES... >( GetGroup( subgroupName ), [&]( auto const * const castedSubGroup )
          {
            lambda( castedSubGroup );
          } );
    }
  }
  ///@}

  /**
   * @name Functor-based wrapper iteration
   *
   * These functions loop over the wrappers contained in this group, and executes a functor that
   * uses the Wrapper as an argument. The functor is only executed if the Wrapper can be casted to
   * a certain type specified by the @p TYPE/S pack. The variadic list consisting of
   * @p TYPE/S will be used recursively to check if the Wrapper is able to be casted to the
   * one of these types. The first type in the @p WRAPPERTYPE/S list will be used to execute the
   * functor, and the next Wrapper will be processed.
   */
  ///@{

  /**
   * @brief Apply the given functor to wrappers.
   * @tparam LAMBDA the type of functor to call
   * @param[in] lambda  the functor to call
   */
  template< typename LAMBDA >
  void forWrappers( LAMBDA lambda )
  {
    for( auto & wrapperIter : m_wrappers )
    {
      lambda( *wrapperIter.second );
    }
  }

  /**
   * @copydoc forWrappers(LAMBDA)
   */
  template< typename LAMBDA >
  void forWrappers( LAMBDA lambda ) const
  {
    for( auto const & wrapperIter : m_wrappers )
    {
      lambda( *wrapperIter.second );
    }
  }

  /**
   * @brief Apply the given functor to wrappers that can be cast to one of specified types.
   * @tparam TYPE   the first type that will be used in the attempted casting of Wrapper
   * @tparam TYPES  a variadic list of types that will be used in the attempted casting of Wrapper
   * @tparam LAMBDA the type of functor to call
   * @param[in] lambda  the functor to call
   */
  template< typename TYPE, typename ... TYPES, typename LAMBDA >
  void forWrappers( LAMBDA lambda )
  {
    for( auto & wrapperIter : m_wrappers )
    {
      applyLambdaToContainer< WrapperBase, Wrapper< TYPE >, Wrapper< TYPES >... >( wrapperIter.second,
                                                                                   std::forward< LAMBDA >( lambda ));
    }
  }

  /**
   * @brief Apply the given functor to wrappers that can be cast to one of specified types.
   * @tparam TYPE   the first type that will be used in the attempted casting of Wrapper
   * @tparam TYPES  a variadic list of types that will be used in the attempted casting of Wrapper
   * @tparam LAMBDA the type of functor to call
   * @param[in] lambda  the functor to call
   */
  template< typename TYPE, typename ... TYPES, typename LAMBDA >
  void forWrappers( LAMBDA lambda ) const
  {
    for( auto const & wrapperIter : m_wrappers )
    {
      applyLambdaToContainer< WrapperBase, Wrapper< TYPE >, Wrapper< TYPES >... >( wrapperIter.second,
                                                                                   std::forward< LAMBDA >( lambda ));
    }
  }

  ///@}
  //END_SPHINX_INCLUDE_LOOP_INTERFACE

  /**
   * @name Initialization and data registration
   */
  ///@{

  /**
   * @brief Run initialization functions on this and all subgroups.
   * @param[in] group A group that is passed in to the initialization functions
   *                  in order to facilitate the initialization.
   *
   * This function will first call InitializePreSubGroups() on this Group, then
   * loop over all subgroups and call Initialize() on them, then
   * call InitializePostSubGroups() on this Group.
   *
   * @note The order in which the sub-Groups are iterated over is defined by
   * InitializationOrder().
   */
  void Initialize( Group * const group );

  /**
   * @brief Sets the initialization order for sub-Groups.
   * @param[out] order An array of strings that define the iteration order.
   *
   * This function will fill the @p order array that is used to specify the
   * order in which the Initialize() function loops over sub-Groups. If a
   * custom order is required by a derived type, this function should be
   * overridden with a implementation that specifies the desired order.
   */
  virtual void InitializationOrder( string_array & order );


  /**
   * @brief Initialization routine to be called after calling
   *        ApplyInitialConditions().
   * @param[in] group A group that is passed in to the initialization functions
   *                  in order to facilitate the initialization.
   *
   * This function provides a capability for post-initial condition problem
   * initialization. First the InitializePostInitialConditions_PreSubGroups()
   * function is called on this Group. Then there is a loop over all subgroups
   * and InitializePostInitialConditions() is called on them. Finally, the
   * InitializePostInitialConditions_PostSubGroups() function is called this
   * Group.
   *
   * @note The order in which the sub-Groups are iterated over is defined by
   * InitializationOrder().
   */
  void InitializePostInitialConditions( Group * const group );

  /**
   * @brief Initialization routine to be called after calling reading a restart file.
   * @param domain the physical domain
   *
   * This functions recurses and calls postRestartInitialization() on nested sub-groups
   * at any depth, providing a capability to add custom post-restart initialization.
   */
  void postRestartInitializationRecursive( Group * const domain );

  /**
   * @brief Recursively read values using ProcessInputFile() from the input
   *        file and put them into the wrapped values for this group.
   * @param[in] targetNode the XML node that to extract input values from.
   */
  void ProcessInputFileRecursive( xmlWrapper::xmlNode & targetNode );

  /**
   * @brief Recursively call PostProcessInput() to apply post processing after
   * reading input values.
   */
  void PostProcessInputRecursive();

  ///@}

  //START_SPHINX_INCLUDE_REGISTER_WRAPPER
  /**
   * @name Wrapper registration interface
   */
  ///@{

  /**
   * @brief Create and register a Wrapper around a new object.
   * @tparam     T the type of the wrapped object
   * @tparam TBASE the base type to cast the returned wrapper to
   * @param[in]  name the name of the wrapper to use as a string key
   * @param[out] rkey a pointer to a index type that will be filled with the new
   *             Wrapper index in this Group
   * @return     a pointer to the newly registered/created Wrapper
   */
  template< typename T, typename TBASE=T >
  Wrapper< TBASE > * registerWrapper( std::string const & name,
                                      wrapperMap::KeyIndex::index_type * const rkey = nullptr );

  /**
   * @copybrief registerWrapper(std::string const &,wrapperMap::KeyIndex::index_type * const)
   * @tparam     T the type of the wrapped object
   * @tparam TBASE the base type to cast the returned wrapper to
   * @param[in] viewKey The KeyIndex that contains the name of the new Wrapper.
   * @return            a pointer to the newly registered/created Wrapper
   */
  template< typename T, typename TBASE=T >
  Wrapper< TBASE > * registerWrapper( Group::wrapperMap::KeyIndex & viewKey );

  /**
   * @brief Register a Wrapper around a given object and take ownership.
   * @tparam T the type of the wrapped object
   * @param[in] name      the name of the wrapper to use as a string key
   * @param[in] newObject an owning pointer to the object that is being registered
   * @return              a pointer to the newly registered/created Wrapper
   */
  template< typename T >
  Wrapper< T > * registerWrapper( std::string const & name,
                                  std::unique_ptr< T > newObject );

  /**
   * @brief Register a Wrapper around a given object and conditionally take ownership.
   * @tparam T the type of the wrapped object
   * @param[in] name          the name of the wrapper to use as a string key
   * @param[in] newObject     a pointer to the object that is being registered
   * @param[in] takeOwnership a flag to indicate whether or not the repository should
   *                          take ownership of the group
   * @return                  a pointer to the newly registered/created Wrapper
   */
  template< typename T >
  Wrapper< T > * registerWrapper( std::string const & name,
                                  T * newObject,
                                  bool takeOwnership );

  /**
   * @brief Register and take ownership of an existing Wrapper.
   * @param name    the name of the wrapper to use as a string key
   * @param wrapper a pointer to the an existing wrapper.
   * @return        an un-typed pointer to the newly registered/created wrapper
   */
  WrapperBase * registerWrapper( string const & name,
                                 WrapperBase * const wrapper );

  /**
   * @brief Removes a Wrapper from this group.
   * @param name the name of the Wrapper to remove from this group.
   */
  void deregisterWrapper( string const & name );

  ///@}
  //END_SPHINX_INCLUDE_REGISTER_WRAPPER

  /**
   * @name Schema generation methods
   */
  ///@{

  /**
   * @brief Build a complete datastructure for schema generation.
   * @param level indent level for printing out the structure
   */
  void GenerateDataStructureSkeleton( integer const level )
  {
    ExpandObjectCatalogs();
    std::string indent( level*2, ' ' );

    for( auto const & subGroupIter : m_subGroups )
    {
      std::cout << indent << subGroupIter.second->getName() << std::endl;
      subGroupIter.second->GenerateDataStructureSkeleton( level + 1 );
    }
  }

  /**
   * @brief Expand any catalogs in the data structure.
   */
  virtual void ExpandObjectCatalogs() {}

  /**
   * @brief Inform the schema generator of any deviations between the xml and GEOS data structures.
   * @param schemaRoot        XML node corresponding to the root
   * @param schemaParent      XML node for the parent node
   * @param documentationType type of XML schema generated
   */
  virtual void SetSchemaDeviations( xmlWrapper::xmlNode schemaRoot,
                                    xmlWrapper::xmlNode schemaParent,
                                    integer documentationType );

  ///@}

  /**
   * @name Mesh data registration
   */
  ///@{

  /**
   * @brief Calls RegisterDataOnMesh() recursively.
   * @param[in,out] MeshBodies the group of MeshBody objects to register data on.
   */
  virtual void RegisterDataOnMeshRecursive( Group * const MeshBodies );

  /**
   * @brief Register data on mesh entities.
   * @param[in,out] MeshBody the group of MeshBody objects to register data on.
   *
   * This function is used to register data on mesh entities such as the NodeManager,
   * FaceManager...etc.
   */
  virtual void RegisterDataOnMesh( Group * const MeshBody );

  ///@}

  /**
   * @name Packing/upacking methods
   */
  ///@{

  /**
   * @brief Get the size required to pack a list of wrappers.
   * @param[in] wrapperNames an array that contains the names of the wrappers to pack.
   * @param[in] recursive    whether or not to perform a recursive pack.
   * @param[in] on_device    whether to use device-based packing functions
   *                         (buffer must be either pinned or a device pointer)
   * @return                 the size of the buffer required to pack the wrappers.
   */
  virtual localIndex PackSize( string_array const & wrapperNames,
                               integer const recursive,
                               bool on_device = false ) const;

  /**
   * @brief Get the size required to pack a list of indices within a list of wrappers.
   * @param[in] wrapperNames an array that contains the names of the wrappers to pack.
   * @param[in] packList     the list of indices to pack
   * @param[in] recursive    whether or not to perform a recursive pack.
   * @param[in] on_device    whether to use device-based packing functions
   *                         (buffer must be either pinned or a device pointer)
   * @return                 the size of the buffer required to pack the wrapper indices.
   */
  virtual localIndex PackSize( string_array const & wrapperNames,
                               arrayView1d< localIndex const > const & packList,
                               integer const recursive,
                               bool on_device = false ) const;

  /**
   * @brief Pack a list of wrappers to a buffer.
   * @param[in,out] buffer   the buffer that will be packed.
   * @param[in] wrapperNames an array that contains the names of the wrappers to pack.
   * @param[in] recursive    whether or not to perform a recursive pack.
   * @param[in] on_device    whether to use device-based packing functions
   *                         (buffer must be either pinned or a device pointer)
   * @return                 the size of data packed to the buffer.
   *
   * This function takes in a reference to a pointer @p buffer, and packs data specified by
   * @p wrapperNames, and @p recursive to that pointer location. The
   * pointer is altered and returned to the new location corresponding the
   * original value of @p buffer plus the size of data packed to the buffer.
   *
   */
  virtual localIndex Pack( buffer_unit_type * & buffer,
                           string_array const & wrapperNames,
                           integer const recursive,
                           bool on_device = false ) const;

  /**
   * @brief Pack a list of indices within a list of wrappers.
   * @param[in,out] buffer   the buffer that will be packed.
   * @param[in] wrapperNames an array that contains the names of the wrappers to pack.
   * @param[in] packList     the list of indices to pack
   * @param[in] recursive    whether or not to perform a recursive pack.
   * @param[in] on_device    whether to use device-based packing functions
   *                         (buffer must be either pinned or a device pointer)
   * @return                 the size of data packed to the buffer.
   *
   * This function takes in a reference to a pointer @p buffer, and packs data specified by
   * @p wrapperNames, @p packList, and @p recursive to that pointer location. The
   * pointer is altered and returned to the new location corresponding the
   * original value of @p buffer plus the size of data packed to the buffer.
   */
  virtual localIndex Pack( buffer_unit_type * & buffer,
                           string_array const & wrapperNames,
                           arrayView1d< localIndex const > const & packList,
                           integer const recursive,
                           bool on_device = false ) const;

  /**
   * @brief Unpack a buffer.
   * @param[in,out] buffer   the buffer to unpack
   * @param[in,out] packList the list of indices that will be unpacked.
   * @param[in] recursive    whether or not to perform a recursive unpack.
   * @param[in] on_device    whether to use device-based packing functions
   *                         (buffer must be either pinned or a device pointer)
   * @return                 the number of bytes unpacked.
   *
   * This function takes a reference to a pointer to const buffer type, and
   * unpacks data from that buffer into the current Group. If @p packList
   * is non-empty, then a check is made to ensure that the data that is
   * unpacked matches @p packList. If @p packList is empty, the values
   * of the indices that are unpacked are stored and returned in packList.
   */
  virtual localIndex Unpack( buffer_unit_type const * & buffer,
                             arrayView1d< localIndex > & packList,
                             integer const recursive,
                             bool on_device = false );

  ///@}

  //***********************************************************************************************

  //START_SPHINX_INCLUDE_GET_WRAPPER
  /**
   * @name Untyped wrapper retrieval methods
   *
   * These functions query the collection of Wrapper objects for the given
   * index/name/KeyIndex and returns a WrapperBase pointer to the object if
   * it exists. If it is not found, nullptr is returned.
   */
  ///@{

  /**
   * @brief Retrieve a WrapperBase stored in this group.
   * @param[in] index an integral lookup value used to search the collection of wrappers.
   * @return          a pointer to the WrapperBase that resulted from the lookup.
   */
  WrapperBase const * getWrapperBase( indexType const index ) const
  { return m_wrappers[index]; }

  /**
   * @copydoc getWrapperBase(indexType const) const
   */
  WrapperBase * getWrapperBase( indexType const index )
  { return m_wrappers[index]; }

  /**
   * @brief Retrieve a WrapperBase stored in this group.
   * @param[in] name a string lookup value used to search the collection of wrappers.
   * @return         a pointer to the WrapperBase that resulted from the lookup.
   */
  WrapperBase const * getWrapperBase( std::string const & name ) const
  { return m_wrappers[name]; }

  /**
   * @copydoc getWrapperBase(std::string const &) const
   */
  WrapperBase * getWrapperBase( std::string const & name )
  { return m_wrappers[name]; }

  /**
   * @brief Retrieve a WrapperBase stored in this group.
   * @param[in] keyIndex a KeyIndex lookup value used to search the collection of wrappers.
   * @return             a pointer to the WrapperBase that resulted from the lookup.
   */
  WrapperBase const * getWrapperBase( wrapperMap::KeyIndex const & keyIndex ) const
  { return m_wrappers[keyIndex]; }

  /**
   * @copydoc getWrapperBase(wrapperMap::KeyIndex const &) const
   */
  WrapperBase * getWrapperBase( wrapperMap::KeyIndex const & keyIndex )
  { return m_wrappers[keyIndex]; }

  /**
   * @brief
   * @param name
   * @return
   */
  indexType getWrapperIndex( std::string const & name ) const
  {
    return m_wrappers.getIndex( name );
  }

  /**
   * @brief Get access to the internal wrapper storage.
   * @return a reference to wrapper map
   */
  wrapperMap const & wrappers() const
  {
    return m_wrappers;
  }

  /**
   * @copydoc wrappers() const
   */
  wrapperMap & wrappers()
  {
    return m_wrappers;
  }

  ///@}

  /**
   * @name Typed wrapper retrieval methods
   *
   * These functions query the collection of Wrapper objects for the given
   * index/key and returns a Wrapper<T> pointer to the object if
   * it exists. The template parameter @p T is used to perform a cast
   * on the WrapperBase pointer that is returned by the lookup, into
   * a Wrapper<T> pointer. If the wrapper is not found, or the
   * WrapperBase pointer cannot be cast to a Wrapper<T> pointer, then nullptr
   * is returned.
   */
  ///@{

  /**
   * @brief Check if a wrapper exists
   * @tparam LOOKUP_TYPE the type of key used to perform the lookup.
   * @param[in] lookup a lookup value used to search the collection of wrappers
   * @return @p true if wrapper exists (regardless of type), @p false otherwise
   */
  template< typename LOOKUP_TYPE >
  bool hasWrapper( LOOKUP_TYPE const & lookup ) const
  {
    return (m_wrappers[lookup] != nullptr);
  }

  /**
   * @brief Retrieve a Wrapper stored in this group.
   * @tparam T           the object type contained in the Wrapper
   * @tparam LOOKUP_TYPE the type of key used to perform the lookup
   * @param[in] index    a lookup value used to search the collection of wrappers
   * @return             A pointer to the Wrapper<T> that resulted from the lookup
   */
  template< typename T, typename LOOKUP_TYPE >
  Wrapper< T > const * getWrapper( LOOKUP_TYPE const & index ) const
  {
    return dynamicCast< Wrapper< T > const * >( m_wrappers[index] );
  }

  /**
   * @copydoc getWrapper(LOOKUP_TYPE const &) const
   */
  template< typename T, typename LOOKUP_TYPE >
  Wrapper< T > * getWrapper( LOOKUP_TYPE const & index )
  { return const_cast< Wrapper< T > * >( const_cast< Group const * >(this)->getWrapper< T >( index ) ); }

  /**
   * @brief Retrieve a Wrapper stored in this group.
   * @tparam T      the object type contained in the Wrapper
   * @param[in] key a string lookup value used to search the collection of wrappers
   * @return        a pointer to the Wrapper<T> that resulted from the lookup
   */
  template< typename T >
  Wrapper< T > const * getWrapper( char const * const key ) const
  { return getWrapper< T >( string( key ) ); }

  /**
   * @copydoc getWrapper(char const * const) const
   */
  template< typename T >
  Wrapper< T > * getWrapper( char const * const key )
  { return getWrapper< T >( string( key ) ); }

  ///@}

  /**
   * @name Wrapper data access methods.
   *
   * These functions can be used to get referece/pointer access to the data
   * stored by wrappers in this group. They are essentially just shortcuts for
   * @p Group::getWrapper() and @p Wrapper<T>::getReference()/getPointer().
   * An additional template parameter can be provided to cast the return pointer
   * or reference to a base class pointer or reference (e.g. Array to ArrayView).
   */
  ///@{

  /**
   * @brief Look up a wrapper and get reference to wrapped object.
   * @tparam T           return value type
   * @tparam WRAPPEDTYPE wrapped value type (by default, same as return)
   * @tparam LOOKUP_TYPE type of value used for wrapper lookup
   * @param lookup       value for wrapper lookup
   * @return             reference to @p T
   *
   * @note An error will be raised if wrapper does not exist or type cast is invalid.
   */
  template< typename T, typename WRAPPEDTYPE=T, typename LOOKUP_TYPE >
  typename std::enable_if< std::is_same< T, WRAPPEDTYPE >::value, T const & >::type
  getReference( LOOKUP_TYPE const & lookup ) const
  {
    Wrapper< WRAPPEDTYPE > const * wrapper = getWrapper< WRAPPEDTYPE >( lookup );
    if( wrapper == nullptr )
    {
      if( hasWrapper( lookup ) )
      {
        GEOSX_ERROR( "call to getWrapper results in nullptr but a view exists. Most likely given the incorrect type. lookup : " << lookup );
      }
      GEOSX_ERROR( "call to getWrapper results in nullptr and a view does not exist. lookup : " << lookup );
    }

    return wrapper->reference();
  }

  /**
   * @copydoc getReference(LOOKUP_TYPE const &) const
   */
  template< typename T, typename WRAPPEDTYPE=T, typename LOOKUP_TYPE >
  typename std::enable_if< !std::is_same< T, WRAPPEDTYPE >::value, T const & >::type
  getReference( LOOKUP_TYPE const & lookup ) const
  {
    static_assert( std::is_base_of< WRAPPEDTYPE, T >::value, "incorrect template arguments" );
    Wrapper< WRAPPEDTYPE > const * wrapper = getWrapper< WRAPPEDTYPE >( lookup );
    if( wrapper == nullptr )
    {
      if( hasWrapper( lookup ) )
      {
        GEOSX_ERROR( "call to getWrapper results in nullptr but a view exists. Most likely given the incorrect type. lookup : " << lookup );
      }
      GEOSX_ERROR( "call to getWrapper results in nullptr and a view does not exist. lookup : " << lookup );
    }

    return dynamicCast< T const & >( wrapper->reference() );
  }

  /**
   * @copydoc getReference(LOOKUP_TYPE const &) const
   */
  template< typename T, typename WRAPPEDTYPE=T, typename LOOKUP_TYPE >
  T & getReference( LOOKUP_TYPE const & lookup )
  { return const_cast< T & >( const_cast< const Group * >(this)->template getReference< T, WRAPPEDTYPE, LOOKUP_TYPE >( lookup ) ); }

  /**
   * @copybrief getReference(LOOKUP_TYPE const &) const
   * @tparam T           return value type
   * @tparam WRAPPEDTYPE wrapped value type (by default, same as return)
   * @param name         name of the wrapper
   * @return             reference to @p T
   *
   * @note An error will be raised if wrapper does not exist or type cast is invalid.
   */
  template< typename T, typename WRAPPEDTYPE=T >
  T const & getReference( char const * const name ) const
  { return getReference< T, WRAPPEDTYPE >( string( name ) ); }

  /**
   * @copydoc getReference(char const * const) const
   */
  template< typename T, typename WRAPPEDTYPE=T >
  T & getReference( char const * const name )
  { return const_cast< T & >( const_cast< const Group * >(this)->getReference< T, WRAPPEDTYPE >( name ) ); }

  /**
   * @brief Look up a wrapper and get reference to wrapped object.
   * @tparam T           return value type
   * @tparam LOOKUP_TYPE type of value used for wrapper lookup
   * @param lookup       value for wrapper lookup
   * @return             pointer to @p T
   *
   * @note @p nullptr will be returned if wrapper does not exist or type cast is invalid.
   */
  template< typename T, typename LOOKUP_TYPE >
  T const * getPointer( LOOKUP_TYPE const & lookup ) const
  {
    T const * rval = nullptr;
    Wrapper< T > const * wrapper = getWrapper< T >( lookup );
    if( wrapper != nullptr )
    {
      rval = wrapper->getPointer();
    }
    return rval;
  }

  /**
   * @copydoc getPointer(LOOKUP_TYPE const &) const
   */
  template< typename T, typename LOOKUP_TYPE >
  T * getPointer( LOOKUP_TYPE const & lookup )
  { return const_cast< T * >( const_cast< Group const * >(this)->getPointer< T >( lookup )); }

  /**
   * @copybrief getPointer(LOOKUP_TYPE const &) const
   * @tparam T           return value type
   * @param name         name of the wrapper
   * @return             pointer to @p T
   *
   * @note nullptr will be returned if wrapper does not exist or type cast is invalid.
   */
  template< typename T >
  T const * getPointer( char const * const name ) const
  { return getPointer< T >( string( name ) ); }

  /**
   * @copydoc getPointer(char const * const) const
   */
  template< typename T >
  T * getPointer( char const * const name )
  { return getPointer< T >( string( name ) ); }
  //END_SPHINX_INCLUDE_GET_WRAPPER

  ///@}

  /**
   * @name Size/capacity management
   */
  /// @{

  /**
   * @brief Resize the group and all contained wrappers that resize with parent.
   * @param newsize the new size of the group
   */
  virtual void resize( localIndex const newsize );

  /**
   * @brief Set the new capacity and reserve it in all wrappers that resize with parent.
   * @param newsize the new capacity of the group
   */
  virtual void reserve( indexType const newsize );

  /**
   * @brief Get the "capacity" of the group, which determines the capacity of resizable wrappers.
   * @return capacity of this group
   */
  inline localIndex capacity() const
  {
    return m_capacity;
  }

  /**
   * @brief Get the "size" of the group, which determines the number of elements in resizable wrappers.
   * @return size of this group
   */
  inline localIndex size() const
  {
    return m_size;
  }

  /// @}

  /**
   * @name Basic group properties
   */
  ///@{

  /**
   * @brief Get group name.
   * @return group name
   */
  inline const string getName() const
  {
    return m_name;
  }

  /**
   * @brief Access the group's parent.
   * @return pointer to parent Group
   */
  Group * getParent()             { return m_parent; }

  /**
   * @copydoc getParent()
   */
  Group const * getParent() const { return m_parent; }

  /**
   * @brief Get the group's index withing its parent group
   * @return integral index of current group within its parent
   */
  localIndex getIndexInParent() const
  {
    return m_parent->GetSubGroups().getIndex( this->m_name );
  }

  /**
   * @brief Get flags that control restart output of this group.
   * @return the current value of restart flags
   */
  RestartFlags getRestartFlags() const { return m_restart_flags; }

  /**
   * @brief Set flags that control restart output of this group.
   * @param flags the new value of restart flags
   */
  void setRestartFlags( RestartFlags flags ) { m_restart_flags = flags; }

  /**
   * @brief Get input flags for schema generation.
   * @return the current value of input flags
   */
  InputFlags getInputFlags() const { return m_input_flags; }

  /**
   * @brief Set input flags for schema generation.
   * @param flags the new value of input flags
   */
  void setInputFlags( InputFlags flags ) { m_input_flags = flags; }

  ///@}

  /**
   * @name Restart output methods
   */
  ///@{

  /**
   * @brief Get the Conduit node object associated with this group
   * @return reference to inner conduit::Node member
   */
  conduit::Node & getConduitNode()
  {
    return m_conduitNode;
  }

  /**
   * @brief Register the group and its wrappers with Conduit.
   */
  void prepareToWrite();

  /**
   * @brief Write the group and its wrappers into Conduit.
   */
  void finishWriting();

  /**
   * @brief Read the group and its wrappers from Conduit.
   */
  void loadFromConduit();

  /// Enable verbosity input for object
  void enableLogLevelInput();

  /// @return The verbosity level
  integer getLogLevel() const { return m_logLevel; }
  ///@}

protected:

  /**
   * @name Overridable initialization hooks
   *
   * These methods can be overridden by derived classes to customize
   * input post processing and object initialization.
   */
  ///@{

  /**
   * This function provides capability to post process input values prior to
   * any other initialization operations.
   */
  virtual void PostProcessInput() {}

  /**
   * @brief Called by Initialize() prior to initializing sub-Groups.
   * @param[in] group A group that is passed in to the initialization functions
   *                  in order to facilitate the initialization.
   */
  virtual void InitializePreSubGroups( Group * const group )
  {
    GEOSX_UNUSED_VAR( group );
  }

  /**
   * @brief Called by Initialize() after to initializing sub-Groups.
   * @param[in] group A group that is passed in to the initialization functions
   *                  in order to facilitate the initialization.
   */
  virtual void InitializePostSubGroups( Group * const group )
  {
    GEOSX_UNUSED_VAR( group );
  }

  /**
   * @brief Called by InitializePostInitialConditions() prior to initializing sub-Groups.
   * @param[in] group A group that is passed in to the initialization functions
   *                  in order to facilitate the initialization.
   */
  virtual void InitializePostInitialConditions_PreSubGroups( Group * const group )
  {
    GEOSX_UNUSED_VAR( group );
  }

  /**
   * @brief Called by InitializePostInitialConditions() after to initializing sub-Groups.
   * @param[in] group A group that is passed in to the initialization functions
   *                  in order to facilitate the initialization.
   */
  virtual void InitializePostInitialConditions_PostSubGroups( Group * const group )
  {
    GEOSX_UNUSED_VAR( group );
  }

  /**
   * @brief Performs initialization required after reading from a restart file.
   * @param domain A pointer to the domain partition.
   */
  virtual void postRestartInitialization( Group * const domain )
  {
    GEOSX_UNUSED_VAR( domain );
  }

  ///@}

private:
  /**
   * @brief Read values from the input file and put them into the
   *        wrapped values for this group.
   * @param[in] targetNode the XML node that to extract input values from.
   */
  virtual void ProcessInputFile( xmlWrapper::xmlNode const & targetNode );

  //START_SPHINX_INCLUDE_02
  /// The parent Group that contains "this" Group in its "sub-Group" collection.
  Group * m_parent = nullptr;

  /// The container for the collection of all wrappers continued in "this" Group.
  wrapperMap m_wrappers;

  /// The container for the collection of all sub-groups contained in "this" Group.
  subGroupMap m_subGroups;

  /// The size/length of this Group...and all Wrapper<> that are are specified to have the same size as their
  /// owning group.
  indexType m_size;

  /// The capacity for wrappers in this group...and all Wrapper<> that are specified to have the same size as their
  /// owning group.
  indexType m_capacity;

  /// The name/key of this Group in its parent collection of sub-Groups.
  string m_name;

  /// Verbosity flag for group logs
  integer m_logLevel;

  RestartFlags m_restart_flags; ///< Restart flag for this group...and
                                ///< subsequently all wrappers in this group
  InputFlags m_input_flags;     ///< Input flag for this group

  /// Reference to the conduit::Node that mirrors this group
  conduit::Node & m_conduitNode;

};

/**
 * @brief Type alias for KeyIndexT type used for sub-group lookups.
 */
using GroupKey = Group::subGroupMap::KeyIndex;

/**
 * @brief Type alias for KeyIndexT type used for wrapper lookups.
 */
using ViewKey = Group::wrapperMap::KeyIndex;



template< typename T >
T * Group::RegisterGroup( std::string const & name,
                          std::unique_ptr< Group > newObject )
{
  return dynamicCast< T * >( m_subGroups.insert( name, newObject.release(), true ) );
}


template< typename T >
T * Group::RegisterGroup( std::string const & name,
                          T * newObject,
                          bool const takeOwnership )
{
  return dynamicCast< T * >( m_subGroups.insert( name, newObject, takeOwnership ) );
}

// Doxygen bug - sees this as a separate function
/// @cond DO_NOT_DOCUMENT
template< typename T, typename TBASE >
Wrapper< TBASE > * Group::registerWrapper( std::string const & name,
                                           ViewKey::index_type * const rkey )
{
  m_wrappers.insert( name,
                     (Wrapper< TBASE >::template Factory< T >( name, this ) ).release(),
                     true );

  if( rkey != nullptr )
  {
    *rkey = m_wrappers.getIndex( name );
  }
  Wrapper< TBASE > * const rval = getWrapper< TBASE >( name );
  if( rval->sizedFromParent() == 1 )
  {
    rval->resize( this->size());
  }
  return rval;
}
/// @endcond

template< typename T, typename TBASE >
Wrapper< TBASE > * Group::registerWrapper( ViewKey & viewKey )
{
  ViewKey::index_type index;
  Wrapper< TBASE > * const rval = registerWrapper< T, TBASE >( viewKey.Key(), &index );
  viewKey.setIndex( index );

  return rval;
}


template< typename T >
Wrapper< T > * Group::registerWrapper( std::string const & name,
                                       std::unique_ptr< T > newObject )
{
  m_wrappers.insert( name,
                     new Wrapper< T >( name, this, newObject.release(), true ),
                     true );

  Wrapper< T > * const rval = getWrapper< T >( name );
  if( rval->sizedFromParent() == 1 )
  {
    rval->resize( this->size());
  }
  return rval;
}



template< typename T >
Wrapper< T > * Group::registerWrapper( std::string const & name,
                                       T * newObject,
                                       bool takeOwnership )
{
  m_wrappers.insert( name,
                     new Wrapper< T >( name, this, newObject, takeOwnership ),
                     true );

  Wrapper< T > * const rval = getWrapper< T >( name );
  if( rval->sizedFromParent() == 1 )
  {
    rval->resize( this->size());
  }
  return rval;
}

template< typename T >
T const * Group::GetGroupByPath( string const & path ) const
{
  // needed for getting root correctly with GetGroupByPath("/");
  if( path.empty())
  {
    return group_cast< T const * >( this );
  }

  size_t directoryMarker = path.find( '/' );

  if( directoryMarker == std::string::npos )
  {
    // Target should be a child of this group
    return this->GetGroup< T >( path );
  }
  else
  {
    // Split the path
    string const child = path.substr( 0, directoryMarker );
    string const subPath = path.substr( directoryMarker+1, path.size());

    if( directoryMarker == 0 )            // From root
    {
      if( this->getParent() == nullptr )  // At root
      {
        return this->GetGroupByPath< T >( subPath );
      }
      else                               // Not at root
      {
        return this->getParent()->GetGroupByPath< T >( path );
      }
    }
    else if( child[0] == '.' )
    {
      if( child[1] == '.' )               // '../' = Reverse path
      {
        return this->getParent()->GetGroupByPath< T >( subPath );
      }
      else                               // './' = This path
      {
        return this->GetGroupByPath< T >( subPath );
      }
    }
    else
    {
      return m_subGroups[child]->GetGroupByPath< T >( subPath );
    }
  }
}

} /* end namespace dataRepository */
} /* end namespace geosx */

#endif /* GEOSX_DATAREPOSITORY_GROUP_HPP_ */
