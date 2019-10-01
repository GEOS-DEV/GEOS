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
#include "ObjectCatalog.hpp"
#include "MappedVector.hpp"
#include "RestartFlags.hpp"
#include "Wrapper.hpp"
#include "xmlWrapper.hpp"

#include <iostream>

#ifndef USE_DYNAMIC_CASTING
/// macro definition to specify whether or not to use dynamic_cast
#define USE_DYNAMIC_CASTING 1;
#endif

#ifndef NOCHARTOSTRING_KEYLOOKUP
/// macro definition to enable/disable char * lookups
#define NOCHARTOSTRING_KEYLOOKUP 0
#endif

/* Forward declaration of axom::sidre::Group */
namespace axom
{
namespace sidre
{
class Group;
}
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
   * @name constructors, destructor, copy, move, assignments
   */
  ///@{

  /**
   * @param[in] name the name of this object manager
   * @param[in]
   */
  explicit Group( std::string const & name,
                  Group * const parent );


  /**
   * @brief move constructor
   * @param[in] source source Group
   */
  Group( Group && source );

  /**
   *
   */
  virtual ~Group();



  Group() = delete;
  Group( Group const & source ) = delete;
  Group & operator=( Group const & ) = delete;
  Group & operator=( Group && ) = delete;

  ///@}


  /**
   * @name Static Factory Catalog Functions
   */
  ///@{

  using CatalogInterface = cxx_utilities::CatalogInterface< Group, std::string const &, Group * const >;
  static CatalogInterface::CatalogType & GetCatalog();
  ///@}


  /// returns typeid(*this)
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
    return typeToCheck == get_typeid() ? true : false;
  }

  //START_SPHINX_INCLUDE_REGISTER_GROUP
  /**
   * @name FUNCTION GROUP for functions that add/register a new sub-Group to the current Group
   * @tparam T The type of the Group to add/register. This should be a type that derives from Group.
   * @param[in] name      The name of the group to use as a string key.
   * @param[in] newObject A pointer/unique_ptr to the object that is being registered.
   * @param[in] takeOwnership A flag to indicate whether or not the repository should
   *                          take ownership of the group. This only applies when
   *                          passing a raw pointer. If passing a unique_ptr, then
   *                          the ownership is always transferred.
   * @param keyIndex A KeyIndex object that will be used to specify the name of
   *                 the new group. The index of the KeyIndex will also be set.
   * @param catalogName The catalog name of the new type.
   * @return A pointer to the newly registered/created Group.
   *
   * These registration functions registers (and may create) a Group or class
   * derived from Group as a subgroup of this Group.
   */
  ///@{

  template< typename T = Group >
  T * RegisterGroup( std::string const & name, std::unique_ptr< Group > newObject );

  template< typename T = Group >
  T * RegisterGroup( std::string const & name,
                     T * newObject,
                     bool const takeOwnership );

  template< typename T = Group >
  T * RegisterGroup( std::string const & name )
  {
    return RegisterGroup< T >( name, std::move( std::make_unique< T >( name, this )) );
  }

  template< typename T = Group >
  T * RegisterGroup( subGroupMap::KeyIndex & keyIndex )
  {
    T * rval = RegisterGroup< T >( keyIndex.Key(), std::move( std::make_unique< T >( keyIndex.Key(), this )) );
    keyIndex.setIndex( this->m_subGroups.getIndex( keyIndex.Key()) );
    return rval;
  }

  template< typename T = Group, typename TBASE = Group >
  T * RegisterGroup( std::string const & name, std::string const & catalogName )
  {
    std::unique_ptr< TBASE > newGroup = TBASE::CatalogInterface::Factory( catalogName, name, this );
    return RegisterGroup< T >( name, std::move( newGroup ) );
  }
  ///@}

  //END_SPHINX_INCLUDE_REGISTER_GROUP

  /**
   * @brief       Downcast a Group *
   * @tparam T    Pointer to the type to downcast into
   * @param group A pointer the group to be casted
   * @return a    Pointer to \p T that refers to the downcasted group
   */
  template< typename T >
  static T group_cast( Group * group )
  {
#ifdef USE_DYNAMIC_CASTING
    return dynamic_cast< T >( group );
#else
    return static_cast< T >( group );
#endif
  }

  /**
   * @brief       Downcast a Group const *
   * @tparam T    Pointer to the type to downcast into
   * @param group A pointer the group to be casted
   * @return a    Pointer to \p T that refers to the downcasted group
   */
  template< typename T >
  static T group_cast( Group const * group )
  {
#ifdef USE_DYNAMIC_CASTING
    return dynamic_cast< T >( group );
#else
    return static_cast< T >( group );
#endif
  }

  /**
   * @brief       Downcast this Group
   * @tparam T    Pointer to the type to downcast into
   * @return a    Pointer to \p T that refers to the this
   */
  template< typename T >
  T group_cast()
  {
#ifdef USE_DYNAMIC_CASTING
    return dynamic_cast< T >( this );
#else
    return static_cast< T >( this );
#endif
  }

  /**
   * @brief       Downcast this Group
   * @tparam T    Pointer to the type to downcast into
   * @return a    Pointer to \p T that refers to the this
   */
  template< typename T >
  T group_cast() const
  {
#ifdef USE_DYNAMIC_CASTING
    return dynamic_cast< T >( this );
#else
    return static_cast< T >( this );
#endif
  }

  //START_SPHINX_INCLUDE_GET_GROUP
  /**
   * @name FUNCTION GROUP for GetGroup()
   * @brief Functions to retrieve a sub-Group from the current Group using various lookup methods.
   * @tparam T The type of the Group. This may be a type that derives from Group.
   * @param[in] index The integral index of the group to retrieve.
   * @param[in] name The name/key of the group lookup and retrieve.
   * @param[in] key The KeyIndex to use for the lookup. Note that const-correctness may be
   *                broken if the key is incorrect as @p key will be modified to contain the
   *                correct index.
   * @param[in] path A unix-style string (absolute, relative paths valid) to lookup the Group
   *                 to return.
   * @return A pointer to @p T (or const @p T) that refers to the group retrieved from the
   *         repository.
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
   * cast fails, then a nullptr is returned. If no template parameter is specified then a default
   * type of Group is assumed.
   */
  ///@{

  template< typename T = Group >
  T * GetGroup( localIndex index )
  {
    return group_cast< T * >( m_subGroups[index] );
  }

  template< typename T = Group >
  T const * GetGroup( localIndex index ) const
  {
    return group_cast< T const * >( m_subGroups[index] );
  }

  template< typename T = Group >
  T * GetGroup( string const & name )
  {
    return group_cast< T * >( m_subGroups[name] );
  }

  template< typename T = Group >
  T const * GetGroup( string const & name ) const
  {
    return group_cast< T const * >( m_subGroups[name] );
  }

  template< typename T = Group >
  T * GetGroup( subGroupMap::KeyIndex & key )
  {
    return group_cast< T * >( m_subGroups[key] );
  }

  template< typename T = Group >
  T const * GetGroup( subGroupMap::KeyIndex & key ) const
  {
    return group_cast< T const * >( m_subGroups[key] );
  }

  template< typename T = Group >
  T * GetGroup( subGroupMap::KeyIndex const & key )
  {
    return group_cast< T * >( m_subGroups[key] );
  }

  template< typename T = Group >
  T const * GetGroup( subGroupMap::KeyIndex const & key ) const
  {
    return group_cast< T const * >( m_subGroups[key] );
  }

  template< typename T = Group >
  T const * GetGroupByPath( string const & path ) const;

  template< typename T = Group >
  T * GetGroupByPath( string const & path )
  {
    return const_cast< T * >(const_cast< Group const * >(this)->GetGroupByPath< T >( path ));
  }
  ///@}
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
   * @name FUNCTION GROUP for applyLambdaToGroup()
   * @brief These functions apply the specified lambda function to a container if the container can be
   *        casted to the templated type/s.
   * @tparam CASTTYPE the first type that will be used in the attempted casting of container.
   * @tparam CASTTYPES a variadic list of types that will be used in the attempted casting of
   *                    container.
   * @tparam LAMBDA the type of lambda function to call in the function
   * @param[in] container A pointer to the container which will be passed to the lambda function
   * @param[in] lambda the lambda function to call in the function
   * @return A boolean to indicate whether the lambda was successfully applied to the container.
   *
   * This function is useful when trying to apply a lambda that passes a pointer to an container,
   * but it is desired that the lambda is only executed if the container can be casted to a certain
   * type. The variadic list consisting of CASTTYPE/S will be used recursively to check if the
   * container is able to be casted to the one of these types. The first type in the CASTTYPE/S list
   * will be used to execute the lambda, and the function will return true.
   */
  ///@{

  /** \cond SKIPME */
  template< typename CONTAINERTYPE, typename LAMBDA >
  static bool applyLambdaToContainer( CONTAINERTYPE const * const GEOSX_UNUSED_ARG( group ), LAMBDA && GEOSX_UNUSED_ARG( lambda ) )
  { return false; }

  template< typename CONTAINERTYPE, typename LAMBDA >
  static bool applyLambdaToContainer( CONTAINERTYPE * const GEOSX_UNUSED_ARG( group ), LAMBDA && GEOSX_UNUSED_ARG( lambda ) )
  { return false; }
  /** \endcond */

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
   * @name FUNCTION GROUP for forSubGroups()
   * @brief These functions apply the specified lambda function to a group if the group can be
   *        casted to the templated type/s.
   * @tparam GROUPTYPE The first type that will be used in the attempted casting of group.
   * @tparam GROUPTYPES A variadic list of types that will be used in the attempted casting of
   *                    group.
   * @tparam LAMBDA The type of lambda function to call in the function
   * @param[in] lambda The lambda function to call in the function
   * @param[in] subgroupNames Optional list of subgroup names to apply the lambda to
   *
   * These functions loop over sub-groups and executes a lambda that uses the sub-group as an
   * argument. The lambda is only executed if the group can be casted to a certain type specified
   * by the GROUPTYPE/S pack. The variadic list consisting of GROUPTYPE/S will be used recursively
   * to check if the group is able to be casted to the one of these types. The first type in the
   * GROUPTYPE/S list will be used to execute the lambda, and the next sub-group will be processed.
   */
  ///@{
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
   * @name FUNCTION GROUP for forWrappers()
   * @brief These functions apply the specified lambda function to a Wrapper if the Wrapper can be
   *        casted to the templated type/s.
   * @tparam TYPE The first type that will be used in the attempted casting of Wrapper.
   * @tparam TYPES A variadic list of types that will be used in the attempted casting of Wrapper.
   * @tparam LAMBDA The type of lambda function to call in the function
   * @param[in] lambda The lambda function to call in the function
   *
   * These functions loop over the Wrappers contained in this group, and executes a lambda that
   * uses the Wrapper as an argument. The lambda is only executed if the Wrapper can be casted to
   * a certain type specified by the TYPE/S pack. The variadic list consisting of
   * TYPE/S will be used recursively to check if the Wrapper is able to be casted to the
   * one of these types. The first type in the WRAPPERTYPE/S list will be used to execute the
   * lambda, and the next Wrapper will be processed.
   */
  ///@{
  template< typename LAMBDA >
  void forWrappers( LAMBDA lambda )
  {
    for( auto & wrapperIter : m_wrappers )
    {
      lambda( *wrapperIter.second );
    }
  }

  template< typename LAMBDA >
  void forWrappers( LAMBDA lambda ) const
  {
    for( auto const & wrapperIter : m_wrappers )
    {
      lambda( *wrapperIter.second );
    }
  }

  template< typename TYPE, typename ... TYPES, typename LAMBDA >
  void forWrappers( LAMBDA lambda )
  {
    for( auto & wrapperIter : m_wrappers )
    {
      applyLambdaToContainer< WrapperBase, Wrapper< TYPE >, Wrapper< TYPES >... >( wrapperIter.second,
                                                                                   std::forward< LAMBDA >( lambda ));
    }
  }

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


  void Initialize( Group * const group );

  virtual void InitializationOrder( string_array & order );


  void InitializePostInitialConditions( Group * const group );

  //START_SPHINX_INCLUDE_REGISTER_WRAPPER
  /**
   * @name FUNCTION GROUP for functions that add/register a new Wrapper to this Group.
   * @tparam T The type of the Wrapper to add/register.
   * @tparam TBASE The type of the Wrapper to add/register.
   * @param[in] name      The name of the group to use as a string key.
   * @param[in] rkey      A pointer to a index type that will be filled with the new
   *                      Wrappers index in this Group.
   * @param[in] viewKey The KeyIndex that contains the name of the new Wrapper.
   * @param[in] type The runtime type to wrap in the new Wrapper.
   * @param[in] newObject A pointer/unique_ptr to the object that is being registered.
   * @param[in] takeOwnership A flag to indicate whether or not the repository should
   *                          take ownership of the group. This only applies when
   *                          passing a raw pointer. If passing a unique_ptr, then
   *                          the ownership is always transferred.
   * @param[in] wrapper A pointer to the an existing wrapper.
   * @return A pointer to the newly registered/created Wrapper.
   *
   * These registration functions register (and may create) a Wrapper's around objects.
   */
  ///@{
  template< typename T, typename TBASE=T >
  Wrapper< TBASE > * registerWrapper( std::string const & name,
                                      wrapperMap::KeyIndex::index_type * const rkey = nullptr );

  template< typename T, typename TBASE=T >
  Wrapper< TBASE > * registerWrapper( Group::wrapperMap::KeyIndex & viewKey );


  WrapperBase * registerWrapper( std::string const & name,
                                 rtTypes::TypeIDs const & type );

  template< typename T >
  Wrapper< T > * registerWrapper( std::string const & name,
                                  std::unique_ptr< T > newObject );

  template< typename T >
  Wrapper< T > * registerWrapper( std::string const & name,
                                  T * newObject,
                                  bool takeOwnership );

  WrapperBase * registerWrapper( string const & name,
                                 WrapperBase * const wrapper );

  ///@}
  //END_SPHINX_INCLUDE_REGISTER_WRAPPER

  /**
   * @brief Removes a Wrapper from this group.
   * @param name The name of the Wrapper to remove from this group.
   */
  void deregisterWrapper( string const & name );

  void PrintDataHierarchy( integer indent = 0 );

  virtual Group * CreateChild( string const & childKey, string const & childName );

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

  /**
   * This function is used to build a complete datastructure for schema generation
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
   * This function is used to expand any catalogs in the data structure
   */
  virtual void ExpandObjectCatalogs() {}

  /**
   * This function is used to inform the schema generator of any
   * deviations between the xml and GEOS data structures.
   */
  virtual void SetSchemaDeviations( xmlWrapper::xmlNode GEOSX_UNUSED_ARG( schemaRoot ),
                                    xmlWrapper::xmlNode GEOSX_UNUSED_ARG( schemaParent ),
                                    integer GEOSX_UNUSED_ARG( documentationType ) ) {}

  virtual void RegisterDataOnMeshRecursive( Group * const MeshBodies );

  virtual void RegisterDataOnMesh( Group * const GEOSX_UNUSED_ARG( MeshBody ) ) {}

  virtual localIndex PackSize( string_array const & wrapperNames,
                               integer const recursive ) const;

  virtual localIndex PackSize( string_array const & wrapperNames,
                               arrayView1d< localIndex const > const & packList,
                               integer const recursive ) const;

  virtual localIndex Pack( buffer_unit_type * & buffer,
                           string_array const & wrapperNames,
                           integer const recursive ) const;

  virtual localIndex Pack( buffer_unit_type * & buffer,
                           string_array const & wrapperNames,
                           arrayView1d< localIndex const > const & packList,
                           integer const recursive ) const;

  virtual localIndex Unpack( buffer_unit_type const * & buffer,
                             arrayView1d< localIndex > & packList,
                             integer const recursive );


  //***********************************************************************************************

  //START_SPHINX_INCLUDE_GET_WRAPPER
  WrapperBase const * getWrapperBase( indexType const index ) const
  { return m_wrappers[index]; }

  WrapperBase * getWrapperBase( indexType const index )
  { return m_wrappers[index]; }

  WrapperBase const * getWrapperBase( std::string const & name ) const
  { return m_wrappers[name]; }

  WrapperBase * getWrapperBase( std::string const & name )
  { return m_wrappers[name]; }

  WrapperBase const * getWrapperBase( wrapperMap::KeyIndex const & keyIndex ) const
  { return m_wrappers[keyIndex]; }

  WrapperBase * getWrapperBase( wrapperMap::KeyIndex const & keyIndex )
  { return m_wrappers[keyIndex]; }


  template< typename T, typename LOOKUP_TYPE >
  Wrapper< T > const * getWrapper( LOOKUP_TYPE const & index ) const
  {
#ifdef USE_DYNAMIC_CASTING
    return dynamic_cast< Wrapper< T > const * >( (m_wrappers[index]) );
#else
    return static_cast< Wrapper< T > const * >( (m_wrappers[index]) );
#endif
  }

  template< typename T, typename LOOKUP_TYPE >
  Wrapper< T > * getWrapper( LOOKUP_TYPE const & index )
  { return const_cast< Wrapper< T > * >( const_cast< Group const * >(this)->getWrapper< T >( index ) ); }

  template< typename T >
  Wrapper< T > const * getWrapper( char const * const key ) const
  { return getWrapper< T >( string( key ) ); }

  template< typename T >
  Wrapper< T > * getWrapper( char const * const key )
  { return getWrapper< T >( string( key ) ); }


  indexType getWrapperIndex( std::string const & name ) const
  {
    return m_wrappers.getIndex( name );
  }


  template< typename T, typename WRAPPEDTYPE=T, typename LOOKUP_TYPE >
  typename std::enable_if< std::is_same< T, WRAPPEDTYPE >::value, T const & >::type
  getReference( LOOKUP_TYPE const & lookup ) const
  {
    Wrapper< WRAPPEDTYPE > const * wrapper = getWrapper< WRAPPEDTYPE >( lookup );
    if( wrapper == nullptr )
    {
      if( hasView( lookup ) )
      {
        GEOS_ERROR( "call to getWrapper results in nullptr but a view exists. Most likely given the incorrect type. lookup : " << lookup );
      }
      GEOS_ERROR( "call to getWrapper results in nullptr and a view does not exist. lookup : " << lookup );
    }

    return wrapper->reference();
  }

  template< typename T, typename WRAPPEDTYPE=T, typename LOOKUP_TYPE >
  typename std::enable_if< !std::is_same< T, WRAPPEDTYPE >::value, T const & >::type
  getReference( LOOKUP_TYPE const & lookup ) const
  {
    static_assert( std::is_base_of< WRAPPEDTYPE, T >::value, "incorrect template arguments" );
    Wrapper< WRAPPEDTYPE > const * wrapper = getWrapper< WRAPPEDTYPE >( lookup );
    if( wrapper == nullptr )
    {
      if( hasView( lookup ) )
      {
        GEOS_ERROR( "call to getWrapper results in nullptr but a view exists. Most likely given the incorrect type. lookup : " << lookup );
      }
      GEOS_ERROR( "call to getWrapper results in nullptr and a view does not exist. lookup : " << lookup );
    }

#ifdef USE_DYNAMIC_CASTING
    return dynamic_cast< T const & >( wrapper->reference() );
#else
    return static_cast< T const & >( wrapper->reference() );
#endif
  }


  template< typename T, typename WRAPPEDTYPE=T, typename LOOKUP_TYPE >
  T & getReference( LOOKUP_TYPE const & lookup )
  { return const_cast< T & >( const_cast< const Group * >(this)->template getReference< T, WRAPPEDTYPE, LOOKUP_TYPE >( lookup ) ); }

  template< typename T, typename WRAPPEDTYPE=T >
  T const & getReference( char const * const name ) const
  { return getReference< T, WRAPPEDTYPE >( string( name ) ); }

  template< typename T, typename WRAPPEDTYPE=T >
  T & getReference( char const * const name )
  { return const_cast< T & >( const_cast< const Group * >(this)->getReference< T, WRAPPEDTYPE >( name ) ); }



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

  template< typename T, typename LOOKUP_TYPE >
  T * getPointer( LOOKUP_TYPE const & lookup )
  { return const_cast< T * >( const_cast< Group const * >(this)->getPointer< T >( lookup )); }

  template< typename T >
  T const * getPointer( char const * const name ) const
  { return getPointer< T >( string( name ) ); }

  template< typename T >
  T * getPointer( char const * const name )
  { return getPointer< T >( string( name ) ); }
  //END_SPHINX_INCLUDE_GET_WRAPPER


  bool hasGroup( std::string const & name ) const
  {
    return (m_subGroups[name] != nullptr);
  }

  template< typename LOOKUP_TYPE >
  bool hasView( LOOKUP_TYPE const & lookup ) const
  {
    return (m_wrappers[lookup] != nullptr);
  }

  inline const string getName() const
  {
    return m_name;
  }


  virtual void resize( localIndex const newsize );

  virtual void reserve( indexType const newsize );


  inline localIndex capacity() const
  {
    return m_capacity;
  }

  inline localIndex size() const
  {
    return m_size;
  }

  axom::sidre::Group * getSidreGroup()
  {
#ifdef GEOSX_USE_ATK
    return m_sidreGroup;
#else
    return nullptr;
#endif
  }

  axom::sidre::Group const * getSidreGroup() const
  {
#ifdef GEOSX_USE_ATK
    return m_sidreGroup;
#else
    return nullptr;
#endif
  }

  static axom::sidre::Group * setSidreGroup( string const & name,
                                             Group * const parent );

  Group * getParent()             { return m_parent; }
  Group const * getParent() const { return m_parent; }

  Group * setParent( Group * const parent )
  {
    m_parent = parent;
#ifdef GEOSX_USE_ATK
    m_sidreGroup = m_parent->getSidreGroup();
#endif

    return m_parent;
  }

  localIndex getIndexInParent() const
  {
    return m_parent->GetSubGroups().getIndex( this->m_name );
  }


  wrapperMap const & wrappers() const
  {
    return m_wrappers;
  }

  wrapperMap & wrappers()
  {
    return m_wrappers;
  }

  RestartFlags getRestartFlags() const { return m_restart_flags; }

  void setRestartFlags( RestartFlags flags ) { m_restart_flags = flags; }

  InputFlags getInputFlags() const { return m_input_flags; }

  void setInputFlags( InputFlags flags ) { m_input_flags = flags; }

  void prepareToWrite();

  void finishWriting() const;

  void prepareToRead();

  void finishReading();

  void postRestartInitializationRecursive( Group * const domain );


protected:
  /**
   * @brief Post processing of the input values.
   */
  virtual void PostProcessInput() {}

  virtual void InitializePreSubGroups( Group * const GEOSX_UNUSED_ARG( group ) ) {}

  virtual void InitializePostSubGroups( Group * const GEOSX_UNUSED_ARG( group ) ) {}

  virtual void InitializePostInitialConditions_PreSubGroups( Group * const GEOSX_UNUSED_ARG( group ) ) {}

  virtual void InitializePostInitialConditions_PostSubGroups( Group * const GEOSX_UNUSED_ARG( group ) ) {}

  virtual void postRestartInitialization( Group * const GEOSX_UNUSED_ARG( domain ) ) {}

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

  RestartFlags m_restart_flags; ///< Restart flag for this group...and
                                ///< subsequently all wrappers in this group
  InputFlags m_input_flags;     ///< Input flag for this group

#ifdef GEOSX_USE_ATK
  /// Pointer to the sidre group that mirrors this group
  axom::sidre::Group * m_sidreGroup;
#endif


};

using GroupKey = Group::subGroupMap::KeyIndex;
using ViewKey = Group::wrapperMap::KeyIndex;



template< typename T >
T * Group::RegisterGroup( std::string const & name,
                          std::unique_ptr< Group > newObject )
{
#ifdef USE_DYNAMIC_CASTING
  return dynamic_cast< T * >( m_subGroups.insert( name, newObject.release(), true ) );
#else
  return static_cast< T * >( m_subGroups.insert( name, newObject.release(), true ) );
#endif
}


template< typename T >
T * Group::RegisterGroup( std::string const & name,
                          T * newObject,
                          bool const takeOwnership )
{
#ifdef USE_DYNAMIC_CASTING
  return dynamic_cast< T * >( m_subGroups.insert( name, newObject, takeOwnership ) );
#else
  return static_cast< T * >( m_subGroups.insert( name, newObject, takeOwnership ) );
#endif
}


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
  if( rval->sizedFromParent() == 1 && rval->shouldResize())
  {
    rval->resize( this->size());
  }
  return rval;
}

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
  if( rval->sizedFromParent() == 1 && rval->shouldResize())
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
