/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/** @file */

#ifndef MANAGEDGROUP_H_
#define MANAGEDGROUP_H_

#include <iostream>
#include <mpi.h>


#include "ObjectCatalog.hpp"
#include "ViewWrapper.hpp"
#include "RestartFlags.hpp"

#include "depricated/Common.h"

#include "MappedVector.hpp"

#include "fileIO/xmlWrapper.hpp"
#include "dataRepository/InputFlags.hpp"

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
 * namespace to encapsulate functions in simulation tools
 */
namespace geosx
{
namespace dataRepository
{

/// the default key type
using keyType = string;

/// the default index type
using indexType = localIndex;

/**
 * class that encapsulates and manages a collection of DataObjects. Can be
 * considered a "node" in a hierarchy of managers that represent physical groupings of data/
 */
class ManagedGroup
{
public:
  /// the type of MappedVector to use for the subGroup collection
  using subGroupMap = MappedVector< ManagedGroup, ManagedGroup*, keyType, indexType  >;

  /// type of the MappedVector to use for the collection of wrappers.
  using viewWrapperMap = MappedVector< ViewWrapperBase, ViewWrapperBase*, keyType, indexType  >;

  /**
   * @name constructors, destructor, copy, move, assignments
   */
  ///@{

  /**
   * @param[in] name the name of this object manager
   * @param[in]
   */
  explicit ManagedGroup( std::string const & name,
                         ManagedGroup * const parent );


  /**
   * @brief move constructor
   * @param[in] source source ManagedGroup
   */
  ManagedGroup( ManagedGroup&& source );

  /**
   *
   */
  virtual ~ManagedGroup();



  ManagedGroup() = delete;
  ManagedGroup( ManagedGroup const & source ) = delete;
  ManagedGroup& operator=( ManagedGroup const & ) = delete;
  ManagedGroup& operator=(ManagedGroup&&) = delete;

  ///@}


  /**
   * @name Static Factory Catalog Functions
   */
  ///@{

  using CatalogInterface = cxx_utilities::CatalogInterface< ManagedGroup, std::string const &, ManagedGroup * const >;
  static CatalogInterface::CatalogType& GetCatalog();
  ///@}


  /// returns typeid(*this)
  virtual const std::type_info& get_typeid() const
  {
    return typeid(*this);
  }

  /**
   * @brief Check a type_info against the type_info of this ManagedGroup
   * @param typeToCheck value to check against
   * @return true of types are the same, false if not
   */
  bool CheckTypeID( std::type_info const & typeToCheck ) const
  {
    return typeToCheck == get_typeid() ? true : false;
  }


  /**
   * @brief Register a new subgroup
   * @tparam T        The type of the group to register
   * @param name      The repository name of the group
   * @param newObject A unique_ptr to the object that is being registered.
   * @return pointer to the subgroup
   *
   * This registration function takes an ManagedGroup or class derived from
   * ManagedGroup and registers it as a subgroup of this ManagedGroup. The
   * data repository takes ownership of the object that is passed in through
   * a unique_ptr.
   */
  template< typename T = ManagedGroup >
  T * RegisterGroup( std::string const & name, std::unique_ptr<ManagedGroup> newObject );

  /**
   *
   * @brief Register a new subgroup
   * @tparam T        The type of the group to register
   * @param name      The repository name of the group
   * @param newObject Pointer to the object tp register
   * @param takeOwnership flag to indicate whether or not the repository should
   *                      take ownership of the group.
   * @return pointer to the subgroup
   *
   * This registration function creates and registers a ManagedGroup or class
   * derived from ManagedGroup and registers it as a subgroup of this
   * ManagedGroup. The data repository takes ownership of the new object if
   * \p takeOwnership is true.
   */
  template< typename T = ManagedGroup >
  T * RegisterGroup( std::string const & name,
                     T * newObject,
                     bool const takeOwnership );

  /**
   * @brief      Register a new subgroup
   * @tparam T   The type of the group to register
   * @param name The repository name of the group
   * @return     pointer to the subgroup
   *
   * This registration function creates and registers a ManagedGroup or class
   * derived from ManagedGroup and registers it as a subgroup of this
   * ManagedGroup. The data repository takes ownership of the new object.
   */
  template< typename T = ManagedGroup >
  T * RegisterGroup( std::string const & name )
  {
    return RegisterGroup<T>( name, std::move(std::make_unique< T >( name, this )) );
  }

  /**
   * @brief          Register a new subgroup
   * @tparam T       The type of the group to register
   * @param keyIndex A KeyIndex object that will be used to specify the name of
   *                 the new group.
   * @return         pointer to the subgroup
   *
   * This registration function creates and registers a ManagedGroup or class
   * derived from ManagedGroup and registers it as a subgroup of this
   * ManagedGroup. The data repository takes ownership of the new object.
   */
  template< typename T = ManagedGroup >
  T * RegisterGroup( subGroupMap::KeyIndex & keyIndex )
  {
    T * rval = RegisterGroup<T>( keyIndex.Key(), std::move(std::make_unique< T >( keyIndex.Key(), this )) );
    keyIndex.setIndex( this->m_subGroups.getIndex(keyIndex.Key()) );
    return rval;
  }

  /**
   * @brief             Register a new subgroup
   * @tparam T          The type of the group to register
   * @tparam TBASE      The type of the group that contains an object catalog
   *                    for the creation of a new \p T
   * @param name        The repository name of the new group
   * @param catalogName The catalog name of the type
   * @return            A pointer to the subgroup
   *
   * This registration function creates and registers a ManagedGroup or class
   * derived from ManagedGroup and registers it as a subgroup of this
   * ManagedGroup. The \p TBASE type provide the object catalog that is used
   * to create the new group. The data repository takes ownership of the new
   * object.
   *
   */
  template< typename T = ManagedGroup, typename TBASE = ManagedGroup >
  T * RegisterGroup( std::string const & name, std::string const & catalogName )
  {
    std::unique_ptr<TBASE> newGroup = TBASE::CatalogInterface::Factory(catalogName, name, this );
    return RegisterGroup<T>( name, std::move(newGroup) );
  }

  /**
   * @brief       Downcast a ManagedGroup *
   * @tparam T    Pointer to the type to downcast into
   * @param group A pointer the group to be casted
   * @return a    Pointer to \p T that refers to the downcasted group
   */
  template< typename T >
  static T group_cast( ManagedGroup * group )
  {
#ifdef USE_DYNAMIC_CASTING
    return dynamic_cast<T>( group );
#else
    return static_cast<T>( group );
#endif
  }

  /**
   * @brief       Downcast a ManagedGroup const *
   * @tparam T    Pointer to the type to downcast into
   * @param group A pointer the group to be casted
   * @return a    Pointer to \p T that refers to the downcasted group
   */
  template< typename T >
  static T group_cast( ManagedGroup const * group )
  {
#ifdef USE_DYNAMIC_CASTING
    return dynamic_cast<T>( group );
#else
    return static_cast<T>( group );
#endif
  }

  /**
   * @brief       Downcast this ManagedGroup
   * @tparam T    Pointer to the type to downcast into
   * @return a    Pointer to \p T that refers to the this
   */
  template< typename T >
  T group_cast()
  {
#ifdef USE_DYNAMIC_CASTING
    return dynamic_cast<T>( this );
#else
    return static_cast<T>( this );
#endif
  }

  /**
   * @brief       Downcast this ManagedGroup
   * @tparam T    Pointer to the type to downcast into
   * @return a    Pointer to \p T that refers to the this
   */
  template< typename T >
  T group_cast() const
  {
#ifdef USE_DYNAMIC_CASTING
    return dynamic_cast<T>( this );
#else
    return static_cast<T>( this );
#endif
  }


  /**
   * @brief Get a ManagedGroup from the repository using direct index
   * @param index integral index to use for as a lookup
   * @return A pointer to \p T that refers to the group retrieved from the
   *         repository.
   */
  template< typename T = ManagedGroup >
  T * GetGroup( localIndex index )
  {
    return group_cast<T*>(m_subGroups[index]);
  }

  /**
   * @brief Get a ManagedGroup from the repository using direct index
   * @param index integral index to use for as a lookup
   * @return A pointer to const \p T that refers to the group retrieved from the
   *         repository.
   */
  template< typename T = ManagedGroup >
  T const * GetGroup( localIndex index ) const
  {
    return group_cast<T const *>(m_subGroups[index]);
  }

  /**
   * @brief Get a ManagedGroup from the repository using string lookup
   * @param name Name of the group to retrieve
   * @return A pointer to \p T that refers to the group retrieved from the
   *         repository.
   */
  template< typename T = ManagedGroup >
  T * GetGroup( string const & name )
  {
    return group_cast<T *>(m_subGroups[name]);
  }

  /**
   * @brief Get a ManagedGroup from the repository using string lookup
   * @param name Name of the group to retrieve
   * @return A pointer to const \p T that refers to the group retrieved from the
   *         repository.
   */
  template< typename T = ManagedGroup >
  T const * GetGroup( string const & name ) const
  {
    return group_cast<T const *>(m_subGroups[name]);
  }

  /**
   * @brief Get a ManagedGroup from the repository using KeyIndex lookup
   * @param key the KeyIndex to use for the lookup
   * @return A pointer to \p T that refers to the group retrieved from the
   *         repository.
   */
  template< typename T = ManagedGroup >
  T * GetGroup( subGroupMap::KeyIndex & key )
  {
    return group_cast<T *>(m_subGroups[key]);
  }

  /**
   * @brief Get a ManagedGroup from the repository using KeyIndex lookup
   * @param key the KeyIndex to use for the lookup
   * @return A pointer to const \p T that refers to the group retrieved from
   *         the repository.
   */
  template< typename T = ManagedGroup >
  T const * GetGroup( subGroupMap::KeyIndex & key ) const
  {
    return group_cast<T const *>(m_subGroups[key]);
  }

  /**
   * @brief Get a ManagedGroup from the repository using KeyIndex lookup
   * @param key the KeyIndex to use for the lookup. Note that
   *            const-correctness may be broken if the key is incorrect.
   * @return A pointer to \p T that refers to the group retrieved from the
   *         repository.
   */
  template< typename T = ManagedGroup >
  T * GetGroup( subGroupMap::KeyIndex const & key )
  {
    return group_cast<T *>(m_subGroups[key]);
  }

  /**
   * @brief Get a ManagedGroup from the repository using KeyIndex lookup
   * @param key the KeyIndex to use for the lookup. Note that
   *            const-correctness may be broken if the key is incorrect.
   * @return A pointer to const \p T that refers to the group retrieved from
   *         the repository.
   */
  template< typename T = ManagedGroup >
  T const * GetGroup( subGroupMap::KeyIndex const & key ) const
  {
    return group_cast<T const *>(m_subGroups[key]);
  }

  /**
   * @brief This will grab the pointer to an object in the data structure
   * @param[in] path a unix-style string (absolute, relative paths valid)
   * @return A pointer to const \p T that refers to the target group.
   */
  template< typename T = ManagedGroup >
  T const * GetGroupByPath( string const & path ) const
  {
    // needed for getting root correctly with GetGroupByPath("/");
    if (path.empty())
    {
      return group_cast<T const *>( this );
    }

    size_t directoryMarker = path.find('/');

    if (directoryMarker == std::string::npos)
    {
      // Target should be a child of this group
      return this->GetGroup<T>(path);
    }
    else
    {
      // Split the path
      string const child = path.substr(0, directoryMarker);
      string const subPath = path.substr(directoryMarker+1, path.size());

      if (directoryMarker == 0)            // From root
      {
        if (this->getParent() == nullptr)  // At root
        {
          return this->GetGroupByPath<T>(subPath);
        }
        else                               // Not at root
        {
          return this->getParent()->GetGroupByPath<T>(path);
        }
      }
      else if (child[0] == '.')
      {
        if (child[1] == '.')               // '../' = Reverse path
        {
          return this->getParent()->GetGroupByPath<T>(subPath);
        }
        else                               // './' = This path
        {
          return this->GetGroupByPath<T>(subPath);
        }
      }
      else
      {
        return m_subGroups[child]->GetGroupByPath<T>(subPath);
      }
    }
  }

  /**
   * @brief This will grab the pointer to an object in the data structure
   * @param[in] path a unix-style string (absolute, relative paths valid)
   * @return A pointer to \p T that refers to the target group.
   */
  template< typename T = ManagedGroup >
  T * GetGroupByPath( string const & path )
  {
    return const_cast<T *>(const_cast< ManagedGroup const * >(this)->GetGroupByPath<T>(path));
  }

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
   * @brief return the number of sub groups in this ManagedGroup
   * @return number of sub groups in this ManagedGroup
   */
  localIndex numSubGroups() const { return m_subGroups.size(); }

  template< typename T = ManagedGroup, typename LAMBDA >
  void forSubGroups( LAMBDA lambda )
  {
    for( auto& subGroupIter : m_subGroups )
    {
#ifdef USE_DYNAMIC_CASTING
      T * subGroup = dynamic_cast<T *>( subGroupIter.second );
#else
      T * subGroup = static_cast<T *>( subGroupIter.second );
#endif
      lambda( subGroup );
    }
  }

  /**
   *
   * @param lambda
   */
  template< typename T = ManagedGroup, typename LAMBDA >
  void forSubGroups( LAMBDA lambda ) const
  {
    for( auto const & subGroupIter : m_subGroups )
    {
#ifdef USE_DYNAMIC_CASTING
      T const * subGroup = dynamic_cast<T const *>( subGroupIter.second );
#else
      T const * subGroup = static_cast<T const *>( subGroupIter.second );
#endif
      lambda( subGroup );
    }
  }

  template< typename T = ViewWrapperBase, typename LAMBDA >
  void forViewWrappers( LAMBDA lambda )
  {
    for( auto& wrapperIter : m_wrappers )
    {
#ifdef USE_DYNAMIC_CASTING
      T & wrapper = dynamic_cast<T &>( *wrapperIter.second );
#else
      T & wrapper = static_cast<T &>( *wrapperIter.second );
#endif
      lambda( wrapper );
    }
  }

  template< typename T = ViewWrapperBase, typename LAMBDA >
  void forViewWrappers( LAMBDA lambda ) const
  {
    for( auto const & wrapperIter : m_wrappers )
    {
#ifdef USE_DYNAMIC_CASTING
      T const & wrapper = dynamic_cast<T const &>( *wrapperIter.second );
#else
      T const & wrapper = static_cast<T const &>( *wrapperIter.second );
#endif
      lambda( wrapper );
    }
  }

  template< typename Wrapped, typename LAMBDA >
  void forViewWrappersByType(LAMBDA lambda)
  {
    for( auto & wrapperIter : m_wrappers )
    {
      ViewWrapper<Wrapped> * const wrapper = ViewWrapper<Wrapped>::cast(wrapperIter.second);
      if ( wrapper != nullptr )
      {
        lambda(*wrapper);
      }
    }
  }

  template< typename Wrapped, typename LAMBDA >
  void forViewWrappersByType(LAMBDA lambda) const
  {
    for( auto const & wrapperIter : m_wrappers )
    {
      ViewWrapper<Wrapped> const * const wrapper = ViewWrapper<Wrapped>::cast(wrapperIter.second);
      if ( wrapper != nullptr )
      {
        lambda(*wrapper);
      }
    }
  }

  void Initialize( ManagedGroup * const group );

  virtual void InitializationOrder( string_array & order );


  void InitializePostInitialConditions( ManagedGroup * const group );

  template< typename T, typename TBASE=T >
  ViewWrapper<TBASE> *
  RegisterViewWrapper( std::string const & name,
                       viewWrapperMap::KeyIndex::index_type * const rkey = nullptr );

  template< typename T, typename TBASE=T >
  ViewWrapper<TBASE> *
  RegisterViewWrapper( ManagedGroup::viewWrapperMap::KeyIndex & viewKey );


  ViewWrapperBase * RegisterViewWrapper( std::string const & name,
                                         rtTypes::TypeIDs const & type );

  template< typename T >
  ViewWrapper<T> * RegisterViewWrapper( std::string const & name,
                                        std::unique_ptr<T> newObject );

  template< typename T >
  ViewWrapper<T> * RegisterViewWrapper( std::string const & name,
                                        T * newObject,
                                        bool takeOwnership );

  /**
   * @brief Register a ViewWrapper into this ManagedGroup
   * @param[in] name the key name to use for this new wrapper
   * @param[in] wrapper a pointer to the new wrapper
   * @return a ViewWrapperBase pointer that holds the address of the new wrapper
   */
  ViewWrapperBase * RegisterViewWrapper( string const & name,
                                         ViewWrapperBase * const wrapper );

  void DeregisterViewWrapper( string const & name );


  void PrintDataHierarchy(integer indent = 0);

  virtual ManagedGroup * CreateChild( string const & childKey, string const & childName );

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
  void GenerateDataStructureSkeleton(integer const level)
  {
    ExpandObjectCatalogs();
    std::string indent( level*2, ' ');

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
  virtual void SetSchemaDeviations(xmlWrapper::xmlNode schemaRoot,
                                   xmlWrapper::xmlNode schemaParent) {}

  virtual void RegisterDataOnMeshRecursive( ManagedGroup * const MeshBodies );

  virtual void RegisterDataOnMesh( ManagedGroup * const MeshBody ) {}
  
  virtual localIndex PackSize( string_array const & wrapperNames,
                               integer const recursive ) const;

  virtual localIndex PackSize( string_array const & wrapperNames,
                               arrayView1d<localIndex const> const & packList,
                               integer const recursive ) const;

  virtual localIndex Pack( buffer_unit_type * & buffer,
                           string_array const & wrapperNames,
                           integer const recursive ) const;

  virtual localIndex Pack( buffer_unit_type * & buffer,
                           string_array const & wrapperNames,
                           arrayView1d<localIndex const> const & packList,
                           integer const recursive ) const;

  virtual localIndex Unpack( buffer_unit_type const *& buffer,
                             arrayView1d<localIndex> & packList,
                             integer const recursive );


  //***********************************************************************************************

  ViewWrapperBase const * getWrapperBase( indexType const index ) const
  { return m_wrappers[index]; }

  ViewWrapperBase * getWrapperBase( indexType const index )
  { return m_wrappers[index]; }

  ViewWrapperBase const * getWrapperBase( std::string const & name ) const
  { return m_wrappers[name]; }

  ViewWrapperBase * getWrapperBase( std::string const & name )
  { return m_wrappers[name]; }

  ViewWrapperBase const * getWrapperBase( viewWrapperMap::KeyIndex const & keyIndex ) const
  { return m_wrappers[keyIndex]; }

  ViewWrapperBase * getWrapperBase( viewWrapperMap::KeyIndex const & keyIndex )
  { return m_wrappers[keyIndex]; }


  template< typename T, typename LOOKUP_TYPE >
  ViewWrapper<T> const * getWrapper( LOOKUP_TYPE const & index ) const
  {
#ifdef USE_DYNAMIC_CASTING
    return dynamic_cast< ViewWrapper<T> const * >( (m_wrappers[index]) );
#else
    return static_cast< ViewWrapper<T> const * >( (m_wrappers[index]) );
#endif
  }

  template< typename T, typename LOOKUP_TYPE >
  ViewWrapper<T> * getWrapper( LOOKUP_TYPE const & index )
  { return const_cast<ViewWrapper<T> *>( const_cast< ManagedGroup const *>(this)->getWrapper<T>( index ) ); }

  template< typename T >
  ViewWrapper<T> const * getWrapper( char const * const key ) const
  { return getWrapper<T>( string( key ) ); }

  template< typename T >
  ViewWrapper<T> * getWrapper( char const * const key )
  { return getWrapper<T>( string( key ) ); }


  indexType getWrapperIndex( std::string const & name ) const
  {
    return m_wrappers.getIndex(name);
  }


  template< typename T, typename WRAPPEDTYPE=T, typename LOOKUP_TYPE>
  typename std::enable_if< std::is_same<T,WRAPPEDTYPE>::value, T const & >::type
  getReference( LOOKUP_TYPE const & lookup ) const
  {
    ViewWrapper<WRAPPEDTYPE> const * wrapper = getWrapper<WRAPPEDTYPE>(lookup);
    if( wrapper == nullptr )
    {
      if ( hasView(lookup) )
      {
        GEOS_ERROR( "call to getWrapper results in nullptr but a view exists. Most likely given the incorrect type. lookup : " << lookup );
      }
      GEOS_ERROR( "call to getWrapper results in nullptr and a view does not exist. lookup : " << lookup );
    }

    return wrapper->reference();
  }

  template< typename T, typename WRAPPEDTYPE=T, typename LOOKUP_TYPE >
  typename std::enable_if< !std::is_same<T,WRAPPEDTYPE>::value, T const & >::type
  getReference( LOOKUP_TYPE const & lookup ) const
  {
    static_assert( std::is_base_of<WRAPPEDTYPE,T>::value,"incorrect template arguments");
    ViewWrapper<WRAPPEDTYPE> const * wrapper = getWrapper<WRAPPEDTYPE>(lookup);
    if( wrapper == nullptr )
    {
      if ( hasView(lookup) )
      {
        GEOS_ERROR( "call to getWrapper results in nullptr but a view exists. Most likely given the incorrect type. lookup : " << lookup );
      }
      GEOS_ERROR( "call to getWrapper results in nullptr and a view does not exist. lookup : " << lookup );
    }

#ifdef USE_DYNAMIC_CASTING
    return dynamic_cast<T const &>( wrapper->reference() );
#else
    return static_cast<T const &>( wrapper->reference() );
#endif
  }


  template< typename T, typename WRAPPEDTYPE=T, typename LOOKUP_TYPE >
  T & getReference( LOOKUP_TYPE const & lookup )
  { return const_cast<T&>( const_cast<const ManagedGroup*>(this)->template getReference<T,WRAPPEDTYPE,LOOKUP_TYPE>( lookup ) ); }

  template< typename T, typename WRAPPEDTYPE=T >
  T const & getReference( char const * const name ) const
  { return getReference<T, WRAPPEDTYPE>( string(name) ); }

  template< typename T, typename WRAPPEDTYPE=T >
  T & getReference( char const * const name )
  { return const_cast<T&>( const_cast<const ManagedGroup*>(this)->getReference<T,WRAPPEDTYPE>( name ) ); }



  template< typename T, typename LOOKUP_TYPE >
  T const * getPointer( LOOKUP_TYPE const & lookup ) const
  {
    T const * rval = nullptr;
    ViewWrapper<T> const * wrapper = getWrapper<T>(lookup);
    if( wrapper != nullptr )
    {
      rval = wrapper->getPointer();
    }
    return rval;
  }

  template< typename T, typename LOOKUP_TYPE >
  T * getPointer( LOOKUP_TYPE const & lookup )
  { return const_cast<T *>( const_cast<ManagedGroup const *>(this)->getPointer<T>(lookup)); }

  template< typename T >
  T const * getPointer( char const * const name ) const
  { return getPointer<T>( string(name) ); }

  template< typename T >
  T * getPointer( char const * const name )
  { return getPointer<T>( string(name) ); }


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

  static axom::sidre::Group * setSidreGroup( string const& name,
                                             ManagedGroup * const parent );

  ManagedGroup * getParent()             { return m_parent; }
  ManagedGroup const * getParent() const { return m_parent; }

  ManagedGroup * setParent( ManagedGroup * const parent )
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


  viewWrapperMap const & wrappers() const
  {
    return m_wrappers;
  }

  viewWrapperMap & wrappers()
  {
    return m_wrappers;
  }

  RestartFlags getRestartFlags() const { return m_restart_flags; }

  void setRestartFlags( RestartFlags flags ) { m_restart_flags = flags; } 

  InputFlags getInputFlags() const { return m_input_flags; }

  void setInputFlags( InputFlags flags ) { m_input_flags = flags; }

  void prepareToWrite() const;

  void finishWriting() const;

  void prepareToRead();

  void finishReading();

protected:
  /**
   * @brief Post processing of the input values.
   */
  virtual void PostProcessInput() {}

  virtual void InitializePreSubGroups( ManagedGroup * const group ) {}

  virtual void InitializePostSubGroups( ManagedGroup * const group ) {}

  virtual void InitializePostInitialConditions_PreSubGroups( ManagedGroup * const group ) {}

  virtual void InitializePostInitialConditions_PostSubGroups( ManagedGroup * const group ) {}


private:  
  /**
   * @brief Read values from the input file and put them into the
   *        wrapped values for this group.
   * @param[in] targetNode the XML node that to extract input values from.
   */
  virtual void ProcessInputFile( xmlWrapper::xmlNode const & targetNode );

  /// the parent of this group
  ManagedGroup* m_parent = nullptr;

  /// the container for all wrappers
  viewWrapperMap m_wrappers;

  /// The container for all sub-groups
  subGroupMap m_subGroups;

#ifdef GEOSX_USE_ATK
  /// Pointer to the sidre group that mirrors this group
  axom::sidre::Group* m_sidreGroup;
#endif

  indexType m_size;             ///< The size/length wrappers in this group
  indexType m_capacity;         ///< The capacity for wrappers in this group
  RestartFlags m_restart_flags; ///< Restart flag for this group...and
                                ///< subsequently all wrappers in this group
  InputFlags m_input_flags;     ///< Input flag for this group
  string m_name;                ///< the repository name of this group. This
                                ///< is the key in the parent group.

};

using GroupKey = ManagedGroup::subGroupMap::KeyIndex;
using ViewKey = ManagedGroup::viewWrapperMap::KeyIndex;



template < typename T >
T * ManagedGroup::RegisterGroup( std::string const & name,
                                 std::unique_ptr<ManagedGroup> newObject )
{
#ifdef USE_DYNAMIC_CASTING
  return dynamic_cast<T*>( m_subGroups.insert( name, newObject.release(), true ) );
#else
  return static_cast<T*>( m_subGroups.insert( name, newObject.release(), true ) );
#endif
}


template< typename T >
T * ManagedGroup::RegisterGroup( std::string const & name,
                                 T * newObject,
                                 bool const takeOwnership )
{
#ifdef USE_DYNAMIC_CASTING
  return dynamic_cast<T*>( m_subGroups.insert( name, newObject, takeOwnership ) );
#else
  return static_cast<T*>( m_subGroups.insert( name, newObject, takeOwnership ) );
#endif
}


template< typename T, typename TBASE >
ViewWrapper<TBASE> * ManagedGroup::RegisterViewWrapper( std::string const & name,
                                                        ViewKey::index_type * const rkey )
{
  m_wrappers.insert( name,
                     (ViewWrapper<TBASE>::template Factory<T>(name,this) ).release(),
                     true );

  if( rkey != nullptr )
  {
    *rkey = m_wrappers.getIndex(name);
  }
  ViewWrapper<TBASE> * const rval = getWrapper<TBASE>(name);
  if( rval->sizedFromParent() == 1 && rval->shouldResize())
  {
    rval->resize(this->size());
  }
  return rval;
}

template< typename T, typename TBASE >
ViewWrapper<TBASE> * ManagedGroup::RegisterViewWrapper( ViewKey & viewKey )
{
  ViewKey::index_type index;
  ViewWrapper<TBASE> * const rval = RegisterViewWrapper<T,TBASE>( viewKey.Key(), &index );
  viewKey.setIndex(index);

  return rval;
}


template < typename T >
ViewWrapper<T> * ManagedGroup::RegisterViewWrapper( std::string const & name,
                                                    std::unique_ptr<T> newObject )
{
  m_wrappers.insert( name,
                     new ViewWrapper<T>( name, this, newObject.release(), true ),
                     true );

  ViewWrapper<T> * const rval = getWrapper<T>(name);
  if( rval->sizedFromParent() == 1 )
  {
    rval->resize(this->size());
  }
  return rval;
}



template< typename T >
ViewWrapper<T> * ManagedGroup::RegisterViewWrapper( std::string const & name,
                                                    T * newObject,
                                                    bool takeOwnership )
{
  m_wrappers.insert( name,
                     new ViewWrapper<T>( name, this, newObject, takeOwnership ),
                     true );

  ViewWrapper<T> * const rval = getWrapper<T>(name);
  if( rval->sizedFromParent() == 1 && rval->shouldResize())
  {
    rval->resize(this->size());
  }
  return rval;
}


} /* end namespace dataRepository */
} /* end namespace geosx */


//typedef geosx::dataRepository::ManagedGroup ObjectDataStructureBaseT;

#endif /* MANAGEDGROUP_H_ */
