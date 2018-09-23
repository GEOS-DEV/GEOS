/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
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

/**
 * @file DataObjectManager.h
 * @date created on Nov 21, 2014
 * @author Randolph R. Settgast
 */


#ifndef MANAGEDGROUP_H_
#define MANAGEDGROUP_H_

#include <iostream>
#include <mpi.h>


#include "ObjectCatalog.hpp"
#include "ViewWrapper.hpp"
#include "RestartFlags.hpp"

#include "depricated/Common.h"
#include "DocumentationNode.hpp"

#include "MappedVector.hpp"

#include "fileIO/xmlWrapper.hpp"
//#include "CodingUtilities/ANSTexception.hpp"

#ifndef USE_DYNAMIC_CASTING
#define USE_DYNAMIC_CASTING 1;
#endif

#ifndef NOCHARTOSTRING_KEYLOOKUP
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

using keyType = string;
using indexType = localIndex;
//using DataKey = DataKeyT<keyType,indexType>;

/**
 * @author Randolph R. Settgast
 *
 * class that encapsulates and manages a collection of DataObjects. Can be
 * considered a "node" in a hierarchy of managers that represent physical groupings of data/
 *
 */
class ManagedGroup
{
public:
  using subGroupMap = MappedVector< ManagedGroup, ManagedGroup*, keyType, indexType  >;
  using viewWrapperMap = MappedVector< ViewWrapperBase, ViewWrapperBase*, keyType, indexType  >;
  /**
   * @name constructors, destructor, copy, move, assignments
   */
  ///@{

  /**
   * @author Randolph R. Settgast
   * @param name the name of this object manager
   */
  explicit ManagedGroup( std::string const & name,
                         ManagedGroup * const parent );

//  explicit ManagedGroup( std::string const & name,
//                         ManagedGroup * const parent,
//                         cxx_utilities::DocumentationNode * docNode );


  /**
   * @brief move constructor
   * @param source source ManagedGroup
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



  virtual const std::type_info& get_typeid() const
  {
    return typeid(*this);
  }

  bool CheckTypeID( std::type_info const & typeToCheck ) const
  {
    return typeToCheck == get_typeid() ? true : false;
  }


  template< typename T = ManagedGroup >
  T * RegisterGroup( std::string const & name, std::unique_ptr<ManagedGroup> newObject );

  template< typename T = ManagedGroup >
  T * RegisterGroup( std::string const & name,
                     T * newObject,
                     bool const takeOwnership );

  template< typename T = ManagedGroup >
  T * RegisterGroup( std::string const & name )
  {
    return RegisterGroup<T>( name, std::move(std::make_unique< T >( name, this )) );
  }

  template< typename T = ManagedGroup >
  T * RegisterGroup( subGroupMap::KeyIndex & keyIndex )
  {
    T * rval = RegisterGroup<T>( keyIndex.Key(), std::move(std::make_unique< T >( keyIndex.Key(), this )) );
    keyIndex.setIndex( this->m_subGroups.getIndex(keyIndex.Key()) );
    return rval;
  }

  template< typename T = ManagedGroup, typename TBASE = ManagedGroup >
  T * RegisterGroup( std::string const & name, std::string const & catalogName )
  {
    std::unique_ptr<TBASE> newGroup = TBASE::CatalogInterface::Factory(catalogName, name, this );
    return RegisterGroup<T,TBASE>( name, std::move(newGroup) );
  }


  template< typename T >
  static T group_cast( ManagedGroup * group )
  {
#ifdef USE_DYNAMIC_CASTING
    return dynamic_cast<T>( group );
#else
    return static_cast<T>( group );
#endif
  }

  template< typename T >
  static T group_cast( ManagedGroup const * group )
  {
#ifdef USE_DYNAMIC_CASTING
    return dynamic_cast<T>( group );
#else
    return static_cast<T>( group );
#endif
  }

  template< typename T >
  T group_cast()
  {
#ifdef USE_DYNAMIC_CASTING
    return dynamic_cast<T>( this );
#else
    return static_cast<T>( this );
#endif
  }

  template< typename T >
  T group_cast() const
  {
#ifdef USE_DYNAMIC_CASTING
    return dynamic_cast<T>( this );
#else
    return static_cast<T>( this );
#endif
  }


  template< typename T = ManagedGroup >
  T * GetGroup( localIndex index )
  {
    return group_cast<T*>(m_subGroups[index]);
  }

  template< typename T = ManagedGroup >
  T const * GetGroup( localIndex index ) const
  {
    return group_cast<T const *>(m_subGroups[index]);
  }

  template< typename T = ManagedGroup >
  T * GetGroup( string const & name )
  {
    return group_cast<T *>(m_subGroups[name]);
  }

  template< typename T = ManagedGroup >
  T const * GetGroup( string const & name ) const
  {
    return group_cast<T const *>(m_subGroups[name]);
  }


  template< typename T = ManagedGroup >
  T * GetGroup( subGroupMap::KeyIndex & key )
  {
    return group_cast<T *>(m_subGroups[key]);
  }

  template< typename T = ManagedGroup >
  T const * GetGroup( subGroupMap::KeyIndex & key ) const
  {
    return group_cast<T const *>(m_subGroups[key]);
  }

  template< typename T = ManagedGroup >
  T * GetGroup( subGroupMap::KeyIndex const & key )
  {
    return group_cast<T *>(m_subGroups[key]);
  }

  template< typename T = ManagedGroup >
  T const * GetGroup( subGroupMap::KeyIndex const & key ) const
  {
    return group_cast<T const *>(m_subGroups[key]);
  }

  /**
   * @brief This will grab the pointer to an object in the data structure
   * @param path a unix-style string (absolute, relative paths valid)
   */
  template< typename T = ManagedGroup >
  T const * GetGroupByPath( string const & path ) const
  {
    size_t directoryMarker = path.find("/");

    if (directoryMarker == std::string::npos)
    {
      // Target should be a child of this group
      return m_subGroups[path];
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
          return this->GetGroupByPath(subPath);
        }
        else                               // Not at root
        {
          return this->getParent()->GetGroupByPath(path);
        }
      }
      else if (child[0] == '.')
      {
        if (child[1] == '.')               // '../' = Reverse path
        {
          return this->getParent()->GetGroupByPath(subPath); 
        }
        else                               // './' = This path
        {
          return this->GetGroupByPath(subPath);
        }
      }
      else
      {
        return m_subGroups[child]->GetGroupByPath(subPath);
      }
    }
  }
  
  template< typename T = ManagedGroup >
  T * GetGroupByPath( string const & path )
  {
    return const_cast<T *>(const_cast< ManagedGroup const * >(this)->GetGroupByPath(path));
  }

  subGroupMap & GetSubGroups()
  {
    return m_subGroups;
  }

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
      if ( wrapperIter.second->get_typeid() == typeid(Wrapped) )
      {
        auto & wrapper = ViewWrapper<Wrapped>::cast(*wrapperIter.second);
        lambda(wrapper);
      }
    }
  }

  template< typename Wrapped, typename LAMBDA >
  void forViewWrappersByType(LAMBDA lambda) const
  {
    for( auto const & wrapperIter : m_wrappers )
    {
      if( wrapperIter.second->get_typeid() == typeid(Wrapped) )
      {
        auto const & wrapper = ViewWrapper<Wrapped>::cast(*wrapperIter.second);
        lambda(wrapper);
      }
    }
  }

  virtual void Initialize( ManagedGroup * const group );

  virtual void InitializationOrder( string_array & order );

  virtual void InitializePreSubGroups( ManagedGroup * const group ) {}

  virtual void InitializePostSubGroups( ManagedGroup * const group ) {}

  virtual void FinalInitializationRecursive( ManagedGroup * const group );
  virtual void FinalInitialization( ManagedGroup * const group ){}


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
   * @param name the key name to use for this new wrapper
   * @param wrapper a pointer to the new wrapper
   * @return a ViewWrapperBase pointer that holds the address of the new wrapper
   */
  ViewWrapperBase * RegisterViewWrapper( string const & name,
                                         ViewWrapperBase * const wrapper );

//  template< typename T >
//  void RegisterViewWrapperRecursive( string const & name );
//
//  template< typename T >
//  void RegisterViewWrapperRecursive( string const & name, string const & targetGroupName );

  ///@}


  /**
   * @name Self Documentation Functions
   */
  ///@{

  cxx_utilities::DocumentationNode * getDocumentationNode()
  {
    return m_docNode;
  }

  void RegisterDocumentationNodes();


  ///@}
  ///


  void PrintDataHierarchy(integer indent = 0);

  virtual void AddChildren( xmlWrapper::xmlNode const & targetNode );

  virtual void CreateChild( string const & childKey, string const & childName );

  virtual void ReadXML( xmlWrapper::xmlNode const & targetNode );

  virtual void ReadXMLsub( xmlWrapper::xmlNode const & targetNode );

  /**
   * This function provides a mechanism by which to post process any values that were read into the
   * xml file prior to initialization.
   */
  virtual void ReadXML_PostProcess() {}

  virtual void BuildDataStructure( dataRepository::ManagedGroup * const rootGroup );

  void SetDocumentationNodes();

  /**
   * Function to generate documentation nodes for each variable in this object. The documentation
   * nodes are then used to register variables and read xml input into variables.
   */
  virtual void FillDocumentationNode();


  /**
   * @param rootGroup The group for which to register new documentation node to.
   * Function to generate documentation nodes for each variable in this an object. The documentation
   * nodes are then used to register variables and read xml input into variables.
   */
  void SetOtherDocumentationNodes(dataRepository::ManagedGroup * const rootGroup);

  virtual void FillOtherDocumentationNodes( dataRepository::ManagedGroup * const group );
  
  virtual localIndex PackSize( string_array const & wrapperNames,
                        integer const recursive ) const;

  virtual localIndex PackSize( string_array const & wrapperNames,
                        localIndex_array const & packList,
                        integer const recursive ) const;

  virtual localIndex Pack( buffer_unit_type * & buffer,
                    string_array const & wrapperNames,
                    integer const recursive ) const;

  virtual localIndex Pack( buffer_unit_type * & buffer,
                    string_array const & wrapperNames,
                    localIndex_array const & packList,
                    integer const recursive ) const;

  virtual localIndex Unpack( buffer_unit_type const *& buffer,
                      localIndex_array & packList,
                      integer const recursive );


  //***********************************************************************************************

  // user defined conversion doesn't work. can't infer template argument
//  class GetDataClass
//  {
//  public:
//    GetDataClass( ManagedGroup & parent ): m_parent( parent ) {}
//
//    inline GetDataClass& operator() ( std::string const & name )
//    {
//      m_name = name;
//      return *this;
//    }
//
//    template< typename T>
//    operator typename ViewWrapper<T>::rtype ()
//    {
//      return m_parent.getData<T>( m_name );
//    }
//  private:
//    ManagedGroup & m_parent;
//    std::string m_name;
//  };
//  GetDataClass GetData = {*this};


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






  template< typename T, typename LOOKUP_TYPE >
  view_rtype_const<T> getData( LOOKUP_TYPE const index ) const
  { return getWrapper<T>(index)->data(); }

  template< typename T, typename LOOKUP_TYPE >
  view_rtype<T> getData( LOOKUP_TYPE const index )
  { return getWrapper<T>(index)->data(); }

  template< typename T >
  view_rtype_const<T> getData( char const * const name ) const
  { return getWrapper<T>( string(name) )->data(); }

  template< typename T >
  view_rtype<T> getData( char const * const name )
  { return getWrapper<T>( string(name) )->data(); }


  template< typename T, typename LOOKUP_TYPE >
  T const & getReference( LOOKUP_TYPE const & lookup ) const
  {
    ViewWrapper<T> const * wrapper = getWrapper<T>(lookup);
    if( wrapper == nullptr )
    {
      GEOS_ERROR( "ManagedGroup::getReferenceT(): call to getWrapper results in nullptr" );
    }
    return wrapper->reference();
  }

  template< typename T, typename LOOKUP_TYPE >
  T & getReference( LOOKUP_TYPE const & lookup )
  { return const_cast<T&>( const_cast<const ManagedGroup*>(this)->getReference<T>( lookup ) ); }

  template< typename T >
  T const & getReference( char const * const name ) const
  { return getReference<T>( string(name) ); }

  template< typename T >
  T & getReference( char const * const name )
  { return const_cast<T&>( const_cast<const ManagedGroup*>(this)->getReference<T>( name ) ); }



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

  bool hasView( std::string const & name ) const
  {
    return (m_wrappers[name] != nullptr);
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

  void prepareToWrite() const;

  void finishWriting() const;

  void prepareToRead();

  void finishReading();


protected:
  cxx_utilities::DocumentationNode * m_docNode = nullptr;

private:  

  ManagedGroup* m_parent = nullptr;
  viewWrapperMap m_wrappers;
  subGroupMap m_subGroups;

#ifdef GEOSX_USE_ATK
  axom::sidre::Group* m_sidreGroup;
#endif

  indexType m_size;
  indexType m_capacity;
  RestartFlags m_restart_flags;
  string m_name;

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

//template< typename T >
//void ManagedGroup::RegisterViewWrapperRecursive( string const & name )
//{
//  this->RegisterViewWrapper<T>(name);
//  forSubGroups( [&] ( ManagedGroup & group ) -> void
//  {
//    group.RegisterViewWrapperRecursive<T>(name);
//  });
//}
//
//template< typename T >
//void ManagedGroup::RegisterViewWrapperRecursive( string const & name, string const & targetGroupName )
//{
//  if( this->m_name == targetGroupName )
//  {
//    forSubGroups( [&] ( ManagedGroup & group ) -> void
//    {
//      this->RegisterViewWrapperRecursive<T>(name);
//    });
//  }
//  else
//  {
//    forSubGroups( [&] ( ManagedGroup & group ) -> void
//    {
//      group.RegisterViewWrapperRecursive<T>(name, targetGroupName);
//    });
//  }
//}






} /* end namespace dataRepository */
} /* end namespace geosx */


//typedef geosx::dataRepository::ManagedGroup ObjectDataStructureBaseT;

#endif /* MANAGEDGROUP_H_ */
