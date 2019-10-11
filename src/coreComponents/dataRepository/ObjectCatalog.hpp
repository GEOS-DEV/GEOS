#ifndef GEOSX_DATAREPOSITORY_OBJECTCATALOG_HPP_
#define GEOSX_DATAREPOSITORY_OBJECTCATALOG_HPP_
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
 * @file ObjectCatalog.hpp
 * The ObjectCatalog acts as a statically initialized factory. It functions in
 * a similar manner to classic virtual factory method, except that it is no
 * maintained list of derived objects that is required to create new objects.
 * Instead, the ``ObjectCatalog`` creates a "catalog" of derived objects using
 * a ``std::unordered_map``.
 * This ``std::unordered_map`` is then statically initialized through the declaration
 * of a
 */

#include "Logger.hpp"
#include "StringUtilities.hpp"

#include <unordered_map>
#include <string>
#include <iostream>
#include <memory>

#ifndef OBJECTCATALOGVERBOSE
#define OBJECTCATALOGVERBOSE 0
#endif


#ifndef BASEHOLDSCATALOG
#define BASEHOLDSCATALOG 1
#endif

namespace geosx
{
namespace dataRepository
{

/**
 *  This class provides the base class/interface for the catalog value objects
 *  @tparam BASETYPE This is the base class of the objects that the factory
 * produces.
 *  @tparam ARGS  variadic template pack to hold the parameters needed for the
 * constructor of the BASETYPE
 */
//START_SPHINX_0
template< typename BASETYPE, typename ... ARGS >
class CatalogInterface
{
public:
  /// This is the type that will be used for the catalog. The catalog is
  /// actually instantiated in the BASETYPE
  //START_SPHINX_1
  typedef std::unordered_map< std::string,
                              std::unique_ptr< CatalogInterface< BASETYPE, ARGS... > > > CatalogType;
  //STOP_SPHINX

  /// default constructor.
  CatalogInterface()
  {
#if OBJECTCATALOGVERBOSE > 1
    GEOS_LOG_RANK( "Calling constructor for CatalogInterface< " << cxx_utilities::demangle( typeid(BASETYPE).name())
                                                                << " , ... >" );
#endif
  }

  ///default destructor
  virtual ~CatalogInterface()
  {
#if OBJECTCATALOGVERBOSE > 1
    GEOS_LOG_RANK( "Calling destructor for CatalogInterface< "<< cxx_utilities::demangle( typeid(BASETYPE).name())
                                                              <<" , ... >" );
#endif
  }

  explicit CatalogInterface( CatalogInterface const & ) = default;
  CatalogInterface( CatalogInterface && ) = default;
  CatalogInterface & operator=( CatalogInterface const & ) = default;
  CatalogInterface & operator=( CatalogInterface && ) = default;

  /**
   * get the catalog from that is stored in the target base class.
   * @return returns the catalog for this
   */
  static CatalogType & GetCatalog()
  {
#if BASEHOLDSCATALOG == 1
    return BASETYPE::GetCatalog();
#else
    static CatalogType catalog;
    return catalog;
#endif
  }

  /**
   * pure virtual to create a new object that derives from BASETYPE
   * @param args these are the arguments to the constructor of the target type
   * @return passes a unique_ptr<BASETYPE> to the newly allocated class.
   */
  virtual std::unique_ptr< BASETYPE > Allocate( ARGS... args ) const = 0;


  static bool hasKeyName( std::string const & objectTypeName )
  {
    return GetCatalog().count( objectTypeName );
  }

  /**
   * static method to create a new object that derives from BASETYPE
   * @param objectTypeName The key to the catalog entry that is able to create
   * the correct type.
   * @param args these are the arguments to the constructor of the target type
   * @return passes a unique_ptr<BASETYPE> to the newly allocated class.
   */
  //START_SPHINX_2
  static std::unique_ptr< BASETYPE > Factory( std::string const & objectTypeName, ARGS... args )
  {
    return GetCatalog().at( objectTypeName ).get()->Allocate( args ... );
  }
  //STOP_SPHINX

  template< typename TYPE >
  static TYPE & catalog_cast( BASETYPE & object )
  {
    std::string castedName = TYPE::CatalogName();
    std::string objectName = object.getName();

    if( castedName != objectName )
    {
#if OBJECTCATALOGVERBOSE > 1
      GEOS_LOG_RANK( "Invalid Cast of " << objectName << " to " << castedName );
#endif
    }

    return static_cast< TYPE & >(object);
  }

};

/**
 * class to hold allocation capability for specific target derived types
 * @tparam TYPE this is the derived type
 * @tparam BASETYPE this is the base class that TYPE derives from
 * @tparam ARGS constructor arguments
 */
//START_SPHINX_3
template< typename BASETYPE, typename TYPE, typename ... ARGS >
class CatalogEntry : public CatalogInterface< BASETYPE, ARGS... >
{
public:
  /// default constructor
  CatalogEntry():
    CatalogInterface< BASETYPE, ARGS... >()
  {
#if OBJECTCATALOGVERBOSE > 1
    GEOS_LOG_RANK( "Calling constructor for CatalogEntry< " << cxx_utilities::demangle( typeid(TYPE).name())
                                                            << " , " << cxx_utilities::demangle( typeid(BASETYPE).name())
                                                            << " , ... >" );
#endif
  }

  /// default destructor
  ~CatalogEntry() override final
  {
#if OBJECTCATALOGVERBOSE > 1
    GEOS_LOG_RANK( "Calling destructor for CatalogEntry< " << cxx_utilities::demangle( typeid(TYPE).name())
                                                           << " , " << cxx_utilities::demangle( typeid(BASETYPE).name())
                                                           << " , ... >" );
#endif

  }

  CatalogEntry( CatalogEntry const & source ):
    CatalogInterface< BASETYPE, ARGS... >( source )
  {}

  CatalogEntry( CatalogEntry && source ):
    CatalogInterface< BASETYPE, ARGS... >( std::move( source ))
  {}

  CatalogEntry & operator=( CatalogEntry const & source )
  {
    CatalogInterface< BASETYPE, ARGS... >::operator=( source );
  }

  CatalogEntry & operator=( CatalogEntry && source )
  {
    CatalogInterface< BASETYPE, ARGS... >::operator=( std::move(source));
  }

  /**
   * inherited virtual to create a new object that derives from BASETYPE
   * @param args these are the arguments to the constructor of the target type
   * @return passes a unique_ptr<BASETYPE> to the newly allocated class.
   */
  //START_SPHINX_4
  virtual std::unique_ptr< BASETYPE > Allocate( ARGS... args ) const override final
  {
#if OBJECTCATALOGVERBOSE > 0
    GEOS_LOG_RANK( "Creating type " << cxx_utilities::demangle( typeid(TYPE).name())
                                    << " from catalog of " << cxx_utilities::demangle( typeid(BASETYPE).name()));
#endif
#if ( __cplusplus >= 201402L )
    return std::make_unique< TYPE >( args ... );
#else
    return std::unique_ptr< BASETYPE >( new TYPE( args ... ) );
#endif
  }
  //STOP_SPHINX

};


/**
 * a class to generate the catalog entry
 */
template< typename BASETYPE, typename TYPE, typename ... ARGS >
class CatalogEntryConstructor
{
public:
  /**
   * Constructor creates a catalog entry using the key defined by
   * TYPE::CatalogName(), and value of CatalogEntry<TYPE,BASETYPE,ARGS...>.
   * After the constructor is executed, this object may be destroyed without
   * consequence.
   */
  CatalogEntryConstructor()
  {
#if OBJECTCATALOGVERBOSE > 1
    GEOS_LOG_RANK( "Calling constructor for CatalogEntryConstructor< " << cxx_utilities::demangle( typeid(TYPE).name())
                                                                       << " , " << cxx_utilities::demangle( typeid(BASETYPE).name())
                                                                       << " , ... >" );
#endif

    std::string name = TYPE::CatalogName();
#if ( __cplusplus >= 201402L )
    std::unique_ptr< CatalogEntry< BASETYPE, TYPE, ARGS... > > temp = std::make_unique< CatalogEntry< BASETYPE, TYPE, ARGS... > >();
#else
    std::unique_ptr< CatalogEntry< BASETYPE, TYPE, ARGS... > > temp = std::unique_ptr< CatalogEntry< BASETYPE, TYPE, ARGS... > >( new CatalogEntry< BASETYPE,
                                                                                                                                                    TYPE,
                                                                                                                                                    ARGS... >()  );
#endif
    ( CatalogInterface< BASETYPE, ARGS... >::GetCatalog() ).insert( std::move( std::make_pair( name, std::move( temp ) ) ) );

#if OBJECTCATALOGVERBOSE > 0
    GEOS_LOG_RANK( "Registered " << cxx_utilities::demangle( typeid(BASETYPE).name())
                                 << " catalog component of derived type "
                                 << cxx_utilities::demangle( typeid(TYPE).name())
                                 << " where " << cxx_utilities::demangle( typeid(TYPE).name())
                                 << "::CatalogName() = " << TYPE::CatalogName());
#endif
  }

  /// default destuctor
  ~CatalogEntryConstructor()
  {
#if OBJECTCATALOGVERBOSE > 1
    GEOS_LOG_RANK( "Calling destructor for CatalogEntryConstructor< " << cxx_utilities::demangle( typeid(TYPE).name())
                                                                      << " , " << cxx_utilities::demangle( typeid(BASETYPE).name())
                                                                      << " , ... >" );
#endif
  }

  CatalogEntryConstructor( CatalogEntryConstructor const & ) = delete;
  CatalogEntryConstructor( CatalogEntryConstructor && ) = delete;
  CatalogEntryConstructor & operator=( CatalogEntryConstructor const & ) = delete;
  CatalogEntryConstructor & operator=( CatalogEntryConstructor && ) = delete;

};

/// Specialization for constructors with empty argument list
template< typename BASETYPE >
class CatalogInterface< BASETYPE >
{
public:
  /// This is the type that will be used for the catalog. The catalog is
  // actually instantiated in the BASETYPE
  typedef std::unordered_map< std::string, std::unique_ptr< CatalogInterface< BASETYPE > > > CatalogType;

  /// default constructor.
  CatalogInterface()
  {
#if OBJECTCATALOGVERBOSE > 1
    GEOS_LOG_RANK( "Calling constructor for CatalogInterface< " << cxx_utilities::demangle( typeid(BASETYPE).name())
                                                                << " , ... >" );
#endif
  }

  ///default destructor
  virtual ~CatalogInterface()
  {
#if OBJECTCATALOGVERBOSE > 1
    GEOS_LOG_RANK( "Calling destructor for CatalogInterface< " << cxx_utilities::demangle( typeid(BASETYPE).name())
                                                               << " , ... >" );
#endif
  }

  explicit CatalogInterface( CatalogInterface const & ) = default;
  CatalogInterface( CatalogInterface && ) = default;
  CatalogInterface & operator=( CatalogInterface const & ) = default;
  CatalogInterface & operator=( CatalogInterface && ) = default;

  /**
   * get the catalog from that is stored in the target base class.
   * @return returns the catalog for this
   */
  static CatalogType & GetCatalog()
  {
#if BASEHOLDSCATALOG == 1
    return BASETYPE::GetCatalog();
#else
    static CatalogType catalog;
    return catalog;
#endif
  }

  /**
   * pure virtual to create a new object that derives from BASETYPE
   * @param args these are the arguments to the constructor of the target type
   * @return passes a unique_ptr<BASETYPE> to the newly allocated class.
   */
  virtual std::unique_ptr< BASETYPE > Allocate(  ) const = 0;

  /**
   * static method to create a new object that derives from BASETYPE
   * @param objectTypeName The key to the catalog entry that is able to create
   * the correct type.
   * @param args these are the arguments to the constructor of the target type
   * @return passes a unique_ptr<BASETYPE> to the newly allocated class.
   */
  static std::unique_ptr< BASETYPE > Factory( std::string const & objectTypeName )
  {
    CatalogInterface< BASETYPE > const * const entry = GetCatalog().at( objectTypeName ).get();
    return entry->Allocate();
  }

  template< typename TYPE >
  static TYPE & catalog_cast( BASETYPE & object )
  {
    std::string castedName = TYPE::CatalogName();
    std::string objectName = object.getName();

    if( castedName != objectName )
    {
#if OBJECTCATALOGVERBOSE > 1
      GEOS_LOG_RANK( "Invalid Cast of " << objectName << " to " << castedName );
#endif
    }

    return static_cast< TYPE & >(object);
  }

};

template< typename BASETYPE, typename TYPE >
class CatalogEntry< BASETYPE, TYPE > : public CatalogInterface< BASETYPE >
{
public:
  /// default constructor
  CatalogEntry():
    CatalogInterface< BASETYPE >()
  {
#if OBJECTCATALOGVERBOSE > 1
    GEOS_LOG_RANK( "Calling constructor for CatalogEntry< " << cxx_utilities::demangle( typeid(TYPE).name())
                                                            << " , " << cxx_utilities::demangle( typeid(BASETYPE).name())
                                                            << " , ... >" );
#endif
  }

  /// default destructor
  ~CatalogEntry() override final
  {
#if OBJECTCATALOGVERBOSE > 1
    GEOS_LOG_RANK( "Calling destructor for CatalogEntry< " << cxx_utilities::demangle( typeid(TYPE).name())
                                                           << " , " << cxx_utilities::demangle( typeid(BASETYPE).name())
                                                           << " , ... >" );
#endif

  }

  CatalogEntry( CatalogEntry const & source ):
    CatalogInterface< BASETYPE >( source )
  {}

  CatalogEntry( CatalogEntry && source ):
    CatalogInterface< BASETYPE >( std::move( source ))
  {}

  CatalogEntry & operator=( CatalogEntry const & source )
  {
    CatalogInterface< BASETYPE >::operator=( source );
  }

  CatalogEntry & operator=( CatalogEntry && source )
  {
    CatalogInterface< BASETYPE >::operator=( std::move(source));
  }

  virtual std::unique_ptr< BASETYPE > Allocate(  ) const override final
  {
#if OBJECTCATALOGVERBOSE > 0
    GEOS_LOG_RANK( "Creating type " << cxx_utilities::demangle( typeid(TYPE).name())
                                    << " from catalog of " << cxx_utilities::demangle( typeid(BASETYPE).name()));
#endif
#if ( __cplusplus >= 201402L )
    return std::make_unique< TYPE >(  );
#else
    return std::unique_ptr< BASETYPE >( new TYPE(  ) );
#endif
  }
};



template< typename BASETYPE, typename TYPE >
class CatalogEntryConstructor< BASETYPE, TYPE >
{
public:
  CatalogEntryConstructor()
  {
#if OBJECTCATALOGVERBOSE > 1
    GEOS_LOG_RANK( "Calling constructor for CatalogEntryConstructor< " << cxx_utilities::demangle( typeid(TYPE).name())
                                                                       << " , " << cxx_utilities::demangle( typeid(BASETYPE).name())
                                                                       << " , ... >" );
#endif

    std::string name = TYPE::CatalogName();
#if ( __cplusplus >= 201402L )
    std::unique_ptr< CatalogEntry< BASETYPE, TYPE > > temp = std::make_unique< CatalogEntry< BASETYPE, TYPE > >();
#else
    std::unique_ptr< CatalogEntry< BASETYPE, TYPE > > temp = std::unique_ptr< CatalogEntry< BASETYPE, TYPE > >( new CatalogEntry< BASETYPE, TYPE >()  );
#endif
    ( CatalogInterface< BASETYPE >::GetCatalog() ).insert( std::move( std::make_pair( name, std::move( temp ) ) ) );

#if OBJECTCATALOGVERBOSE > 0
    GEOS_LOG_RANK( "Registered " << cxx_utilities::demangle( typeid(BASETYPE).name())
                                 << " catalog component of derived type "
                                 << cxx_utilities::demangle( typeid(TYPE).name())
                                 << " where " << cxx_utilities::demangle( typeid(TYPE).name())
                                 << "::CatalogName() = " << TYPE::CatalogName());
#endif
  }

  /// default destuctor
  ~CatalogEntryConstructor()
  {
#if OBJECTCATALOGVERBOSE > 1
    GEOS_LOG_RANK( "Calling destructor for CatalogEntryConstructor< " << cxx_utilities::demangle( typeid(TYPE).name())
                                                                      << " , " << cxx_utilities::demangle( typeid(BASETYPE).name()) << " , ... >" );
#endif
  }

  CatalogEntryConstructor( CatalogEntryConstructor const & ) = delete;
  CatalogEntryConstructor( CatalogEntryConstructor && ) = delete;
  CatalogEntryConstructor & operator=( CatalogEntryConstructor const & ) = delete;
  CatalogEntryConstructor & operator=( CatalogEntryConstructor && ) = delete;

};


}
}


/**
 * Macro that takes in the base class of the catalog, the derived class, and the
 * argument types for the constructor of
 * the derived class/base class, and create an object of type
 * CatalogEntryConstructor<ClassName,BaseType,__VA_ARGS__> in
 * an anonymous namespace. This should be called from the source file for the
 * derived class, which will result in the
 * generation of a CatalogEntry<BaseType,ClassName,...> prior to main().
 */
#define REGISTER_CATALOG_ENTRY( BaseType, DerivedType, ... ) \
  namespace { geosx::dataRepository::CatalogEntryConstructor< BaseType, DerivedType, __VA_ARGS__ > catEntry_ ## DerivedType; }

#define REGISTER_CATALOG_ENTRY0( BaseType, DerivedType ) \
  namespace { geosx::dataRepository::CatalogEntryConstructor< BaseType, DerivedType > catEntry_ ## DerivedType; }

#endif /* GEOSX_DATAREPOSITORY_OBJECTCATALOG_HPP_ */
