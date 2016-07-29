/**
 * @file   ObjectCatalog.hpp
 * @author Randolph Settgast
 */

#ifndef OBJECTCATALOG_HPP_
#define OBJECTCATALOG_HPP_
#include <unordered_map>
#include <string>
#include <iostream>
#include "StringUtilities.hpp"

#ifndef OBJECTCATALOGVERBOSE
#define OBJECTCATALOGVERBOSE 0
#endif

/**
 * namespace to hold the object catalog classes
 */
namespace objectcatalog
{

/**
 *  This class provides the base class/interface for the catalog value objects
 *  @tparam BASETYPE This is the base class of the objects that the factory produces.
 *  @tparam ...ARGS  variadic template pack to hold the parameters needed for the constructor of the BASETYPE
 */
template<typename BASETYPE, typename ... ARGS>
class CatalogInterface
{
public:
  /// This is the type that will be used for the catalog. The catalog is actually instantiated in the BASETYPE
  typedef std::unordered_map<std::string, std::unique_ptr< CatalogInterface<BASETYPE, ARGS ...> > > CatalogType;

  /// default constructor.
  CatalogInterface()
  {
#if OBJECTCATALOGVERBOSE > 1
    std::cout << "Calling constructor for CatalogInterface< "<< geosx::stringutilities::demangle(typeid(BASETYPE).name())
              <<" , ... >"<<std::endl;
#endif
  }

  ///default destructor
  virtual ~CatalogInterface()
  {
#if OBJECTCATALOGVERBOSE > 1
    std::cout << "Calling destructor for CatalogInterface< "<< geosx::stringutilities::demangle(typeid(BASETYPE).name())
              <<" , ... >"<<std::endl;
#endif
  }

  explicit CatalogInterface(CatalogInterface const &) = default;
  CatalogInterface(CatalogInterface &&) = default;
  CatalogInterface& operator=(CatalogInterface const &) = default;
  CatalogInterface& operator=(CatalogInterface &&) = default;

  /**
   * get the catalog from that is stored in the target base class.
   * @return returns the catalog for this
   */
  static CatalogType& GetCatalog()
  {
    return BASETYPE::GetCatalog();
  }

  /**
   * pure virtual to create a new object that derives from BASETYPE
   * @param name this is the key that was used to select the correct catalog entry
   * @param args these are the arguments to the constructor of the target type
   * @return passes a unique_ptr<BASETYPE> to the newly allocated class.
   */
  virtual std::unique_ptr<BASETYPE> Allocate( ARGS& ... args ) const = 0;

  /**
   * static method to create a new object that derives from BASETYPE
   * @param objectTypeName The key to the catalog entry that is able to create the correct type.
   * @param args these are the arguments to the constructor of the target type
   * @return passes a unique_ptr<BASETYPE> to the newly allocated class.
   */
  static std::unique_ptr<BASETYPE> Factory( std::string const & objectTypeName, ARGS& ... args )
  {
    CatalogInterface<BASETYPE, ARGS ...> const * const entry = GetCatalog().at( objectTypeName ).get();
    return entry->Allocate( args ... );
  }

};

/**
 * class to hold allocation capability for specific target derived types
 * @tparam TYPE this is the derived type
 * @tparam BASETYPE this is the base class that TYPE derives from
 * @tparam ARGS constructor arguemtns
 */
template<typename BASETYPE, typename TYPE, typename ... ARGS>
class CatalogEntry : public CatalogInterface<BASETYPE, ARGS ...>
{
public:
  /// default constructor
  CatalogEntry() :
    CatalogInterface<BASETYPE, ARGS ...>()
  {
#if OBJECTCATALOGVERBOSE > 1
    std::cout << "Calling constructor for CatalogEntry< "<< geosx::stringutilities::demangle(typeid(TYPE).name())
              <<" , "<<geosx::stringutilities::demangle(typeid(BASETYPE).name())
              <<" , ... >"<<std::endl;
#endif
  }

  /// default destructor
  ~CatalogEntry() override final
  {
#if OBJECTCATALOGVERBOSE > 1
    std::cout << "Calling destructor for CatalogEntry< "<< geosx::stringutilities::demangle(typeid(TYPE).name())
              <<" , "<<geosx::stringutilities::demangle(typeid(BASETYPE).name())
              <<" , ... >"<<std::endl;
#endif

  }

  CatalogEntry(CatalogEntry const & source ) :
    CatalogInterface<BASETYPE, ARGS ...>(source)
  {}

  CatalogEntry(CatalogEntry && source ) :
    CatalogInterface<BASETYPE, ARGS ...>(std::move(source))
  {}

  CatalogEntry& operator=(CatalogEntry const & source )
  {
    CatalogInterface<BASETYPE, ARGS ...>::operator=(source);
  }

  CatalogEntry& operator=(CatalogEntry && source )
  {
    CatalogInterface<BASETYPE, ARGS ...>::operator=(std::move(source));
  }

  /**
   * inherited virtual to create a new object that derives from BASETYPE
   * @param name this is the key that was used to select the correct catalog entry
   * @param args these are the arguments to the constructor of the target type
   * @return passes a unique_ptr<BASETYPE> to the newly allocated class.
   */
  virtual std::unique_ptr<BASETYPE> Allocate( ARGS& ... args ) const override final
  {
#if OBJECTCATALOGVERBOSE > 0
    std::cout << "Creating type "<< geosx::stringutilities::demangle(typeid(TYPE).name())
              <<" from catalog of "<<geosx::stringutilities::demangle(typeid(BASETYPE).name())<<std::endl;
#endif
    return std::unique_ptr<BASETYPE>( new TYPE( args ... ) );
  }
};


/**
 * a class to generate the catalog entry
 */
template<typename BASETYPE, typename TYPE, typename ... ARGS>
class CatalogEntryConstructor
{
public:
  /**
   * Constructor creates a catalog entry using the key defined by TYPE::CatalogName(), and value of CatalogEntry<TYPE,BASETYPE,ARGS...>.
   * After the constructor is executed, this object may be destroyed without consequence.
   */
  CatalogEntryConstructor()
  {
#if OBJECTCATALOGVERBOSE > 1
    std::cout << "Calling constructor for CatalogueEntryConstructor< "<< geosx::stringutilities::demangle(typeid(TYPE).name())
              <<" , "<<geosx::stringutilities::demangle(typeid(BASETYPE).name())
              <<" , ... >"<<std::endl;
#endif

    std::string name = TYPE::CatalogName();
    ( CatalogInterface<BASETYPE, ARGS ...>::GetCatalog() ).insert( std::make_pair( name, std::make_unique< CatalogEntry<BASETYPE, TYPE, ARGS ...> >() ) );

#if OBJECTCATALOGVERBOSE > 0
    std::cout <<"Registered  "
              <<geosx::stringutilities::demangle(typeid(BASETYPE).name())
              <<" catalogue component of derived type "
              <<geosx::stringutilities::demangle(typeid(TYPE).name())
              <<" where "<<geosx::stringutilities::demangle(typeid(TYPE).name())<<"::CatalogueName() = "<<TYPE::CatalogName() << std::endl;
#endif
  }

  /// default destuctor
  ~CatalogEntryConstructor()
  {
#if OBJECTCATALOGVERBOSE > 1
    std::cout << "Calling destructor for CatalogueEntryConstructor< "<< geosx::stringutilities::demangle(typeid(TYPE).name())
              <<" , "<<geosx::stringutilities::demangle(typeid(BASETYPE).name())
              <<" , ... >"<<std::endl;
#endif
  }

  CatalogEntryConstructor(CatalogEntryConstructor const &) = delete;
  CatalogEntryConstructor(CatalogEntryConstructor &&) = delete;
  CatalogEntryConstructor& operator=(CatalogEntryConstructor const &) = delete;
  CatalogEntryConstructor& operator=(CatalogEntryConstructor &&) = delete;

};

}


/// Compiler directive to simplify registration

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wexit-time-destructors"
#endif

/**
 * Macro that takes in the base class of the catalog, the derived class, and the argument types for the constructor of
 * the derived class/base class, and create an object of type CatalogEntryConstructor<ClassName,BaseType,__VA_ARGS__> in
 * an anonymous namespace. This should be called from the source file for the derived class, which will result in the
 * generation of a CatalogEntry<BaseType,ClassName,...> prior to main().
 */
#define REGISTER_CATALOG_ENTRY( BaseType, DerivedType, ...) \
  namespace { objectcatalog::CatalogEntryConstructor<BaseType,DerivedType,__VA_ARGS__> catEntry; }

#ifdef __clang__
#pragma clang diagnostic push
#endif


#endif /* OBJECTCATALOG_HPP_ */
