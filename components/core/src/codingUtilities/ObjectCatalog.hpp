

#ifndef OBJECTCATALOG_HPP_
#define OBJECTCATALOG_HPP_
#include<unordered_map>
#include<string>
#include<iostream>
#include "StringUtilities.hpp"

#ifndef OBJECTCATALOGVERBOSE
#define OBJECTCATALOGVERBOSE 2
#endif

namespace objectcatalog
{

template<typename BASETYPE, typename ...ARGS>
class CatalogInterface
{
public:
  typedef std::unordered_map<std::string, std::unique_ptr< CatalogInterface<BASETYPE, ARGS...> > > CatalogType;

  CatalogInterface()
  {
#if OBJECTCATALOGVERBOSE > 1
    std::cout << "Calling constructor for CatalogInterface< "<< geosx::stringutilities::demangle(typeid(BASETYPE).name())
                                                              <<" , ... >"<<std::endl;
#endif
  }

  virtual ~CatalogInterface()
  {
#if OBJECTCATALOGVERBOSE > 1
    std::cout << "Calling destructor for CatalogInterface< "<< geosx::stringutilities::demangle(typeid(BASETYPE).name())
                                                              <<" , ... >"<<std::endl;
#endif
  }

  CatalogInterface(CatalogInterface const &) = default;
  CatalogInterface(CatalogInterface &&) = default;
  CatalogInterface& operator=(CatalogInterface const &) = default;
  CatalogInterface& operator=(CatalogInterface &&) = default;

  static CatalogType& GetCatalogue()
  {
    return BASETYPE::GetCatalogue();
  }

  virtual std::unique_ptr<BASETYPE> Allocate( std::string const & name, ARGS&... args ) const = 0;

  static std::unique_ptr<BASETYPE> Factory( const std::string& objectTypeName, ARGS&...args )
  {
    CatalogInterface<BASETYPE, ARGS...> const * const entry = GetCatalogue().at( objectTypeName ).get();
    return entry->Allocate( objectTypeName, args... );
  }

};

template<typename TYPE, typename BASETYPE, typename ...ARGS>
class CatalogEntry : public CatalogInterface<BASETYPE, ARGS...>
{
public:
  CatalogEntry() :
    CatalogInterface<BASETYPE, ARGS...>()
  {
#if OBJECTCATALOGVERBOSE > 1
    std::cout << "Calling constructor for CatalogEntry< "<< geosx::stringutilities::demangle(typeid(TYPE).name())
                                                          <<" , "<<geosx::stringutilities::demangle(typeid(BASETYPE).name())
                                                          <<" , ... >"<<std::endl;
#endif
  }

  ~CatalogEntry() override final
  {
#if OBJECTCATALOGVERBOSE > 1
    std::cout << "Calling destructor for CatalogEntry< "<< geosx::stringutilities::demangle(typeid(TYPE).name())
                                                          <<" , "<<geosx::stringutilities::demangle(typeid(BASETYPE).name())
                                                          <<" , ... >"<<std::endl;
#endif

  }

  CatalogEntry(CatalogEntry const & source ) :
    CatalogInterface<BASETYPE, ARGS...>(source)
  {}

  CatalogEntry(CatalogEntry && source ):
    CatalogInterface<BASETYPE, ARGS...>(std::move(source))
  {}

  CatalogEntry& operator=(CatalogEntry const & source )
  {
    CatalogInterface<BASETYPE, ARGS...>::operator=(source);
  }

  CatalogEntry& operator=(CatalogEntry && source )
  {
    CatalogInterface<BASETYPE, ARGS...>::operator=(std::move(source));
  }

  virtual std::unique_ptr<BASETYPE> Allocate(  std::string const & name, ARGS&... args ) const override final
  {
#if OBJECTCATALOGVERBOSE > 0
    std::cout << "Creating "<< name <<" of type "<< geosx::stringutilities::demangle(typeid(TYPE).name())
              <<" from catalog of "<<geosx::stringutilities::demangle(typeid(BASETYPE).name())<<std::endl;
#endif
    return std::unique_ptr<BASETYPE>( new TYPE( args... ) );
  }
};



template<typename TYPE, typename BASETYPE, typename ...ARGS>
class CatalogEntryConstructor
{
public:
  CatalogEntryConstructor()
  {
#if OBJECTCATALOGVERBOSE > 1
    std::cout << "Calling constructor for CatalogueEntryConstructor< "<< geosx::stringutilities::demangle(typeid(TYPE).name())
                                                                      <<" , "<<geosx::stringutilities::demangle(typeid(BASETYPE).name())
                                                                      <<" , ... >"<<std::endl;
#endif

    std::string name = TYPE::CatalogName();
    ( CatalogInterface<BASETYPE, ARGS...>::GetCatalogue() ).insert( std::make_pair( name, std::make_unique< CatalogEntry<TYPE,BASETYPE, ARGS...> >() ) );

#if OBJECTCATALOGVERBOSE > 0
    std::cout <<"Registered  "
              <<geosx::stringutilities::demangle(typeid(BASETYPE).name())
              <<" catalogue component of derived type "
              <<geosx::stringutilities::demangle(typeid(TYPE).name())
              <<" where "<<geosx::stringutilities::demangle(typeid(TYPE).name())<<"::CatalogueName() = "<<TYPE::CatalogName() << std::endl;
#endif
  }

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

#define REGISTER_CATALOG_ENTRY( BaseType, ClassName, ...) \
namespace { objectcatalog::CatalogEntryConstructor<ClassName,BaseType,__VA_ARGS__> catEntry; }

#ifdef __clang__
#pragma clang diagnostic push
#endif


#endif /* OBJECTCATALOG_HPP_ */
