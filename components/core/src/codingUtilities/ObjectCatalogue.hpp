/*
 * SingletonFactory.hpp
 *
 *  Created on: Nov 30, 2014
 *      Author: rrsettgast
 */

#ifndef OBJECTCATALOGUE_HPP_
#define OBJECTCATALOGUE_HPP_
#include<unordered_map>
#include<string>
#include<map>
#include<iostream>
#include<memory>
#include <cxxabi.h>
#include "macros.hpp"
#include "StringUtilities.hpp"

#ifndef OBJECTCATALOGVERBOSE
#define OBJECTCATALOGVERBOSE 2
#endif
#ifndef OBJECTCATALOGCONSTRUCTOR
#define OBJECTCATALOGCONSTRUCTOR 1
#endif

namespace objectcatalogue
{

template<typename BASETYPE, typename ...ARGS>
class CatalogInterface
{
public:
#if OBJECTCATALOGCONSTRUCTOR == 1
  typedef std::unordered_map<std::string, std::unique_ptr< CatalogInterface<BASETYPE, ARGS...> > > CatalogType;
#else
  typedef std::unordered_map<std::string, CatalogInterface<BASETYPE, ARGS...> * const  > CatalogType;
#endif

  CatalogInterface()
  {
#if OBJECTCATALOGVERBOSE > 1
    std::cout << "Calling constructor for CatalogueEntryBase< "<< geosx::stringutilities::demangle(typeid(BASETYPE).name())
                                                              <<" , ... >"<<std::endl;
#endif
  }

  virtual ~CatalogInterface()
  {
#if OBJECTCATALOGVERBOSE > 1
    std::cout << "Calling destructor for CatalogueEntryBase< "<< geosx::stringutilities::demangle(typeid(BASETYPE).name())
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
#if OBJECTCATALOGCONSTRUCTOR == 1
    CatalogInterface<BASETYPE, ARGS...> const * const entry = GetCatalogue().at( objectTypeName ).get();
#else
    CatalogInterface<BASETYPE, ARGS...> const * const entry = GetCatalogue().at( objectTypeName );
#endif
    return entry->Allocate( objectTypeName, args... );
  }

};

template<typename TYPE, typename BASETYPE, typename ...ARGS>
class CatalogueEntry : public CatalogInterface<BASETYPE, ARGS...>
{
public:
  CatalogueEntry() :
    CatalogInterface<BASETYPE, ARGS...>()
  {
#if OBJECTCATALOGVERBOSE > 1
    std::cout << "Calling constructor for CatalogueEntry< "<< geosx::stringutilities::demangle(typeid(TYPE).name())
                                                          <<" , "<<geosx::stringutilities::demangle(typeid(BASETYPE).name())
                                                          <<" , ... >"<<std::endl;
#endif

#if OBJECTCATALOGCONSTRUCTOR == 0
    std::string name = TYPE::CatalogueName();
    ( CatalogInterface<BASETYPE, ARGS...>::GetCatalogue() ).insert( std::make_pair( name, this ) );
#endif
  }

  ~CatalogueEntry() override final
  {
#if OBJECTCATALOGVERBOSE > 1
    std::cout << "Calling destructor for CatalogueEntry< "<< geosx::stringutilities::demangle(typeid(TYPE).name())
                                                          <<" , "<<geosx::stringutilities::demangle(typeid(BASETYPE).name())
                                                          <<" , ... >"<<std::endl;
#endif

  }

  CatalogueEntry(CatalogueEntry const & source ) :
    CatalogInterface<BASETYPE, ARGS...>(source)
  {}

  CatalogueEntry(CatalogueEntry && source ):
    CatalogInterface<BASETYPE, ARGS...>(std::move(source))
  {}

  CatalogueEntry& operator=(CatalogueEntry const & source )
  {
    CatalogInterface<BASETYPE, ARGS...>::operator=(source);
  }

  CatalogueEntry& operator=(CatalogueEntry && source )
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



#if OBJECTCATALOGCONSTRUCTOR == 1
template<typename TYPE, typename BASETYPE, typename ...ARGS>
class CatalogueEntryConstructor
{
public:
  CatalogueEntryConstructor()
  {
#if OBJECTCATALOGVERBOSE > 1
    std::cout << "Calling constructor for CatalogueEntryConstructor< "<< geosx::stringutilities::demangle(typeid(TYPE).name())
                                                                      <<" , "<<geosx::stringutilities::demangle(typeid(BASETYPE).name())
                                                                      <<" , ... >"<<std::endl;
#endif

    std::string name = TYPE::CatalogueName();
    ( CatalogInterface<BASETYPE, ARGS...>::GetCatalogue() ).insert( std::make_pair( name, std::make_unique< CatalogueEntry<TYPE,BASETYPE, ARGS...> >() ) );

#if OBJECTCATALOGVERBOSE > 0
    std::cout <<"Registered  "
              <<geosx::stringutilities::demangle(typeid(BASETYPE).name())
              <<" catalogue component of derived type "
              <<geosx::stringutilities::demangle(typeid(TYPE).name())
              <<" where "<<geosx::stringutilities::demangle(typeid(TYPE).name())<<"::CatalogueName() = "<<TYPE::CatalogueName() << std::endl;
#endif
  }

  ~CatalogueEntryConstructor()
  {
#if OBJECTCATALOGVERBOSE > 1
    std::cout << "Calling destructor for CatalogueEntryConstructor< "<< geosx::stringutilities::demangle(typeid(TYPE).name())
                                                                      <<" , "<<geosx::stringutilities::demangle(typeid(BASETYPE).name())
                                                                      <<" , ... >"<<std::endl;
#endif
  }

  CatalogueEntryConstructor(CatalogueEntryConstructor const &) = delete;
  CatalogueEntryConstructor(CatalogueEntryConstructor &&) = delete;
  CatalogueEntryConstructor& operator=(CatalogueEntryConstructor const &) = delete;
  CatalogueEntryConstructor& operator=(CatalogueEntryConstructor &&) = delete;

};
#endif

}


/// Compiler directive to simplify registration

#if OBJECTCATALOGCONSTRUCTOR == 1
#define REGISTER_CATALOGUE_ENTRY( BaseType, ClassName, ...) \
_Pragma("clang diagnostic push")\
namespace \
{ \
objectcatalogue::CatalogueEntryConstructor<ClassName,BaseType,__VA_ARGS__> catEntry; }
#else
#define REGISTER_CATALOGUE_ENTRY( BaseType, ClassName, ...) \
namespace \
{ \
objectcatalogue::CatalogueEntry<ClassName,BaseType,__VA_ARGS__> catEntry; }
#endif
//_Pragma("clang diagnostic ignored \"-Wglobal-constructors\"")\
//_Pragma("clang diagnostic ignored \"-Wexit-time-destructors\"")\
//_Pragma("clang diagnostic pop")



#endif /* OBJECTCATALOGUE_HPP_ */
