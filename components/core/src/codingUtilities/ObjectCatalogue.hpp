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

namespace objectcatalogue
{

template<typename BASETYPE, typename ...ARGS>
class CatalogueEntryBase
{
public:
  CatalogueEntryBase() = default;
  virtual ~CatalogueEntryBase() = default;
  CatalogueEntryBase(CatalogueEntryBase const &) = default;
  CatalogueEntryBase(CatalogueEntryBase &&) = default;
  CatalogueEntryBase& operator=(CatalogueEntryBase const &) = default;
  CatalogueEntryBase& operator=(CatalogueEntryBase &&) = default;

  virtual std::unique_ptr<BASETYPE> Allocate( std::string const & name, ARGS&... args ) const = 0;
};

template<typename BASETYPE, typename ...ARGS>
class CatalogInterface
{
public:
  typedef std::unordered_map<std::string, std::unique_ptr< CatalogueEntryBase<BASETYPE, ARGS...> > > CatalogueType;

  static CatalogueType& GetCatalogue()
  {
    return BASETYPE::GetCatalogue();
  }

  static std::unique_ptr<BASETYPE> Factory( const std::string& objectTypeName, ARGS&...args )
  {
    CatalogueEntryBase<BASETYPE, ARGS...> const * const entry = GetCatalogue().at( objectTypeName ).get();
    return entry->Allocate( objectTypeName, args... );
  }

  CatalogInterface() = delete;
  ~CatalogInterface() = delete;
  CatalogInterface(CatalogInterface const &) = delete;
  CatalogInterface(CatalogInterface &&) = delete;
  CatalogInterface& operator=(CatalogInterface const &) = delete;
  CatalogInterface& operator=(CatalogInterface &&) = delete;

};

template<typename TYPE, typename BASETYPE, typename ...ARGS>
class CatalogueEntry : public CatalogueEntryBase<BASETYPE, ARGS...>
{
public:
  CatalogueEntry() :
    CatalogueEntryBase<BASETYPE, ARGS...>()
  {}

  ~CatalogueEntry() override final = default;

  CatalogueEntry(CatalogueEntry const & source ) :
    CatalogueEntryBase<BASETYPE, ARGS...>(source)
  {}

  CatalogueEntry(CatalogueEntry && source ):
    CatalogueEntryBase<BASETYPE, ARGS...>(std::move(source))
  {}

  CatalogueEntry& operator=(CatalogueEntry const & source )
  {
    CatalogueEntryBase<BASETYPE, ARGS...>::operator=(source);
  }

  CatalogueEntry& operator=(CatalogueEntry && source )
  {
    CatalogueEntryBase<BASETYPE, ARGS...>::operator=(std::move(source));
  }

  virtual std::unique_ptr<BASETYPE> Allocate(  std::string const & name, ARGS&... args ) const override final
  {
    std::cout << "Creating "<< name <<" of type "<< geosx::stringutilities::demangle(typeid(TYPE).name())
              <<" from catalog of "<<geosx::stringutilities::demangle(typeid(BASETYPE).name())<<std::endl;

    return std::unique_ptr<BASETYPE>( new TYPE( args... ) );
  }
};

template<typename TYPE, typename BASETYPE, typename ...ARGS>
class CatalogueEntryConstructor
{
public:
  CatalogueEntryConstructor()
  {
    std::string name = TYPE::CatalogueName();
    ( CatalogInterface<BASETYPE, ARGS...>::GetCatalogue() ).insert( std::make_pair( name, std::make_unique< CatalogueEntry<TYPE,BASETYPE, ARGS...> >() ) );

    std::cout <<"Registered  "
              <<geosx::stringutilities::demangle(typeid(BASETYPE).name())
              <<" catalogue component of derived type "
              <<geosx::stringutilities::demangle(typeid(TYPE).name())
              <<" where "<<geosx::stringutilities::demangle(typeid(TYPE).name())<<"::CatalogueName() = "<<TYPE::CatalogueName() << std::endl;
  }

  ~CatalogueEntryConstructor() = default;
  CatalogueEntryConstructor(CatalogueEntryConstructor const &) = delete;
  CatalogueEntryConstructor(CatalogueEntryConstructor &&) = delete;
  CatalogueEntryConstructor& operator=(CatalogueEntryConstructor const &) = delete;
  CatalogueEntryConstructor& operator=(CatalogueEntryConstructor &&) = delete;

};

/// Compiler directive to simplify registration
#define REGISTER_CATALOGUE_ENTRY( BaseType, ClassName, ...) \
namespace \
{ \
objectcatalogue::CatalogueEntryConstructor<ClassName,BaseType,__VA_ARGS__> catEntry; \
}
//objectcatalogue::CatalogueEntry<ClassName,BaseType,__VA_ARGS__> catEntry; \
//objectcatalogue::CatalogueEntry<ClassName,BaseType,__VA_ARGS__> * catEntry = new objectcatalogue::CatalogueEntry<ClassName,BaseType,__VA_ARGS__>(); }
}





#define CATALOGUE( BASETYPE, PARAMS, ARGS )\
class CatalogueEntryBase\
{\
public:\
  typedef std::unordered_map<std::string, CatalogueEntryBase*> CatalogueType;\
\
  CatalogueEntryBase() {}\
  virtual std::unique_ptr<BASETYPE> Allocate( PARAMS ) = 0;\
  virtual ~CatalogueEntryBase() {}\
\
  static CatalogueType& GetCatalogue()\
  {\
    static CatalogueType * const catalogue = new CatalogueType();\
    return *catalogue;\
  }\
\
  static std::unique_ptr<BASETYPE> Factory( const std::string& objectTypeName, PARAMS )\
  {\
    std::cout << "Creating solver of type: " << objectTypeName << std::endl;\
    CatalogueEntryBase* const entry = GetCatalogue().at( objectTypeName );\
    return entry->Allocate( ARGS );\
  }\
};\
\
template<typename TYPE>\
class CatalogueEntry : public CatalogueEntryBase\
{\
public:\
  CatalogueEntry() : CatalogueEntryBase()\
  {\
    std::string name = TYPE::CatalogueName();\
    ( CatalogueEntryBase::GetCatalogue() )[name] = this;\
    std::cout << "Registered Solver: " << name << std::endl;\
  }\
\
  ~CatalogueEntry() final {}\
\
  virtual std::unique_ptr<BASETYPE> Allocate( PARAMS ) final\
  {\
    return std::make_unique<TYPE>( ARGS );\
  }\
};\
\

//#define REGISTER_CATALOGUE_ENTRY( BASETYPE, TYPE ) namespace{ BASETYPE::CatalogueEntry< TYPE > reg_##TYPE; }




#endif /* OBJECTCATALOGUE_HPP_ */
