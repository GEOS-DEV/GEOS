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
  typedef std::unordered_map<std::string, CatalogueEntryBase<BASETYPE, ARGS...>*> CatalogueType;

  CatalogueEntryBase()
  {
  }

  virtual std::unique_ptr<BASETYPE> Allocate( ARGS&... args ) = 0;
  virtual ~CatalogueEntryBase()
  {
  }

//  static CatalogueType& GetCatalogue()
//  {
//    static CatalogueType * const catalogue = new CatalogueType();
//    return *catalogue;
//  }
  static CatalogueType& GetCatalogue()
  {
    return BASETYPE::GetCatalogue();
  }


  static std::unique_ptr<BASETYPE> Factory( const std::string& objectTypeName, ARGS&...args )
  {
//    std::cout << "Creating solver of type: " << objectTypeName << std::endl;
    std::cout << "Creating "<< geosx::stringutilities::demangle(typeid(BASETYPE).name()) <<" of type: " << objectTypeName << std::endl;

    CatalogueEntryBase<BASETYPE, ARGS...>* const entry = GetCatalogue().at( objectTypeName );
    return entry->Allocate( args... );
  }

};

template<typename TYPE, typename BASETYPE, typename ...ARGS>
class CatalogueEntry : public CatalogueEntryBase<BASETYPE, ARGS...>
{
public:
  CatalogueEntry() :
    CatalogueEntryBase<BASETYPE, ARGS...>()
  {
    std::string name = TYPE::CatalogueName();
    ( CatalogueEntryBase<BASETYPE, ARGS...>::GetCatalogue() )[name] = this;
    std::cout <<"Registered  "
              <<geosx::stringutilities::demangle(typeid(BASETYPE).name())
              <<" catalogue component of derived type "
              <<geosx::stringutilities::demangle(typeid(TYPE).name())
              <<" named "<< name << std::endl;

  }

  ~CatalogueEntry() final
  {
  }

  virtual std::unique_ptr<BASETYPE> Allocate( ARGS&... args ) final
  {
    return std::unique_ptr<BASETYPE>( new TYPE( args... ) );
  }

};

/// Compiler directive to simplify registration
#define REGISTER_CATALOGUE_ENTRY( BaseType, ClassName, ...) namespace{ objectcatalogue::CatalogueEntry<ClassName,BaseType,__VA_ARGS__> reg_; }
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
