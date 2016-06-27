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

  template< typename BASETYPE, typename ...ARGS >
  class ObjectCatalogueEntryBase
  {
  public:
    typedef std::map< std::string, ObjectCatalogueEntryBase<BASETYPE,ARGS... >* > CatalogueType;

    ObjectCatalogueEntryBase(){}
    virtual std::unique_ptr<BASETYPE> Allocate( ARGS&... args ) = 0;
    virtual ~ObjectCatalogueEntryBase(){}

    static CatalogueType& GetCatalogue()
    {
      static CatalogueType catalogue;
      return catalogue;
    }

    static std::unique_ptr<BASETYPE> Factory( const std::string& objectTypeName, ARGS&...args )
    {
      std::cout<<"Creating solver of type: "<<objectTypeName<<std::endl;
      ObjectCatalogueEntryBase<BASETYPE,ARGS...>* const  entry = GetCatalogue().at(objectTypeName);
      return entry->Allocate(args...);
    }

  };


  template< typename TYPE, typename BASETYPE, typename...ARGS >
  class ObjectCatalogueEntry : public ObjectCatalogueEntryBase<BASETYPE,ARGS...>
  {
  public:
    ObjectCatalogueEntry():
      ObjectCatalogueEntryBase<BASETYPE,ARGS...>()
    {
      std::string name = TYPE::CatalogueName();
      (ObjectCatalogueEntryBase<BASETYPE,ARGS...>::GetCatalogue())[name] = this;
      std::cout<<"Registered Solver: "<<name<<std::endl;
    }

    ~ObjectCatalogueEntry() final{}

    virtual std::unique_ptr<BASETYPE> Allocate( ARGS&... args) final
    {
      return std::unique_ptr<BASETYPE>(new TYPE(args...));
    }


};




/// Compiler directive to simplify autoregistration
#define REGISTER_FACTORY( ClassName, BaseType, ARGS...) namespace{ ObjectCatalogueEntry<ClassName,BaseType,ARGS> reg_; }

#endif /* OBJECTCATALOGUE_HPP_ */
