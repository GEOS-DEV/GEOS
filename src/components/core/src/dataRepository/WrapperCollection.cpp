/*
 * DataObjectManager.cpp
 *
 *  Created on: Nov 21, 2014
 *      Author: rrsettgast
 */

#include "WrapperCollection.hpp"
#include "dataRepository/SidreWrapper.hpp"

namespace geosx
{
namespace dataRepository
{

WrapperCollection::WrapperCollection( std::string const & name,
                                      WrapperCollection * const parent ) :
  m_keyLookup(),
  m_wrappers(),
  m_parent(parent),
  m_subObjectManagers(),
  m_sidreGroup(nullptr)
{
  asctoolkit::sidre::DataGroup * sidreParent = nullptr;
  if( m_parent==nullptr )
  {
    sidreParent = SidreWrapper::dataStore().getRoot();
  }
  else
  {
    sidreParent = parent->m_sidreGroup;
  }

  if( sidreParent->hasGroup(name) )
  {
    m_sidreGroup = sidreParent->getGroup(name);
  }
  else
  {
    m_sidreGroup = sidreParent->createGroup(name);
  }

  *(RegisterWrapper<std_size_t>( "size" ).data()) = 0;
  RegisterWrapper<std::string>( "name" ).reference() = name;
  RegisterWrapper<std::string>( "path" );


}

WrapperCollection::~WrapperCollection()
{
  // TODO Auto-generated destructor stub
}

//DataObjectManager::DataObjectManager( DataObjectManager const & source ):
//    m_size( source.m_size ),
//    m_name( source.m_name ),
//    m_path( source.m_path ),
//    m_keyLookup( source.m_keyLookup ),
//    m_dataObjects( source.m_dataObjects ),
//    m_parent( source.m_parent )
//{}

WrapperCollection::WrapperCollection( WrapperCollection&& source ) :
  m_keyLookup( std::move(source.m_keyLookup) ),
  m_wrappers( std::move(source.m_wrappers) ),
  m_parent( std::move(source.m_parent) )
{}

WrapperCollection::CatalogInterface::CatalogType& WrapperCollection::GetCatalog()
{
  static WrapperCollection::CatalogInterface::CatalogType catalog;
  return catalog;
}

WrapperBase& WrapperCollection::RegisterWrapper( std::string const & name, rtTypes::TypeIDs const & type )
{
  return *( rtTypes::ApplyTypeLambda( type,
                                      [this, &name]( auto a ) -> WrapperBase*
      {
        return &( this->RegisterWrapper<decltype(a)>(name) );
      } ) );
}

void WrapperCollection::resize( std::size_t const newsize )
{
  for( auto&& i : this->m_wrappers )
  {
    i->resize(newsize);
  }
  *(this->getWrapper<std_size_t>("size").data())=newsize;
}


}
} /* namespace ODS */
