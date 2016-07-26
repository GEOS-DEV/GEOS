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
  if( parent==nullptr )
  {
    m_sidreGroup = SidreWrapper::dataStore().getRoot();
  }
  else
  {
    if( parent->m_sidreGroup->hasGroup(name) )
    {
      m_sidreGroup = parent->m_sidreGroup->getGroup(name);
    }
    else
    {
      m_sidreGroup = parent->m_sidreGroup->createGroup(name);
    }
  }

  *(RegisterWrapper<std_size_t>( "size" ).data()) = 0;
  std::string& temp = RegisterWrapper<std::string>( "name" ).dataRef();
  temp = name;
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
