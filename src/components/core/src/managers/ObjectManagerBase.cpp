/*
 * ObjectManagerBase.cpp
 *
 *  Created on: Sep 15, 2016
 *      Author: settgast1
 */

#include "ObjectManagerBase.hpp"

namespace geosx
{
using namespace dataRepository;

ObjectManagerBase::ObjectManagerBase( std::string const & name,
                                      ManagedGroup * const parent ):
    ManagedGroup(name,parent),
    m_localToGlobalMap( RegisterViewWrapper< gArray1d >("localToGlobal").reference() )
{


  this->RegisterGroup<ManagedGroup>("Sets");
  this->RegisterViewWrapper< iArray1d >("isExternal");
}
//ObjectManagerBase::ObjectManagerBase( std::string const & name,
//                                      ManagedGroup * const parent,
//                                      cxx_utilities::DocumentationNode * docNode ):
//    ManagedGroup(name,parent,docNode),
//    m_localToGlobalMap( RegisterViewWrapper< gArray1d >("localToGlobal").reference() )
//{
//
//
//  this->RegisterGroup<ManagedGroup>("Sets");
//  this->RegisterViewWrapper< iArray1d >("isExternal");
//}


ObjectManagerBase::~ObjectManagerBase()
{}



ObjectManagerBase::CatalogInterface::CatalogType& ObjectManagerBase::GetCatalog()
{
  static ObjectManagerBase::CatalogInterface::CatalogType catalog;
  return catalog;
}


void ObjectManagerBase::ConstructSetFromSetAndMap( const lSet& inputSet,
                                                          const lArray2d& map,
                                                          const std::string& newSetName )
{

  ManagedGroup& sets = GetGroup(std::string("Sets"));
  lSet& newset = sets.RegisterViewWrapper<lSet>(newSetName).reference();
  newset.clear();

  int mapSize = map.Dimension(1);
  for( localIndex ka=0 ; ka<size() ; ++ka )
  {
    const localIndex* const sublist = map[ka];
    int addToSet = 0;
    for( int a=0 ; a<mapSize ; ++a )
    {
      if( inputSet.count( sublist[a] ) == 1 )
      {
        ++addToSet;
      }
    }
    if( addToSet == mapSize )
    {
      newset.insert( ka );
    }
  }
}

void ObjectManagerBase::ConstructSetFromSetAndMap( const lSet& inputSet,
                                                          const Array1dT<lArray1d>& map,
                                                          const std::string& newSetName )
{

  ManagedGroup& sets = GetGroup(std::string("Sets"));
  lSet& newset = sets.RegisterViewWrapper<lSet>(newSetName).reference();
  newset.clear();

  for( localIndex ka=0 ; ka<size() ; ++ka )
  {
    int addToSet = 0;
    int mapSize = map[ka].size();
    for( int a=0 ; a<mapSize ; ++a )
    {
      if( inputSet.count( map[ka][a] ) == 1 )
      {
        ++addToSet;
      }
    }
    if( addToSet == mapSize )
    {
      newset.insert( ka );
    }
  }
}

void ObjectManagerBase::ConstructListOfBoundaryObjects( lArray1d& objectList ) const
{
  const iArray1d& isDomainBoundary = this->GetFieldData<int>("isDomainBoundary");
  for( localIndex k=0 ; k<size() ; ++k )
  {
    if( isDomainBoundary[k] == 1 )
    {
      objectList.push_back( k );
    }
  }
}

void ObjectManagerBase::ConstructListOfBoundaryObjects( gArray1d& objectList ) const
{
  const iArray1d& isDomainBoundary = this->GetFieldData<int>("isDomainBoundary");
  for( localIndex k=0 ; k<size() ; ++k )
  {
    if( isDomainBoundary[k] == 1 )
    {
      objectList.push_back( this->m_localToGlobalMap[k] );
    }
  }
  std::sort( objectList.begin(), objectList.end() );
}

} /* namespace geosx */
