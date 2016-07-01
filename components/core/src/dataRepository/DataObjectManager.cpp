/*
 * DataObjectManager.cpp
 *
 *  Created on: Nov 21, 2014
 *      Author: rrsettgast
 */

#include "DataObjectManager.hpp"

namespace geosx
{

DataObjectManager::DataObjectManager( std::string const & name ):
m_size(0),
m_name(name),
m_path(),
m_keyLookup(),
m_dataObjects(),
m_parent(nullptr),
m_subObjectManagers()
{}

DataObjectManager::~DataObjectManager()
{
  // TODO Auto-generated destructor stub
}


DataObjectManager::DataObjectManager( DataObjectManager&& other )
{
  m_name = std::move( other.m_name );
  m_keyLookup = std::move( other.m_keyLookup );
  m_dataObjects = std::move( other.m_dataObjects );
}


void DataObjectManager::resize( const std::size_t newsize )
{

  m_size = newsize;
}



} /* namespace ODS */
