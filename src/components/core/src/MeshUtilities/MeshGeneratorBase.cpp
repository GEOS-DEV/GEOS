/*
 * MeshGeneratorBase.cpp
 *
 *  Created on: Oct 17, 2017
 *      Author: sherman
 */

#include "MeshGeneratorBase.hpp"


namespace geosx
{
using namespace dataRepository;

MeshGeneratorBase::MeshGeneratorBase( string const & name, ManagedGroup * const parent ) :
  ManagedGroup( name, parent )
{}

MeshGeneratorBase::~MeshGeneratorBase()
{}

MeshGeneratorBase::CatalogInterface::CatalogType& MeshGeneratorBase::GetCatalog()
{
  static MeshGeneratorBase::CatalogInterface::CatalogType catalog;
  return catalog;
}

}
