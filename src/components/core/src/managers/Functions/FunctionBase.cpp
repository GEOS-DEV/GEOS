/*
 * FunctionBase.cpp
 *
 *  Created on: July 6, 2017
 *      Author: sherman
 */

#include "FunctionBase.hpp"
#include "common/DataTypes.hpp"
#include <mathpresso/mathpresso.h>

namespace geosx
{

using namespace dataRepository;



FunctionBase::FunctionBase( const std::string& name,
                            ManagedGroup * const parent ) :
  ManagedGroup( name, parent )
{}

FunctionBase::~FunctionBase()
{
  // TODO Auto-generated destructor stub
}


void FunctionBase::FillDocumentationNode( dataRepository::ManagedGroup * const domain )
{
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();
  
  docNode->setName(this->CatalogName());
  docNode->setSchemaType("Node");
  docNode->setShortDescription("Function Base");

}



REGISTER_CATALOG_ENTRY( FunctionBase, FunctionBase, std::string const &, ManagedGroup * const )

} /* namespace ANST */
