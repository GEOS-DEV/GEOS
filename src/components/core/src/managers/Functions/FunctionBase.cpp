/*
 * FunctionBase.cpp
 *
 *  Created on: July 6, 2017
 *      Author: sherman
 */

#include "FunctionBase.hpp"
#include "common/DataTypes.hpp"

namespace geosx
{

using namespace dataRepository;



FunctionBase::FunctionBase( const std::string& name,
                            ManagedGroup * const parent ):
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

  docNode->AllocateChildNode( keys::inputVarNames,
                              keys::inputVarNames,
                              -1,
                              "string_array",
                              "string_array",
                              "Name of fields are input to function.",
                              "",
                              "REQUIRED",
                              "",
                              0,
                              1,
                              0 );

  docNode->AllocateChildNode( keys::inputVarTypes,
                              keys::inputVarTypes,
                              -1,
                              "string_array",
                              "string_array",
                              "Name of fields are input to function.",
                              "",
                              "REQUIRED",
                              "",
                              0,
                              1,
                              0 );


}

integer FunctionBase::isFunctionOfTime() const
{
  integer rval=0;
  string_array const & inputVarNames = this->getReference<string_array>( dataRepository::keys::inputVarNames );
  localIndex numVars = inputVarNames.size();

  if( numVars==1 )
  {
    if( inputVarNames[0]=="time" )
    {
      rval = 2;
    }
  }
  else
  {
    for( auto varIndex=0 ; varIndex<numVars ; ++varIndex )
    {
      if( inputVarNames[varIndex]=="time")
      {
        rval = 1;
        break;
      }
    }
  }
  return rval;
}



//REGISTER_CATALOG_ENTRY( FunctionBase, FunctionBase, std::string const &,
// ManagedGroup * const )

} /* namespace ANST */
