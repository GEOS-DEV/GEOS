/*
 * NewFunctionManager.cpp
 *
 *  Created on: June 16, 2017
 *      Author: sherman
 */

#include "NewFunctionManager.hpp"
#include "dataRepository/ManagedGroup.hpp"
#include "FunctionBase.hpp"
#include "common/DataTypes.hpp"

namespace geosx
{

using namespace dataRepository;


NewFunctionManager::NewFunctionManager( const std::string& name,
                                        ManagedGroup * const parent ) :
  ManagedGroup( name, parent )
{}

NewFunctionManager::~NewFunctionManager()
{
  // TODO Auto-generated destructor stub
}

void NewFunctionManager::FillDocumentationNode( dataRepository::ManagedGroup * const /*group*/ )
{
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();
  
  docNode->setName("Functions");
  docNode->setSchemaType("Node");
  docNode->setShortDescription("Function manager");
}

void NewFunctionManager::ReadXML( dataRepository::ManagedGroup * domain,
                                  xmlWrapper::xmlNode const & problemNode )
{
  xmlWrapper::xmlNode topLevelNode = problemNode.child("Functions");
  if (topLevelNode != NULL)
  {
    std::cout << "Functions:" << std::endl;

    for (xmlWrapper::xmlNode functionNode=topLevelNode.first_child(); functionNode; functionNode=functionNode.next_sibling())
    {
      // Register the new function
      FunctionBase * newFunction = CreateFunction(functionNode.name(), functionNode.attribute("name").value());

      // Set the documentation node, register, and read xml
      newFunction->SetDocumentationNodes( this );
      newFunction->BuildDataStructure( this );
      newFunction->ReadXML(functionNode );
      newFunction->InitializeFunction();
    }
  }
}

FunctionBase * NewFunctionManager::CreateFunction( string const & functionCatalogKey,
                                                   string const & functionName )
{
  std::cout << "   " << functionCatalogKey << ": " << functionName << std::endl;
  std::unique_ptr<FunctionBase> function = FunctionBase::CatalogInterface::Factory( functionCatalogKey, functionName, this );
  FunctionBase * rval = this->RegisterGroup<FunctionBase>( functionName, std::move(function) );

  return rval;
}

} /* namespace ANST */
