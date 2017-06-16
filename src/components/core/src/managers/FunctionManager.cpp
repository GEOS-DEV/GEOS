/*
 * FunctionManager.cpp
 *
 *  Created on: June 16, 2017
 *      Author: sherman
 */

#include "FunctionManager.hpp"
#include "dataRepository/ManagedGroup.hpp"
#include "common/DataTypes.hpp"
#include <mathpresso/mathpresso.h>

namespace geosx
{

namespace dataRepository
{
namespace keys
{
  std::string const variableNames = "variableNames";
  std::string const expression = "expression";
}
}

using namespace dataRepository;


// JIT Compiled Functions
JIT_Function::JIT_Function( const std::string& name,
                            ManagedGroup * const parent ) :
  ManagedGroup( name, parent ),
  parserContext(),
  parserExpression()
{}

JIT_Function::~JIT_Function()
{
  // TODO Auto-generated destructor stub
}

void JIT_Function::FillDocumentationNode( dataRepository::ManagedGroup * const domain )
{
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();
  
  docNode->setName(this->CatalogName());
  docNode->setSchemaType("Node");
  docNode->setShortDescription("JIT function");

  docNode->AllocateChildNode( keys::variableNames,
                              keys::variableNames,
                              -1,
                              "string_array",
                              "string_array",
                              "List of variables in expression",
                              "List of variables in expression.  The order must match the evaluate argment.",
                              "",
                              "",
                              1,
                              1,
                              0 );

  docNode->AllocateChildNode( keys::expression,
                              keys::expression,
                              -1,
                              "string",
                              "string",
                              "Symbolic math expression",
                              "Symbolic math expression",
                              "default",
                              "",
                              1,
                              1,
                              0 );

}

void JIT_Function::BuildDataStructure( ManagedGroup * const domain )
{
  JIT_Function::BuildDataStructure( domain );
  RegisterDocumentationNodes();
}

void JIT_Function::Initialize( ManagedGroup& functionManager )
{
  // Register variables
  view_rtype<string_array> variables = nodes.getData<string_array>(keys::variableNames);
  for (int ii=0; ii<variables.size(); ++ii)
  {
    parserContext.addVariable(variables[ii].c_str(), ii * sizeof(double));
  }

  // Add built in constants/functions (PI, E, sin, cos, ceil, exp, etc.), compile
  parserContext.addBuiltIns();
  std::string expression = getData<std::string>(keys::JIT_Expression);
  mathpresso::Error err = parserExpression.compile(parserContext, expression.c_str(), mathpresso::kNoOptions);
  if (err != mathpresso::kErrorOk) 
  {
    throw std::invalid_argument("JIT compile error: " + err);
  }
}



// Function Manager
FunctionManager::FunctionManager( const std::string& name,
                                  ManagedGroup * const parent ) :
  ManagedGroup( name, parent )
{}

FunctionManager::~FunctionManager()
{
  // TODO Auto-generated destructor stub
}

void FunctionManager::FillDocumentationNode( dataRepository::ManagedGroup * const /*group*/ )
{
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();
  SolverBase::FillDocumentationNode( domain );
  
  docNode->setName(this->CatalogName());
  docNode->setSchemaType("Node");
  docNode->setShortDescription("JIT function manager");
}

void FunctionManager::ReadXML( dataRepository::ManagedGroup& domain,
                               xmlWrapper::xmlNode const & problemNode )
{
  xmlWrapper::xmlNode topLevelNode = problemNode.child("Functions");
  if (topLevelNode != NULL)
  {
    std::cout << "Functions:" << std::endl;

    for (xmlWrapper::xmlNode functionNode=topLevelNode.first_child(); functionNode; functionNode=functionNode.next_sibling())
    {
      // Register the new solver
      std::cout << "   " << functionNode.name() << std::endl;
      std::string functionID = functionNode.attribute("name").value();
      JIT_Function & newFunction = CreateSolver( functionNode.name(), functionID );

      // Set the documentation node, register, and read xml
      newFunction.SetDocumentationNodes( &domain );
      newFunction.BuildDataStructure( &domain );
      newFunction.ReadXML(functionNode );
    }
  }
}

JIT_Function & FunctionManager::CreateFunction( string const & functionCatalogKey, string const & functionName )
{
  std::unique_ptr<JIT_Function> function = JIT_Function::CatalogInterface::Factory( functionCatalogKey, functionName, this );
  JIT_Function & rval = this->RegisterGroup<JIT_Function>( functionName, std::move(function) );

  return rval;
}

REGISTER_CATALOG_ENTRY( ObjectManagerBase, FunctionManager, std::string const &, ManagedGroup * const )
} /* namespace ANST */
