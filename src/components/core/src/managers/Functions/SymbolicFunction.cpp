/*
 * SymbolicFunction.cpp
 *
 *  Created on: July 6, 2017
 *      Author: sherman
 */

#include "SymbolicFunction.hpp"
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



SymbolicFunction::SymbolicFunction( const std::string& name,
                                    ManagedGroup * const parent ) :
  FunctionBase( name, parent ),
  parserContext(),
  parserExpression()
{}

SymbolicFunction::~SymbolicFunction()
{
  // TODO Auto-generated destructor stub
}


void SymbolicFunction::FillDocumentationNode( dataRepository::ManagedGroup * const domain )
{
  FunctionBase::FillDocumentationNode(domain);
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

void SymbolicFunction::BuildDataStructure( ManagedGroup * const domain )
{
  RegisterDocumentationNodes();
}

void SymbolicFunction::InitializeFunction()
{
  // Register variables
  string_array & variables = getReference<string_array>(keys::variableNames);
  for (int ii=0; ii<variables.size(); ++ii)
  {
    parserContext.addVariable(variables[ii].c_str(), ii * sizeof(double));
  }

  // Add built in constants/functions (PI, E, sin, cos, ceil, exp, etc.), compile
  parserContext.addBuiltIns();
  std::string expression = getData<std::string>(keys::expression);
  mathpresso::Error err = parserExpression.compile(parserContext, expression.c_str(), mathpresso::kNoOptions);
  if (err != mathpresso::kErrorOk) 
  {
    throw std::invalid_argument("JIT Compiler Error");
  }
}


REGISTER_CATALOG_ENTRY( FunctionBase, SymbolicFunction, std::string const &, ManagedGroup * const )

} /* namespace ANST */
