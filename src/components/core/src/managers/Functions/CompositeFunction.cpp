/*
 * CompositeFunction.cpp
 *
 *  Created on: August 17, 2017
 *      Author: sherman
 */

#include "NewFunctionManager.hpp"
#include "CompositeFunction.hpp"
#include "common/DataTypes.hpp"
#include <mathpresso/mathpresso.h>

namespace geosx
{

namespace dataRepository
{
namespace keys
{
  std::string const functionNames = "functionNames";
  std::string const variableNames = "variableNames";
  std::string const expression = "expression";
}
}

using namespace dataRepository;



CompositeFunction::CompositeFunction( const std::string& name,
                                    ManagedGroup * const parent ) :
  FunctionBase( name, parent ),
  parserContext(),
  parserExpression(),
  m_numSubFunctions(),
  m_subFunctions()
{}

CompositeFunction::~CompositeFunction()
{
  // TODO Auto-generated destructor stub
}


void CompositeFunction::FillDocumentationNode( dataRepository::ManagedGroup * const domain )
{
  FunctionBase::FillDocumentationNode(domain);
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();
  
  docNode->setName(this->CatalogName());
  docNode->setSchemaType("Node");
  docNode->setShortDescription("Composite function");

  docNode->AllocateChildNode( keys::functionNames,
                              keys::functionNames,
                              -1,
                              "string_array",
                              "string_array",
                              "List of source functions",
                              "List of source functions.  The order must match the variableNames argment.",
                              "",
                              "",
                              1,
                              1,
                              0 );

  docNode->AllocateChildNode( keys::variableNames,
                              keys::variableNames,
                              -1,
                              "string_array",
                              "string_array",
                              "List of variables in expression",
                              "List of variables in expression",
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
                              "Composite math expression",
                              "Composite math expression",
                              "default",
                              "",
                              1,
                              1,
                              0 );

}

void CompositeFunction::BuildDataStructure( ManagedGroup * const domain )
{
  RegisterDocumentationNodes();
}

void CompositeFunction::InitializeFunction()
{
  // Register variables
  view_rtype<string_array> variables = getData<string_array>(keys::variableNames);
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

  // Grab pointers to sub functions
  NewFunctionManager& functionManager = NewFunctionManager::Instance();
  view_rtype<string_array> functionNames = getData<string_array>(keys::functionNames);
  m_numSubFunctions = functionNames.size();
  for (localIndex ii=0; ii<m_numSubFunctions; ++ii)
  {
    m_subFunctions.push_back(&(functionManager.GetGroup<FunctionBase>(functionNames[ii])));
  }
}


void CompositeFunction::Evaluate( dataRepository::ManagedGroup const * const group,
                                  real64 const time,
                                  lSet const & set,
                                  real64_array & result ) const
{
  // Evaluate each of the subFunctions independently and place the results into a temporary field
  Array1dT<real64_array> subFunctionResults;
  for (localIndex ii=0; ii<m_numSubFunctions; ++ii)
  {
    real64_array tmp(result.size());
    m_subFunctions[ii]->Evaluate(group, time, set, tmp);
    subFunctionResults.emplace_back(std::move(tmp));
  }

  // Evaluate the symbolic math
  real64 functionResults[m_numSubFunctions];
  for( auto const & ii : set )
  {
    for (localIndex jj=0; jj<m_numSubFunctions; ++jj)
    {
      functionResults[jj] = subFunctionResults[jj][ii];
    }
    result[ii] = parserExpression.evaluate( reinterpret_cast<void*>( functionResults ));
  }
}


real64 CompositeFunction::Evaluate( real64 const * const input ) const
{
  real64 functionResults[m_numSubFunctions];

  for (localIndex ii=0; ii<m_numSubFunctions; ++ii)
  {
    functionResults[ii] = m_subFunctions[ii]->Evaluate(input);
  }

  return parserExpression.evaluate( reinterpret_cast<void*>( functionResults ));
}


REGISTER_CATALOG_ENTRY( FunctionBase, CompositeFunction, std::string const &, ManagedGroup * const )

} /* namespace ANST */
