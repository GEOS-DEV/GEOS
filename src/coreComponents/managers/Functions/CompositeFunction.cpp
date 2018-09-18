/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file CompositeFunction.cpp
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
                                      ManagedGroup * const parent ):
  FunctionBase( name, parent ),
  parserContext(),
  parserExpression(),
  m_numSubFunctions(),
  m_subFunctions()
{}


CompositeFunction::~CompositeFunction()
{}


void CompositeFunction::FillDocumentationNode()
{
  FunctionBase::FillDocumentationNode();
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

  docNode->getChildNode(keys::inputVarNames)->setDefault("");
  docNode->getChildNode(keys::inputVarTypes)->setDefault("");

}


void CompositeFunction::InitializeFunction()
{
  // Register variables
  string_array const & variables = getReference<string_array>(keys::variableNames);
  for (string_array::size_type ii=0 ; ii<variables.size() ; ++ii)
  {
    parserContext.addVariable(variables[ii].c_str(), static_cast<int>(ii * sizeof(double)));
  }

  // Add built in constants/functions (PI, E, sin, cos, ceil, exp, etc.),
  // compile
  parserContext.addBuiltIns();
  std::string expression = getData<std::string>(keys::expression);
  mathpresso::Error err = parserExpression.compile(parserContext, expression.c_str(), mathpresso::kNoOptions);
  if (err != mathpresso::kErrorOk)
  {
    throw std::invalid_argument("JIT Compiler Error");
  }

  // Grab pointers to sub functions
  NewFunctionManager * functionManager = NewFunctionManager::Instance();
  string_array const & functionNames = getReference<string_array>(keys::functionNames);
  m_numSubFunctions = integer_conversion<localIndex>(functionNames.size());
  for (localIndex ii=0 ; ii<m_numSubFunctions ; ++ii)
  {
    m_subFunctions.push_back(functionManager->GetGroup<FunctionBase>(functionNames[ii]));
  }
}


void CompositeFunction::Evaluate( dataRepository::ManagedGroup const * const group,
                                  real64 const time,
                                  set<localIndex> const & set,
                                  real64_array & result ) const
{
  // Evaluate each of the subFunctions independently and place the results into
  // a temporary field
  array1d<real64_array> subFunctionResults;
  for (localIndex ii=0 ; ii<m_numSubFunctions ; ++ii)
  {
    real64_array tmp(result.size());
    m_subFunctions[ii]->Evaluate(group, time, set, tmp);
    subFunctionResults.push_back(std::move(tmp));
  }

  // Evaluate the symbolic math
  real64 functionResults[m_maxNumSubFunctions];
  for( auto const & ii : set )
  {
    for (localIndex jj=0 ; jj<m_numSubFunctions ; ++jj)
    {
      functionResults[jj] = subFunctionResults[jj][ii];
    }
    result[ii] = parserExpression.evaluate( reinterpret_cast<void*>( functionResults ));
  }
}


real64 CompositeFunction::Evaluate( real64 const * const input ) const
{
  real64 functionResults[m_maxNumSubFunctions];

  for (localIndex ii=0 ; ii<m_numSubFunctions ; ++ii)
  {
    functionResults[ii] = m_subFunctions[ii]->Evaluate(input);
  }

  return parserExpression.evaluate( reinterpret_cast<void*>( functionResults ));
}


REGISTER_CATALOG_ENTRY( FunctionBase, CompositeFunction, std::string const &, ManagedGroup * const )

} /* namespace ANST */
