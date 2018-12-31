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
 * @file SymbolicFunction.cpp
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
                                    ManagedGroup * const parent ):
  FunctionBase( name, parent ),
  parserContext(),
  parserExpression()
{
  RegisterViewWrapper<string_array>(keys::variableNames)->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("List of variables in expression.  The order must match the evaluate argument");

  RegisterViewWrapper<string>(keys::expression)->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Symbolic math expression");
}


SymbolicFunction::~SymbolicFunction()
{}

void SymbolicFunction::InitializeFunction()
{
  // Register variables
  string_array & variables = getReference<string_array>(keys::variableNames);
  for ( localIndex ii=0 ; ii<variables.size(); ++ii)
  {
    parserContext.addVariable(variables[ii].c_str(), static_cast<int>(ii * sizeof(double)));
  }

  // Add built in constants/functions (PI, E, sin, cos, ceil, exp, etc.),
  // compile
  parserContext.addBuiltIns();
  std::string const& expression = getReference<std::string>(keys::expression);
  mathpresso::Error err = parserExpression.compile(parserContext, expression.c_str(), mathpresso::kNoOptions);
  GEOS_ERROR_IF(err != mathpresso::kErrorOk, "JIT Compiler Error");
}


REGISTER_CATALOG_ENTRY( FunctionBase, SymbolicFunction, std::string const &, ManagedGroup * const )

} /* namespace ANST */
