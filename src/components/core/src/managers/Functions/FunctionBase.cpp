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
 * @file FunctionBase.cpp
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
{}


void FunctionBase::FillDocumentationNode()
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


real64_array FunctionBase::EvaluateStats( dataRepository::ManagedGroup const * const group,
                                          real64 const time,
                                          lSet const & set) const
{
  localIndex N = set.size();
  real64_array sub(N);
  Evaluate( group, time, set, sub );

  real64_array result(3);
  result[0] = 1e10;   // min
  result[1] = 0.0;    // avg
  result[2] = -1e10;  // max
  for( localIndex ii=0; ii<N; ii++ )
  {
    result[0] = std::min(result[0], sub[ii]);
    result[1] += sub[ii];
    result[2] = std::max(result[0], sub[ii]);
  }
  result[1] /= N;

  return result;
}


} /* namespace ANST */
