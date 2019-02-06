/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
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
  ManagedGroup( name, parent ),
  m_inputVarNames()
{
  setInputFlags(InputFlags::OPTIONAL_NONUNIQUE);

  RegisterViewWrapper( keys::inputVarNames, &m_inputVarNames, 0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setSizedFromParent(0)->
    setDescription("Name of fields are input to function.");
}


FunctionBase::~FunctionBase()
{}

integer FunctionBase::isFunctionOfTime() const
{
  integer rval=0;
  string_array const & inputVarNames = this->getReference<string_array>( dataRepository::keys::inputVarNames );
  localIndex numVars = integer_conversion<localIndex>(inputVarNames.size());

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
                                          set<localIndex> const & set) const
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
