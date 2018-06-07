// Copyright (c) 2018, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-746361. All Rights
// reserved. See file COPYRIGHT for details.
//
// This file is part of the GEOSX Simulation Framework.

//
// GEOSX is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.
/*
 * SymbolicFunction.hpp
 *
 *  Created on: July 6, 2017
 *      Author: sherman
 */

#ifndef SYMBOLICFUNCTION_HPP_
#define SYMBOLICFUNCTION_HPP_

#include "FunctionBase.hpp"
#include <mathpresso/mathpresso.h>

namespace geosx
{


class SymbolicFunction : public FunctionBase
{
public:
  SymbolicFunction( const std::string& name,
                    dataRepository::ManagedGroup * const parent );

  virtual ~SymbolicFunction() override;
  static string CatalogName() { return "SymbolicFunction"; }
  virtual void FillDocumentationNode() override;
  virtual void BuildDataStructure( dataRepository::ManagedGroup * const domain ) override;

  virtual void InitializeFunction() override;

  inline void Evaluate( dataRepository::ManagedGroup const * const group,
                        real64 const time,
                        lSet const & set,
                        real64_array & result ) const override final
  {
    FunctionBase::EvaluateT<SymbolicFunction>( group, time, set, result );
  }


  inline real64 Evaluate( real64 const * const input ) const override final
  {
    return parserExpression.evaluate( reinterpret_cast<void*>( const_cast<real64*>(input) ) );
  }

private:
  mathpresso::Context parserContext;
  mathpresso::Expression parserExpression;
};


} /* namespace geosx */

#endif /* SYMBOLICFUNCTION_HPP_ */
