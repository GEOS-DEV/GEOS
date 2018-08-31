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
 * @file SymbolicFunction.hpp
 */

#ifndef SYMBOLICFUNCTION_HPP_
#define SYMBOLICFUNCTION_HPP_

#include "FunctionBase.hpp"
#include <mathpresso/mathpresso.h>

namespace geosx
{

/**
 * @class SymbolicFunction
 *
 * An interface for an arbitrary symbolic function
 */
class SymbolicFunction : public FunctionBase
{
public:
  /// Main constructor
  SymbolicFunction( const std::string& name,
                    dataRepository::ManagedGroup * const parent );

  /// Destructor
  virtual ~SymbolicFunction() override;

  /// Catalog name interface
  static string CatalogName() { return "SymbolicFunction"; }

  /// Documentation assignment
  virtual void FillDocumentationNode() override;
  
  /// Function initialization
  virtual void InitializeFunction() override;

  /**
   * @brief Method to evaluate a function on a target object
   * @param group a pointer to the object holding the function arguments
   * @param time current time
   * @param set the subset of nodes to apply the function to
   * @param result an array to hold the results of the function
   */
  inline void Evaluate( dataRepository::ManagedGroup const * const group,
                        real64 const time,
                        set<localIndex> const & set,
                        real64_array & result ) const override final
  {
    FunctionBase::EvaluateT<SymbolicFunction>( group, time, set, result );
  }

  /**
   * @brief Method to evaluate a function
   * @param input a scalar input
   */
  inline real64 Evaluate( real64 const * const input ) const override final
  {
    return parserExpression.evaluate( reinterpret_cast<void*>( const_cast<real64*>(input) ) );
  }

private:
  /// Symbolic math driver objects
  mathpresso::Context parserContext;
  mathpresso::Expression parserExpression;
};


} /* namespace geosx */

#endif /* SYMBOLICFUNCTION_HPP_ */
