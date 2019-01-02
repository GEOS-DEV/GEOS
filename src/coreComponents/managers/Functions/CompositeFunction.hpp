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
 * @file CompositeFunction.hpp
 */

#ifndef COMPOSITEFUNCTION_HPP_
#define COMPOSITEFUNCTION_HPP_

#include "FunctionBase.hpp"
#include <mathpresso/mathpresso.h>

namespace geosx
{

/**
 * @class CompositeFunction
 *
 * An interface for combinined functions
 */
class CompositeFunction : public FunctionBase
{
public:
  /// Main constructor
  CompositeFunction( const std::string& name,
                     dataRepository::ManagedGroup * const parent );

  /// Destructor
  virtual ~CompositeFunction() override;

  /// Catalog name interface
  static string CatalogName() { return "CompositeFunction"; }
  
  /// Function initialization
  virtual void InitializeFunction() override;

  /**
   * @brief Method to evaluate a function on a target object
   * @param group a pointer to the object holding the function arguments
   * @param time current time
   * @param set the subset of nodes to apply the function to
   * @param result an array to hold the results of the function
   */
  virtual void Evaluate( dataRepository::ManagedGroup const * const group,
                         real64 const time,
                         set<localIndex> const & sets,
                         real64_array & result ) const override final;

  /**
   * @brief Method to evaluate a function
   * @param input a scalar input
   */
  virtual real64 Evaluate( real64 const * const input) const override final;

private:

  string_array m_functionNames;
  string_array m_variableNames;
  string       m_expression;

  mathpresso::Context parserContext;
  mathpresso::Expression parserExpression;

  localIndex m_numSubFunctions;
  static constexpr localIndex m_maxNumSubFunctions = 10;
  std::vector<FunctionBase*> m_subFunctions;

};


} /* namespace geosx */

#endif /* COMPOSITEFUNCTION_HPP_ */
