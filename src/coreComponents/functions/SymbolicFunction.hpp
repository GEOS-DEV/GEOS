/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */
#ifndef GEOSX_FUNCTIONS_SYMBOLICFUNCTION_HPP_
#define GEOSX_FUNCTIONS_SYMBOLICFUNCTION_HPP_

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
  /// @copydoc geosx::dataRepository::Group::Group( string const & name, Group * const parent )
  SymbolicFunction( const string & name,
                    dataRepository::Group * const parent );

  /**
   * @brief The destructor
   */
  virtual ~SymbolicFunction() override;

  /**
   * @brief The catalog name interface
   * @return name of the TableFunction in the FunctionBase catalog
   */
  static string catalogName() { return "SymbolicFunction"; }

  /**
   * @brief Initialize the table function
   */
  virtual void initializeFunction() override;

  /**
   * @brief Method to evaluate a function on a target object
   * @param group a pointer to the object holding the function arguments
   * @param time current time
   * @param set the subset of nodes to apply the function to
   * @param result an array to hold the results of the function
   */
  inline void evaluate( dataRepository::Group const & group,
                        real64 const time,
                        SortedArrayView< localIndex const > const & set,
                        real64_array & result ) const override final
  {
    FunctionBase::evaluateT< SymbolicFunction >( group, time, set, result );
  }

  /**
   * @brief Method to evaluate a function
   * @param input a scalar input
   * @return the function result
   */
  inline real64 evaluate( real64 const * const input ) const override final
  {
    return parserExpression.evaluate( reinterpret_cast< void * >( const_cast< real64 * >(input) ) );
  }


  /**
   * @brief Set the symbolic variable names
   * @param variableNames An array of variable names used in the expression
   */
  void setSymbolicVariableNames( string_array variableNames ) { m_variableNames = variableNames; }

  /**
   * @brief Set the symbolic expression
   * @param expression A string containing the symbolic expression
   */
  void setSymbolicExpression( string expression ) { m_expression = expression; }



private:
  // Symbolic math driver objects
  mathpresso::Context parserContext;
  mathpresso::Expression parserExpression;

  /// Symbolic expression variable names
  string_array m_variableNames;

  /// Symbolic expression
  string m_expression;
};


} /* namespace geosx */

#endif /* GEOSX_FUNCTIONS_SYMBOLICFUNCTION_HPP_ */
