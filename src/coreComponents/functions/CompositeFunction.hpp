/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file CompositeFunction.hpp
 */

#ifndef GEOSX_FUNCTIONS_COMPOSITEFUNCTION_HPP_
#define GEOSX_FUNCTIONS_COMPOSITEFUNCTION_HPP_

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
  /// @copydoc geosx::dataRepository::Group::Group( string const & name, Group * const parent )
  CompositeFunction( const string & name,
                     dataRepository::Group * const parent );

  /**
   * @brief destructor
   */
  virtual ~CompositeFunction() override;

  /**
   * @brief Static Factory Catalog Functions
   * @return the catalog name
   */
  static string catalogName() { return "CompositeFunction"; }

  /**
   * @brief Function initialization
   */
  virtual void initializeFunction() override;

  /**
   * @brief Method to evaluate a function on a target object
   * @param group a pointer to the object holding the function arguments
   * @param time current time
   * @param set the subset of nodes to apply the function to
   * @param result an array to hold the results of the function
   */
  virtual void evaluate( dataRepository::Group const & group,
                         real64 const time,
                         SortedArrayView< localIndex const > const & set,
                         real64_array & result ) const override final;

  /**
   * @brief Method to evaluate a function
   * @param input a scalar input
   * @return the function evaluation
   */
  virtual real64 evaluate( real64 const * const input ) const override final;

private:

  string_array m_functionNames;
  string_array m_variableNames;
  string m_expression;

  mathpresso::Context parserContext;
  mathpresso::Expression parserExpression;

  localIndex m_numSubFunctions;
  static constexpr localIndex m_maxNumSubFunctions = 10;
  std::vector< FunctionBase * > m_subFunctions;

};

} /* namespace geosx */

#endif /* GEOSX_FUNCTIONS_COMPOSITEFUNCTION_HPP_ */
