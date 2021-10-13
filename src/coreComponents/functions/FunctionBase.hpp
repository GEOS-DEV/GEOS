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
 * @file FunctionBase.hpp
 */

#ifndef GEOSX_FUNCTIONS_FUNCTIONBASE_HPP_
#define GEOSX_FUNCTIONS_FUNCTIONBASE_HPP_

#include "common/DataTypes.hpp"
#include "dataRepository/Group.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

namespace geosx
{

namespace dataRepository
{
namespace keys
{
/**
 * @brief The key for inputVarNames
 * @return the key
 */
string const inputVarNames( "inputVarNames" );
}
}

/**
 * @class FunctionBase
 *
 * An object for interfacing with arbitrary N-dimensional functions.
 */
class FunctionBase : public dataRepository::Group
{
public:
  /// @copydoc geosx::dataRepository::Group::Group( string const & name, Group * const parent )
  FunctionBase( const string & name,
                dataRepository::Group * const parent );

  /**
   * @brief destructor
   */
  virtual ~FunctionBase() override;

  /**
   * @brief Static Factory Catalog Functions
   * @return the catalog name
   */
  static string catalogName() { return "FunctionBase"; }

  /**
   * @brief Function initialization
   */
  virtual void initializeFunction(){}

  /**
   * @brief Test to see if the function is a 1D function of time
   * @return integer value:
   *         - 0 is the function does not have time as parameter
   *         - 1 is the function has time as one of the parameters
   *         - 2 is the function has time as only parameter
   */
  integer isFunctionOfTime() const;

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
                         real64_array & result ) const = 0;

  /**
   * @brief Method to evaluate a function
   * @param input a scalar input
   * @return the function evaluation
   */
  virtual real64 evaluate( real64 const * const input ) const = 0;

  /// Alias for the catalog interface
  using CatalogInterface = dataRepository::CatalogInterface< FunctionBase, string const &, Group * const >;

  /**
   * @brief return the catalog entry for the function
   * @return the catalog entry
   */
  static CatalogInterface::CatalogType & getCatalog()
  {
    static CatalogInterface::CatalogType catalog;
    return catalog;
  }

  /**
   * @brief This generates statistics by applying a function to an object
   * @param group a pointer to the object holding the function arguments
   * @param time current time
   * @param set the subset of nodes to apply the function to
   * @return An array holding the min, average, max values of the results
   */
  real64_array evaluateStats( dataRepository::Group const & group,
                              real64 const time,
                              SortedArray< localIndex > const & set ) const;

  /**
   * @brief Set the input variable names
   * @param inputVarNames A list of input variable names
   */
  void setInputVarNames( string_array inputVarNames ) { m_inputVarNames = inputVarNames; }


protected:
  /// names for the input variables
  string_array m_inputVarNames;

  /**
   * @brief Method to apply an function with an arbitrary type of output
   * @tparam LEAF the return type
   * @param[in] group a pointer to the object holding the function arguments
   * @param[in] time current time
   * @param[in] set the subset of nodes to apply the function to
   * @param[out] result the results
   */
  template< typename LEAF >
  void evaluateT( dataRepository::Group const & group,
                  real64 const time,
                  SortedArrayView< localIndex const > const & set,
                  real64_array & result ) const;

  virtual void postProcessInput() override { initializeFunction(); }

};

template< typename LEAF >
void FunctionBase::evaluateT( dataRepository::Group const & group,
                              real64 const time,
                              SortedArrayView< localIndex const > const & set,
                              real64_array & result ) const
{
  real64 const * input_ptrs[4];
  localIndex varSize[4] = {0, 0, 0, 0};
  int timeVar[4] = {1, 1, 1, 1};

  arrayView1d< string const > const & inputVarNames = this->getReference< string_array >( dataRepository::keys::inputVarNames );
  localIndex const numVars = LvArray::integerConversion< localIndex >( inputVarNames.size());
  localIndex groupSize = group.size();
  localIndex totalVarSize = 0;
  for( auto varIndex=0; varIndex<numVars; ++varIndex )
  {
    string const & varName = inputVarNames[varIndex];

    if( varName=="time" )
    {
      input_ptrs[varIndex] = &time;
      varSize[varIndex] = 1;
      timeVar[varIndex] = 0;
      ++totalVarSize;
    }
    else if( groupSize > 0 )
    {
      // Should we throw a warning if the group is zero-length?
      dataRepository::WrapperBase const & wrapper = group.getWrapperBase( varName );
      input_ptrs[ varIndex ] = reinterpret_cast< double const * >( wrapper.voidPointer() );

      localIndex wrapperSize = LvArray::integerConversion< localIndex >( wrapper.size() );
      varSize[varIndex] = wrapperSize / groupSize;
      totalVarSize += varSize[varIndex];
    }
  }

  // Make sure the inputs do not exceed the maximum length
  GEOSX_ERROR_IF( totalVarSize > 4, "Function input size is: " << totalVarSize );

  // Make sure the result / set size match
  GEOSX_ERROR_IF( result.size() != set.size(), "To apply a function to a set, the size of the result and set must match" );


  forAll< serialPolicy >( set.size(), [&, set]( localIndex const i )
  {
    localIndex const index = set[ i ];
    double input[4];
    int c = 0;
    for( int a=0; a<numVars; ++a )
    {
      for( int b=0; b<varSize[a]; ++b )
      {
        input[c] = input_ptrs[a][(index*varSize[a]+b)*timeVar[a]];
        ++c;
      }
    }

    // Note: we expect that result is the same size as the set
    result[i] = static_cast< LEAF const * >(this)->evaluate( input );
  } );
}
} /* namespace geosx */

#endif /* GEOSX_FUNCTIONS_FUNCTIONBASE_HPP_ */
