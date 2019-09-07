/*
* ------------------------------------------------------------------------------------------------------------
* SPDX-License-Identifier: LGPL-2.1-only
*
* Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
* Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
* Copyright (c) 2018-2019 Total, S.A
* Copyright (c) 2019-     GEOSX Contributors
* All right reserved
*
* See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
* ------------------------------------------------------------------------------------------------------------
*/

/**
* @file TableFunction.hpp
*/

#ifndef TABLEFUNCTION_HPP_
#define TABLEFUNCTION_HPP_

#include "FunctionBase.hpp"

namespace geosx
{

/**
* @class TableFunction
*
* An interface for a dense table-based function
*/
class TableFunction : public FunctionBase
{
public:
/// Main constructor
TableFunction( const std::string& name,
dataRepository::Group * const parent );

/// Destructor
virtual ~TableFunction() override;

/// Catalog name interface
static string CatalogName() { return "TableFunction"; }

/// Initialize the function
virtual void InitializeFunction() override;

void reInitializeFunction();

/**
* @brief Method to evaluate a function on a target object
* @param group a pointer to the object holding the function arguments
* @param time current time
* @param set the subset of nodes to apply the function to
* @param result an array to hold the results of the function
*/
virtual void Evaluate( dataRepository::Group const * const group,
real64 const time,
SortedArrayView<localIndex const> const & set,
real64_array & result ) const override final
{
FunctionBase::EvaluateT<TableFunction>( group, time, set, result );
}

/**
* @brief Method to evaluate a function
* @param input a scalar input
*/
virtual real64 Evaluate( real64 const * const input) const override final;


array1d<real64_array> const & getCoordinates() const { return m_coordinates; }
array1d<real64_array>       & getCoordinates()       { return m_coordinates; }

array1d<real64> const & getValues() const { return m_values; }
array1d<real64>       & getValues()       { return m_values; }

private:
/// An array of table axes
array1d<real64_array> m_coordinates;

/// Table values (in fortran order)
real64_array m_values;

/// Table size indicators
static localIndex constexpr m_maxDimensions = 4;
localIndex m_dimensions;
localIndex_array m_size;
localIndex_array m_indexIncrement;

// m_corners should be of size m_maxDimensions x (2^m_maxDimensions)
localIndex m_corners[m_maxDimensions][16];
localIndex m_numCorners;
};


} /* namespace geosx */

#endif /* TABLEFUNCTION_HPP_ */
