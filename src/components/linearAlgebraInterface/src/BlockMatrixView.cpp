/*
 * BlockMatrixView.cpp
 *
 *  Created on: Aug 24, 2018
 *      Author: Matthias
 */

#include "BlockMatrixView.hpp"

namespace geosx
{

// Empty constructor
template< typename LAI >BlockMatrixView<LAI>::BlockMatrixView()
{
}

template< typename LAI >BlockMatrixView<LAI>::BlockMatrixView(array1d<ParallelMatrix> * Mats)
{
}


// Apply the block matrix to a block vector.
template< typename LAI > void BlockMatrixView<LAI>::apply()
{}

// Set to residual form.
template< typename LAI > void BlockMatrixView<LAI>::residual()
{}

// Clear row and multiply the diagonal entry by <tt>factor</tt>.
template< typename LAI > void BlockMatrixView<LAI>::clearRow(globalIndex rowIndex, real64 factor)
{}

// Accessor for block.
template< typename LAI > void BlockMatrixView<LAI>::getBlock(integer blockRowIndex, integer blockColIndex)
{}

// Accessor for block.
template< typename LAI > void BlockMatrixView<LAI>::getBlock(std::string blockName)
{}


}

