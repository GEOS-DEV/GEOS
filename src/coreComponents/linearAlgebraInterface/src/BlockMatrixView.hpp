/*
 * BlockMatrixView.hpp
 *
 *  Created on: Aug 24, 2018
 *      Author: Matthias
 */

#ifndef SRC_COMPONENTS_LINEARALGEBRAINTERFACE_SRC_BLOCKMATRIXVIEW_HPP_
#define SRC_COMPONENTS_LINEARALGEBRAINTERFACE_SRC_BLOCKMATRIXVIEW_HPP_

#include "TrilinosInterface.hpp"
//#include "HypreInterface.hpp"
#include "BlockVectorView.hpp"

namespace geosx
{

/**
 * \class BlockMatrixView
 * \brief This class creates and provides basic support for block
 *        matrices objects (templated on the LA interface).
 */

template< typename LAI >
class BlockMatrixView
{

  using ParallelMatrix = typename LAI::ParallelMatrix;
  using ParallelVector = typename LAI::ParallelVector;

public:

  //! @name Constructor/Destructor Methods
  //@{
  /**
   * @brief Empty matrix constructor.
   *
   * Create an empty block matrix.
   */
  BlockMatrixView();

  /**
   * @brief Empty matrix constructor.
   *
   * Create a block matrix of size (<tt>nRows</tt>,<tt>nCols</tt>).
   */
  BlockMatrixView( integer nRows, integer nCols );

  /**
   * @brief Virtual destructor.
   */
  virtual ~BlockMatrixView() = default;
  //@}

  //! @name Linear Algebra Methods
  //@{
  /**
   * @brief Apply the block matrix to a block vector.
   */
  void multiply(BlockVectorView<LAI> &solution,
                BlockVectorView<LAI> &rhs);

  /**
   * @brief Set to residual form.
   */
  //ParallelVector residual();

  /**
   * @brief Clear row and multiply the diagonal entry by <tt>factor</tt>.
   */
  void clearRow( globalIndex rowIndex, real64 factor );
  //@}

  //! @name Accessors/Setters
  //@{
  /**
   * @brief Get the matrix corresponding to block (<tt>blockRowIndex</tt>,<tt>blockColIndex</tt>).
   */
  ParallelMatrix * getBlock( integer blockRowIndex, integer blockColIndex );

  /**
   * @brief Get the matrix corresponding to block <tt>name</tt>.
   */
  ParallelMatrix * getBlock( std::string blockName );

  /**
   * @brief Set block (<tt>i</tt>,<tt>j</tt>) using <tt>matrix</tt>.
   */
  void setBlock( integer blockRowIndex, integer blockColIndex, ParallelMatrix &matrix );

  //@}

private:

  array2d<ParallelMatrix *> m_matrices;

};

// Empty constructor (inlined)
template< typename LAI >
inline
BlockMatrixView<LAI>::BlockMatrixView()
{}

// Constructor with a size (inlined)
template< typename LAI >
inline
BlockMatrixView<LAI>::BlockMatrixView( integer nRows, integer nCols )
{
  m_matrices.resize( nRows, nCols );
}



// Apply the block matrix to a block vector (hard coded to 2 by 2 for now).
template< typename LAI >
void BlockMatrixView<LAI>::multiply( BlockVectorView<LAI> &solution,
                                     BlockVectorView<LAI> &rhs )
{
  for( integer row = 0 ; row < m_matrices.size( 0 ) ; row++ )
  {
    rhs.scale( row, 0. );
    ParallelVector temp( *rhs.getBlock(row) );
    for( integer col = 0 ; col < m_matrices.size( 1 ) ; col++ )
    {
      if (m_matrices[row][col] != nullptr)
      {
        m_matrices[row][col]->multiply( *solution.getBlock(col), temp );
        rhs.update( row, 1.0, temp, 1.0 );
      }
    }
  }
}

//// Set to residual form.
//template< typename LAI >
//typename LAI::ParallelVector BlockMatrixView<LAI>::residual()
//{
//  integer numRows = 1;
//  integer numCols = 2;
//
//  ParallelVector res( *m_rhs[0] );
//
//  for( integer row = 0 ; row < numRows ; row++ )
//  {
//    for( integer col = 0 ; col < numCols - 1 ; col++ )
//    {
//      ParallelVector temp( *m_rhs[row] );
//      m_matrices[row][col]->multiply( *m_solution[col], temp );
//      m_matrices[row][col+1]->multiply( *m_solution[col+1], *m_rhs[row] );
//      m_rhs[row]->update( 1.0, temp, 1.0 );
//      res.update( -1.0, *m_rhs[row], 1.0 );
//    }
//    return res;
//  }
//}

// Clear row and multiply the diagonal entry by <tt>factor</tt>.
template< typename LAI >
void BlockMatrixView<LAI>::clearRow( globalIndex rowIndex, real64 factor )
{}

// Accessor for block.
template< typename LAI >
typename LAI::ParallelMatrix * BlockMatrixView<LAI>::getBlock( integer blockRowIndex, integer blockColIndex )
{
  return m_matrices[blockRowIndex][blockColIndex];
}

// Setter for block.
template< typename LAI >
void BlockMatrixView<LAI>::setBlock( integer blockRowIndex, integer blockColIndex, typename LAI::ParallelMatrix &matrix )
{
  m_matrices[blockRowIndex][blockColIndex] = &matrix;
}

}


#endif /* SRC_COMPONENTS_LINEARALGEBRAINTERFACE_SRC_BLOCKMATRIXVIEW_HPP_ */
