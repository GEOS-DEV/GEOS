/*
 * BlockMatrixView.hpp
 *
 *  Created on: Aug 24, 2018
 *      Author: Matthias
 */

#ifndef SRC_COMPONENTS_LINEARALGEBRAINTERFACE_SRC_BLOCKMATRIXVIEW_HPP_
#define SRC_COMPONENTS_LINEARALGEBRAINTERFACE_SRC_BLOCKMATRIXVIEW_HPP_

#include "TrilinosInterface.hpp"

namespace geosx
{

/**
 * \class BlockMatrixView
 * \brief This class creates and provides basic support for Trilinos-based block
 *        matrices objects.
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
   * Create a block matrix from an array of matrices.
   */
  //BlockMatrixView(array1d<ParallelMatrix> * Mats);

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
  void apply();

  /**
   * @brief Set to residual form.
   */
  ParallelVector residual();

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
   * @brief Get the solution vector corresponding to block (<tt>j</tt>).
   */
  ParallelVector * getSolution( integer blockColIndex );

  /**
   * @brief Get the rhs vector corresponding to block (<tt>j</tt>).
   */
  ParallelVector * getRhs( integer blockRowIndex );

  /**
   * @brief Get the matrix corresponding to block <tt>name</tt>.
   */
  ParallelMatrix * getBlock( std::string blockName );

  /**
   * @brief Set block (<tt>i</tt>,<tt>j</tt>) using <tt>matrix</tt>.
   */
  void setBlock( integer blockRowIndex, integer blockColIndex, ParallelMatrix * matrix );

  /**
   * @brief Get the solution vector corresponding to block (<tt>j</tt>).
   */
  void setSolution( integer blockColIndex, ParallelVector * vector );

  /**
   * @brief Get the rhs vector corresponding to block (<tt>j</tt>).
   */
  void setRhs( integer blockRowIndex, ParallelVector * vector );

  //@}

private:

  ParallelMatrix * m_matrices[2][2];
  ParallelVector * m_solution[2];
  ParallelVector * m_rhs[2];

};

// Empty constructor (inlined)
template< typename LAI >
inline
BlockMatrixView<LAI>::BlockMatrixView()
{}

// Apply the block matrix to a block vector (hard coded to 2 by 2 for now).
template< typename LAI >
void BlockMatrixView<LAI>::apply()
{
  integer numRows = 1;
  integer numCols = 2;

  for (integer row = 0; row < numRows; row++)
  {
    for (integer col = 0; col < numCols - 1; col++)
    {
      ParallelVector temp( *m_rhs[row] );
      m_matrices[row][col]->multiply( *m_solution[col], temp );
      m_matrices[row][col+1]->multiply( *m_solution[col+1], *m_rhs[row] );
      m_rhs[row]->update( 1.0, temp, 1.0 );
    }
  }
}

// Set to residual form.
template< typename LAI >
typename LAI::ParallelVector BlockMatrixView<LAI>::residual()
{
  integer numRows = 1;
  integer numCols = 2;

  ParallelVector res( *m_rhs[0] );

  for (integer row = 0; row < numRows; row++)
  {
    for (integer col = 0; col < numCols - 1; col++)
    {
      ParallelVector temp( *m_rhs[row] );
      m_matrices[row][col]->multiply( *m_solution[col], temp );
      m_matrices[row][col+1]->multiply( *m_solution[col+1], *m_rhs[row] );
      m_rhs[row]->update( 1.0, temp, 1.0 );
      res.update( -1.0, *m_rhs[row], 1.0 );
    }
    return res;
  }
}

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
void BlockMatrixView<LAI>::setBlock( integer blockRowIndex, integer blockColIndex, typename LAI::ParallelMatrix * matrix )
{
  m_matrices[blockRowIndex][blockColIndex] = matrix;
}

// Accessor for block.
template< typename LAI >
typename LAI::ParallelVector * BlockMatrixView<LAI>::getSolution( integer blockColIndex )
{
  return m_solution[blockColIndex];
}

// Accessor for block.
template< typename LAI >
typename LAI::ParallelVector * BlockMatrixView<LAI>::getRhs( integer blockRowIndex )
{
  return m_rhs[blockRowIndex];
}

// Setter for solution.
template< typename LAI >
void BlockMatrixView<LAI>::setSolution( integer blockColIndex, typename LAI::ParallelVector * vector )
{
  m_solution[blockColIndex] = vector;
}

// Setter for rhs.
template< typename LAI >
void BlockMatrixView<LAI>::setRhs( integer blockRowIndex, typename LAI::ParallelVector * vector )
{
  m_rhs[blockRowIndex] = vector;
}

}


#endif /* SRC_COMPONENTS_LINEARALGEBRAINTERFACE_SRC_BLOCKMATRIXVIEW_HPP_ */
