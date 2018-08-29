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
  BlockMatrixView(array1d<ParallelMatrix> * Mats);

  /**
   * @brief Virtual destructor.
   */
  virtual ~BlockMatrixView() = default;
  //@}

  //! @name Linear Algebra Methods
  //@{
  /**
   * @brief Apply the block matrix to a block vector.
   *
   */
  void apply();

  /**
   * @brief Set to residual form.
   */
  void residual();

  /**
   * @brief Clear row and multiply the diagonal entry by <tt>factor</tt>.
   */
  void clearRow(globalIndex rowIndex, real64 factor);
  //@}

  //! @name Accessors
  //@{
  /**
   * @brief Get the matrix corresponding to block (<tt>i</tt>,<tt>j</tt>).
   *
   */
  void getBlock(integer blockRowIndex, integer blockColIndex);

  /**
   * @brief Get the matrix corresponding to block <tt>name</tt>.
   *
   */
  void getBlock(std::string blockName);
  //@}

};

// Empty constructor (inlined)
template< typename LAI >
inline
BlockMatrixView<LAI>::BlockMatrixView()
{
}

template< typename LAI >
inline
BlockMatrixView<LAI>::BlockMatrixView(array1d<ParallelMatrix> * Mats)
{
}

// Apply the block matrix to a block vector.
template< typename LAI >
inline
void BlockMatrixView<LAI>::apply()
{}

// Set to residual form.
template< typename LAI >
inline
void BlockMatrixView<LAI>::residual()
{}

// Clear row and multiply the diagonal entry by <tt>factor</tt>.
template< typename LAI >
inline
void BlockMatrixView<LAI>::clearRow(globalIndex rowIndex, real64 factor)
{}

// Accessor for block.
template< typename LAI >
inline
void BlockMatrixView<LAI>::getBlock(integer blockRowIndex, integer blockColIndex)
{}

// Accessor for block.
template< typename LAI >
inline
void BlockMatrixView<LAI>::getBlock(std::string blockName)
{}


}


#endif /* SRC_COMPONENTS_LINEARALGEBRAINTERFACE_SRC_BLOCKMATRIXVIEW_HPP_ */
