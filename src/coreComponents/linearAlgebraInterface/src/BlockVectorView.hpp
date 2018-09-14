/*
 * BlockMatrixView.hpp
 *
 *  Created on: Aug 24, 2018
 *      Author: Matthias
 */

#ifndef SRC_COMPONENTS_LINEARALGEBRAINTERFACE_SRC_BLOCKVECTORVIEW_HPP_
#define SRC_COMPONENTS_LINEARALGEBRAINTERFACE_SRC_BLOCKVECTORVIEW_HPP_

#include "TrilinosInterface.hpp"
//#include "HypreInterface.hpp"

namespace geosx
{

/**
 * \class BlockVectorView
 * \brief This class creates and provides basic support for block
 *        vectors objects (templated on the LA interface).
 */

template< typename LAI >
class BlockVectorView
{

  using ParallelVector = typename LAI::ParallelVector;

public:

  //! @name Constructor/Destructor Methods
  //@{
  /**
   * @brief Empty matrix constructor.
   *
   * Create an empty block matrix.
   */
  BlockVectorView();

  /**
   * @brief Empty matrix constructor.
   *
   * Create a block matrix of size (<tt>nRows</tt>,<tt>nCols</tt>).
   */
  BlockVectorView( integer nBlocks );

  /**
   * @brief Virtual destructor.
   */
  virtual ~BlockVectorView() = default;
  //@}

  //! @name Linear Algebra Operations
  //@{
  /**
   * @brief Get the matrix corresponding to block (<tt>blockRowIndex</tt>,<tt>blockColIndex</tt>).
   */
  void scale( integer blockIndex, real64 factor );

  /**
   * @brief 2-norm of the block vector.
   */
  void norm2( real64 &result );

  /**
   * @brief Get the matrix corresponding to block (<tt>blockRowIndex</tt>,<tt>blockColIndex</tt>).
   */
  void update( integer blockIndex, real64 const alpha, EpetraVector const &vec, real64 const beta );

  //@}

  //! @name Accessors/Setters
  //@{
  /**
   * @brief Get the matrix corresponding to block (<tt>blockRowIndex</tt>,<tt>blockColIndex</tt>).
   */
  ParallelVector * getBlock( integer blockIndex );

  /**
   * @brief Get the matrix corresponding to block <tt>name</tt>.
   */
  ParallelVector * getBlock( std::string blockName );

  /**
   * @brief Set block (<tt>i</tt>,<tt>j</tt>) using <tt>matrix</tt>.
   */
  void setBlock( integer blockIndex, ParallelVector &vector );

  //@}

private:

  array1d<ParallelVector *> m_vectors;

};

// Empty constructor (inlined)
template< typename LAI >
inline
BlockVectorView<LAI>::BlockVectorView()
{}

// Constructor with a size (inlined)
template< typename LAI >
inline
BlockVectorView<LAI>::BlockVectorView( integer nBlocks )
{
  m_vectors.resize( nBlocks );
}

// Scale a block.
template< typename LAI >
void BlockVectorView<LAI>::scale( integer blockIndex, real64 factor )
{
  m_vectors[blockIndex]->scale(factor);
}

// Compute the 2 norm of the block vector.
template< typename LAI >
void BlockVectorView<LAI>::norm2( real64 &result )
{
  real64 accum = 0;
  for (integer i = 0; i < m_vectors.size(); i++)
  {
    real64 temp;
    m_vectors[i]->norm2(temp);
    accum = accum + temp*temp;
  }
  result = std::sqrt(accum);
}

// Update a block.
template< typename LAI >
void BlockVectorView<LAI>::update( integer blockIndex, real64 const alpha, EpetraVector const &vec, real64 const beta )
{
  m_vectors[blockIndex]->update( alpha, vec, beta  );
}

// Setter for block.
template< typename LAI >
void BlockVectorView<LAI>::setBlock( integer blockIndex, typename LAI::ParallelVector &vector )
{
  m_vectors[blockIndex] = &vector;
}

// Accessor for block.
template< typename LAI >
typename LAI::ParallelVector * BlockVectorView<LAI>::getBlock( integer blockIndex )
{
  return m_vectors[blockIndex];
}

}


#endif /* SRC_COMPONENTS_LINEARALGEBRAINTERFACE_SRC_BLOCKMATRIXVIEW_HPP_ */
