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
   * @brief Empty vector constructor.
   *
   * Create an empty block vector.
   */
  BlockVectorView();

  /**
   * @brief Sized vector constructor.
   *
   * Create a block vector of size (<tt>nRows</tt>,<tt>nCols</tt>).
   */
  BlockVectorView( typename LAI::laiLID nBlocks );

  /**
   * @brief Copy constructor.
   *
   * Create an empty block vector.
   */
  BlockVectorView( BlockVectorView<LAI> const &in_blockVec );

  /**
   * @brief Virtual destructor.
   */
  virtual ~BlockVectorView() = default;
  //@}

  //! @name Linear Algebra Operations
  //@{

  /**
   * @brief Scale the block vector.
   */
  void scale( real64 factor );

  /**
   * @brief Scale the vector corresponding to block (<tt>blockIndex</tt>).
   */
  void scale( typename LAI::laiLID blockIndex,
              real64 factor );

  /**
   * @brief Dot product.
   */
  void dot( BlockVectorView<LAI> const &blockVec,
            real64 &result );

  /**
   * @brief 2-norm of the block vector.
   */
  void norm2( real64 &result );

  /**
   * @brief 2-norm of the vector corresponding to block (<tt>blockIndex</tt>).
   */
  void norm2( typename LAI::laiLID blockIndex,
              real64 &result );

  /**
   * @brief Update the block vector.
   */
  void update( real64 const alpha,
               BlockVectorView<LAI> const &blockVec,
               real64 const beta );

  /**
   * @brief Update the vector corresponding to block (<tt>blockIndex</tt>).
   */
  void update( typename LAI::laiLID blockIndex,
               real64 const alpha,
               ParallelVector const &vec,
               real64 const beta );

  //@}

  //! @name Accessors/Setters
  //@{

  /**
   * @brief Get block size.
   */
  typename LAI::laiLID blockSize() const;

  /**
   * @brief Get global size.
   */
  TrilinosInterface::laiGID globalSize();

  /**
   * @brief Get the vector corresponding to block (<tt>blockRowIndex</tt>,<tt>blockColIndex</tt>).
   */
  ParallelVector * getBlock( typename LAI::laiLID blockIndex ) const;

  /**
   * @brief Get the vector corresponding to block <tt>name</tt>.
   */
  ParallelVector * getBlock( std::string blockName ) const;

  /**
   * @brief Set block (<tt>i</tt>,<tt>j</tt>) using <tt>matrix</tt>.
   */
  void setBlock( typename LAI::laiLID blockIndex,
                 ParallelVector &vector );

  //@}

private:

  array1d<ParallelVector *> m_vectors;

};

// Empty constructor (inlined)
template< typename LAI >
inline
BlockVectorView<LAI>::BlockVectorView()
{}

// Constructor with a size
template< typename LAI >
BlockVectorView<LAI>::BlockVectorView( typename LAI::laiLID nBlocks )
{
  m_vectors.resize( nBlocks );
}

// Copy constructor
template< typename LAI >
BlockVectorView<LAI>::BlockVectorView( BlockVectorView<LAI> const &in_blockVec )
{
  m_vectors.resize( in_blockVec.blockSize() );
  for( typename LAI::laiLID i = 0 ; i < m_vectors.size() ; i++ )
  {
    //ParallelVector tempVec( *in_blockVec.getBlock( i ) );
    m_vectors[i] = new ParallelVector ( *in_blockVec.getBlock( i ));
  }
}

// Scale a block.
template< typename LAI >
void BlockVectorView<LAI>::scale( real64 factor )
{
  for( typename LAI::laiLID i = 0 ; i < m_vectors.size() ; i++ )
    m_vectors[i]->scale( factor );
}

// Scale a block.
template< typename LAI >
void BlockVectorView<LAI>::scale( typename LAI::laiLID blockIndex,
                                  real64 factor )
{
  m_vectors[blockIndex]->scale( factor );
}

// Dot product.
template< typename LAI >
void BlockVectorView<LAI>::dot( BlockVectorView<LAI> const &blockVec,
                                real64 &result )
{
  real64 accum = 0;
  for( typename LAI::laiLID i = 0 ; i < m_vectors.size() ; i++ )
  {
    real64 temp = 0;
    m_vectors[i]->dot( *blockVec.getBlock( i ), temp );
    accum = accum + temp;
  }
  result = accum;
}

// Compute the 2 norm of one block vector.
template< typename LAI >
void BlockVectorView<LAI>::norm2( typename LAI::laiLID blockIndex,
                                  real64 &result )
{
  m_vectors[blockIndex]->norm2( result );
}

// Compute the 2 norm of the block vector.
template< typename LAI >
void BlockVectorView<LAI>::norm2( real64 &result )
{
  real64 accum = 0;
  for( typename LAI::laiLID i = 0 ; i < m_vectors.size() ; i++ )
  {
    real64 temp;
    m_vectors[i]->norm2( temp );
    accum = accum + temp*temp;
  }
  result = std::sqrt( accum );
}

// Update a the vector.
template< typename LAI >
void BlockVectorView<LAI>::update( real64 const alpha,
                                   BlockVectorView<LAI> const &blockVec,
                                   real64 const beta )
{
  for( typename LAI::laiLID i = 0 ; i < m_vectors.size() ; i++ )
    m_vectors[i]->update( alpha, *blockVec.getBlock( i ), beta );
}

// Update a block.
template< typename LAI >
void BlockVectorView<LAI>::update( typename LAI::laiLID blockIndex,
                                   real64 const alpha,
                                   ParallelVector const &vec,
                                   real64 const beta )
{
  m_vectors[blockIndex]->update( alpha, vec, beta );
}

// Setter for block.
template< typename LAI >
void BlockVectorView<LAI>::setBlock( typename LAI::laiLID blockIndex,
                                     ParallelVector &vector )
{
  m_vectors[blockIndex] = &vector;
}

// Accessor for block size
template< typename LAI >
typename LAI::laiLID BlockVectorView<LAI>::blockSize() const
{
  return m_vectors.size();
}

// Accessor for size.
template< typename LAI >
TrilinosInterface::laiGID BlockVectorView<LAI>::globalSize(  )
{
  TrilinosInterface::laiGID size = 0;
  for( typename LAI::laiLID i = 0 ; i < m_vectors.size() ; i++ )
    size = m_vectors[i]->globalSize() + size;
  return size;
}

// Accessor for block.
template< typename LAI >
typename LAI::ParallelVector * BlockVectorView<LAI>::getBlock( typename LAI::laiLID blockIndex ) const
{
  return m_vectors[blockIndex];
}

}


#endif /* SRC_COMPONENTS_LINEARALGEBRAINTERFACE_SRC_BLOCKMATRIXVIEW_HPP_ */
