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
 * @file BlockMatrixView.hpp
 *
 *  Created on: Aug 24, 2018
 *      Author: Matthias Cremon
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
   * @brief Matrix of (<tt>nRows</tt>,<tt>nCols</tt>) blocks.
   *
   * Create a block matrix of size (<tt>nRows</tt>,<tt>nCols</tt>).
   */
  BlockMatrixView( localIndex const nRows,
                   localIndex const nCols );

  /**
   * @brief Destructor.
   */
  ~BlockMatrixView() = default;
  //@}

  //! @name Linear Algebra Methods
  //@{
  /**
   * @brief Apply the block matrix to a block vector.
   *
   * Computes the matrix-vector product <tt>Ax = b</tt>.
   *
   * \param x Input vector.
   * \param b Output vector.
   *
   */
  void multiply( BlockVectorView<LAI> const &x,
                 BlockVectorView<LAI> &b ) const;

  /**
   * @brief Compute residual <tt>r = b - Ax</tt>.
   *
   * \param x Input solution.
   * \param b Input right hand size.
   * \param r Output residual.
   *
   */
  void residual( BlockVectorView<LAI> const &x,
                 BlockVectorView<LAI> const &b,
                 BlockVectorView<LAI> &r ) const;

  /**
   * @brief Scale matrix using <tt>factor</tt>.
   */
  void scale( real64 const factor );


  //@}
  //! @name Accessors/Setters
  //@{

  /**
   * @brief Get the matrix corresponding to block (<tt>blockRowIndex</tt>,<tt>blockColIndex</tt>).
   */
  ParallelMatrix & block( localIndex const blockRowIndex,
                          localIndex const blockColIndex ) const;

  /**
   * @brief Set block (<tt>i</tt>,<tt>j</tt>) using <tt>matrix</tt>.
   */
  void set( localIndex const blockRowIndex,
            localIndex const blockColIndex,
            ParallelMatrix &matrix );

  //@}

private:

  array2d<ParallelMatrix *> m_matrices;

};

}// end geosx namespace


#endif /* SRC_COMPONENTS_LINEARALGEBRAINTERFACE_SRC_BLOCKMATRIXVIEW_HPP_ */
