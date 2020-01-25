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
 */

#ifndef GEOSX_LINEARALGEBRA_UTILITIES_BLOCKOPERATORVIEW_HPP_
#define GEOSX_LINEARALGEBRA_UTILITIES_BLOCKOPERATORVIEW_HPP_

#include "linearAlgebra/interfaces/LinearOperator.hpp"
#include "linearAlgebra/utilities/BlockVectorView.hpp"

namespace geosx
{

/**
 * @brief This class creates and provides basic support for block operator objects.
 * @tparam Vector type of vector that sub-blocks of this view can operate on
 */

template< typename VECTOR, typename OPERATOR = LinearOperator<VECTOR> >
class BlockOperatorView : public LinearOperator< BlockVectorView<VECTOR> >
{
  using Vector = typename LinearOperator< BlockVectorView<VECTOR> >::Vector;
  using Operator = OPERATOR;

public:

  //! @name Constructor/Destructor Methods
  //@{
  /**
   * @brief Empty matrix constructor.
   *
   * Create an empty block matrix.
   */
  BlockOperatorView();

  /**
   * @brief Matrix of (<tt>nRows</tt>,<tt>nCols</tt>) blocks.
   *
   * Create a block matrix of size (<tt>nRows</tt>,<tt>nCols</tt>).
   */
  BlockOperatorView( localIndex const nRows,
                     localIndex const nCols );

  /**
   * @brief Destructor.
   */
  virtual ~BlockOperatorView() override = default;
  //@}

  //! @name Linear Algebra Methods
  //@{
  /**
   * @brief Apply the block matrix to a block vector.
   *
   * Computes the matrix-vector product <tt>Ax = b</tt>.
   *
   * @param x Input vector.
   * @param b Output vector.
   *
   */
  virtual void multiply( BlockVectorView<VECTOR> const & x,
                         BlockVectorView<VECTOR> & b ) const override;

  //@}
  //! @name Accessors/Setters
  //@{

  /**
   * @brief Get the matrix corresponding to block (<tt>blockRowIndex</tt>,<tt>blockColIndex</tt>).
   */
  OPERATOR & block( localIndex const blockRowIndex,
                    localIndex const blockColIndex ) const;

  /**
   * @brief Set block (<tt>i</tt>,<tt>j</tt>) using <tt>matrix</tt>.
   */
  void set( localIndex const blockRowIndex,
            localIndex const blockColIndex,
            OPERATOR & matrix );

  //@}

private:

  array2d< OPERATOR * > m_matrices;

};

}// end geosx namespace


#endif /*GEOSX_LINEARALGEBRA_UTILITIES_BLOCKOPERATORVIEW_HPP_*/
