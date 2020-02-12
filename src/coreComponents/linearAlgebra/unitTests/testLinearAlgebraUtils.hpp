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
 * @file testLinearAlgebraUtils.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_UNITTESTS_TESTLINEARALGEBRAUTILS_HPP
#define GEOSX_LINEARALGEBRA_UNITTESTS_TESTLINEARALGEBRAUTILS_HPP

#include "common/DataTypes.hpp"
#include "mpiCommunications/MpiWrapper.hpp"

// BEGIN_RST_NARRATIVE testLAOperations.rst

// ==============================
// Compute Identity
// ==============================
// This function computes the identity matrix. It can be used to generate a dummy
// preconditioner.
template<typename LAI>
void computeIdentity( MPI_Comm comm,
                      geosx::globalIndex n,
                      typename LAI::ParallelMatrix &I )
{
  // total dofs = n^2
  geosx::globalIndex N = n * n;

  // Create a matrix of size N with 1 non-zero per row
  I.createWithGlobalSize( N, 1, comm );

  // Loop over rows to fill the matrix
  for( geosx::globalIndex i = I.ilower() ; i < I.iupper() ; i++ )
  {
    // Set the value for element (i,i) to 1
    I.insert( i, i, 1.0 );
  }

  // Close the matrix (make data contiguous in memory)
  I.close();
}

/**
 * @brief Compute the 2D Laplace operator
 *
 * \param comm MPI communicator.
 * \param n size of the nxn mesh for the square 2D Laplace operator matrix. Matrix size will be N=n^2.
 */

// ==============================
// Compute 2D Laplace Operator
// ==============================
// This function computes the matrix corresponding to a 2D Laplace operator. These
// matrices arise from a classical finite volume formulation on a cartesian mesh
// (5-point stencil).  Input is the mesh size, n, from which the total dofs is N = n^2;
template<typename LAI>
void compute2DLaplaceOperator( MPI_Comm comm,
                               geosx::globalIndex n,
			        		             typename LAI::ParallelMatrix &laplace2D )
{
  // total dofs = n^2
  geosx::globalIndex N = n * n;

  // Create a matrix of global size N with 5 non-zeros per row
  laplace2D.createWithGlobalSize( N, 5, comm );

  // Allocate arrays to fill the matrix (values and columns)
  geosx::real64 values[5];
  geosx::globalIndex cols[5];

  // Loop over rows to fill the matrix
  for( geosx::globalIndex i = laplace2D.ilower() ; i < laplace2D.iupper() ; i++ )
  {
    // Re-set the number of non-zeros for row i to 0.
    geosx::localIndex nnz = 0;

    // The left -n: position i-n
    if( i - n >= 0 )
    {
      cols[nnz] = i - n;
      values[nnz] = -1.0;
      nnz++;
    }

    // The left -1: position i-1
    if( i - 1 >= 0 )
    {
      cols[nnz] = i - 1;
      values[nnz] = -1.0;
      nnz++;
    }

    // Set the diagonal: position i
    cols[nnz] = i;
    values[nnz] = 4.0;
    nnz++;

    // The right +1: position i+1
    if( i + 1 < N )
    {
      cols[nnz] = i + 1;
      values[nnz] = -1.0;
      nnz++;
    }

    // The right +n: position i+n
    if( i + n < N )
    {
      cols[nnz] = i + n;
      values[nnz] = -1.0;
      nnz++;
    }

    // Set the values for row i
    laplace2D.insert( i, cols, values, nnz );
  }

  // Close the matrix (make data contiguous in memory)
  laplace2D.close();

}

//@}

#endif //GEOSX_LINEARALGEBRA_UNITTESTS_TESTLINEARALGEBRAUTILS_HPP
