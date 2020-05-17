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

using namespace geosx;

static void Q12d_local( real64 const & hx, real64 const & hy,
                        real64 const & E, real64 const & nu,
                        stackArray2d< real64, 8*8 > & Ke );

/*! @name Utility functions.
 * @brief Functions used to construct useful matrices in the test files.
 */
//@{

/**
 * @brief Compute an identity matrix
 *
 * @param comm MPI communicator.
 * @param N global size of the square identity matrix.
 */

// BEGIN_RST_NARRATIVE testLAOperations.rst

// ==============================
// Compute Identity
// ==============================
// This function computes the identity matrix. It can be used to generate a dummy
// preconditioner.
template< typename MATRIX >
void computeIdentity( MPI_Comm comm,
                      globalIndex N,
                      MATRIX & I )
{
  // Create a matrix of size N with 1 non-zero per row
  I.createWithGlobalSize( N, 1, comm );

  I.open();

  // Loop over rows to fill the matrix
  for( globalIndex i = I.ilower(); i < I.iupper(); i++ )
  {
    // Set the value for element (i,i) to 1
    I.insert( i, i, 1.0 );
  }

  // Close the matrix (make data contiguous in memory)
  I.close();
}

template< typename MATRIX >
void computeZero( MPI_Comm comm,
                  globalIndex N,
                  MATRIX & I )
{
  I.createWithGlobalSize( N, 0, comm );
  I.open();
  I.close();
}

/**
 * @brief Compute the 2D Laplace operator
 *
 * @param comm MPI communicator.
 * @param n size of the nxn mesh for the square 2D Laplace operator matrix. Matrix size will be N=n^2.
 */

// ==============================
// Compute 2D Laplace Operator
// ==============================
// This function computes the matrix corresponding to a 2D Laplace operator. These
// matrices arise from a classical finite volume formulation on a cartesian mesh
// (5-point stencil).  Input is the mesh size, n, from which the total dofs is N = n^2;
template< typename MATRIX >
void compute2DLaplaceOperator( MPI_Comm comm,
                               globalIndex n,
                               MATRIX & laplace2D )
{
  // total dofs = n^2
  globalIndex N = n * n;

  // Create a matrix of global size N with 5 non-zeros per row
  laplace2D.createWithGlobalSize( N, 5, comm );

  // Allocate arrays to fill the matrix (values and columns)
  real64 values[5];
  globalIndex cols[5];

  // Open the matrix
  laplace2D.open();

  // Loop over rows to fill the matrix
  for( globalIndex i = laplace2D.ilower(); i < laplace2D.iupper(); i++ )
  {
    // Re-set the number of non-zeros for row i to 0.
    localIndex nnz = 0;

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

// ==============================
// Compute 2D Laplace Operator
// ==============================
// This function computes the matrix corresponding to a 2D elasticity operator,
// assuming plane strain conditions based on finite element discretization.
// A regular grid grid consisting of quadrilateral element. ...
//
// beam 0 disp on the left boundary
//
// INPUT: Lx, Ly, Nx, Ny
//
template< typename MATRIX >
void compute2DElasticityOperator( MPI_Comm comm,
                                  real64 domainSizeX,
                                  real64 domainSizeY,
                                  globalIndex nCellsX,
                                  globalIndex nCellsY,
                                  real64 youngModulus,
                                  real64 poissonRatio,
                                  MATRIX & elasticity2D )
{
  localIndex const rank  = LvArray::integerConversion< localIndex >( MpiWrapper::Comm_rank( comm ) );
  localIndex const nproc = LvArray::integerConversion< localIndex >( MpiWrapper::Comm_size( comm ) );

  // Compute total number of grid nodes (nNodes) and elements (nCells)
  globalIndex nCells = nCellsX * nCellsY;
  GEOSX_ERROR_IF( nCells < nproc, "less than one cell per processor" );
  globalIndex nNodes = ( nCellsX + 1 ) * ( nCellsY + 1);
  real64 hx = domainSizeX / static_cast< real64 >( nCellsX );
  real64 hy = domainSizeY / static_cast< real64 >( nCellsY );

  // Compute cell partitioning
  localIndex nLocalCells = LvArray::integerConversion< localIndex >( nCells / nproc );
  localIndex nExtraCells = LvArray::integerConversion< localIndex >( nCells ) - nLocalCells * nproc;
  globalIndex iCellLower = rank * nLocalCells + ( rank == 0 ? 0 : nExtraCells );
  globalIndex iCellUpper = iCellLower + nLocalCells + ( rank == 0 ? 0 : nExtraCells ) - 1;

  // Construct local stiffness matrix (same for all cells)
//  array2d< real64 > Ke( 8, 8 );
//  Q12d_local( hx, hy, youngModulus, poissonRatio, Ke );

  stackArray2d< real64, 8*8 > Ke( 8, 8 );
  Q12d_local( hx, hy, youngModulus, poissonRatio, Ke );

  // Create a matrix of global size N with at most 18 non-zeros per row
  elasticity2D.createWithGlobalSize( nNodes*2, 18, comm );

  // Open the matrix
  elasticity2D.open();

  // Loop over grid cells
  stackArray1d< globalIndex, 4 > cellNodes( 4 );
  stackArray1d< globalIndex, 8 > localDofIndex( 8 );

  for( localIndex iCell = iCellLower; iCell <= iCellUpper; ++iCell )
  {

    // Compute local DOF global indeces
    cellNodes( 0 ) = (iCell / nCellsX) + iCell;
    cellNodes( 1 ) = cellNodes( 0 ) + 1;
    cellNodes( 3 ) = cellNodes( 1 ) + nCellsX;
    cellNodes( 2 ) = cellNodes( 3 ) + 1;
    for( localIndex i = 0; i < 4; ++i )
    {
      localDofIndex( 2 * i )     = cellNodes( i ) * 2;
      localDofIndex( 2 * i + 1 ) = localDofIndex( 2 * i ) + 1;
    }

    // Assemble local stiffness matrix and right-hand side
    elasticity2D.insert( localDofIndex.data(), localDofIndex.data(), Ke.data(), 8, 8 );
  }

  // Close the matrix
  elasticity2D.close();
}

static void Q12d_local( real64 const & hx, real64 const & hy,
                        real64 const & E, real64 const & nu,
                        stackArray2d< real64, 8*8 > & Ke )
{
  real64 fac = E / ( 1. - 2. * nu ) / (1. + nu );

  // Populate stiffness matrix

  // --- Fill diagonal entries
  real64 Dxx = ( fac * hx * ( 1. - 2. * nu ) ) / ( 6. * hy )
               - ( fac * hy * ( -1. + nu ) ) / ( 3. * hx );
  real64 Dyy = ( fac * hy * ( 1. - 2. * nu ) ) / ( 6. * hx )
               - ( fac * hx * ( -1. + nu ) ) / ( 3. * hy );
  for( localIndex i = 0; i < 8; i += 2 )
  {
    Ke( i, i ) = Dxx;
    Ke( i + 1, i + 1 ) = Dyy;
  }

  // --- Fill upper triangular part
  // --- --- Ke( 0, 1:7 )
  Ke( 0, 1 ) = fac / 8.;
  Ke( 0, 2 ) = ( fac * hx * ( 1. - 2. * nu ) ) / ( 12. * hy )
               + ( fac * hy * ( -1. + nu ) ) / ( 3. * hx );
  Ke( 0, 3 ) = ( fac * ( -1 + 4. * nu ) ) / 8.;
  Ke( 0, 4 ) = ( fac * hy * ( -1. + nu ) ) / ( 6. * hx )
               + ( fac * hx * (-1. + 2. * nu ) ) / ( 12. * hy );
  Ke( 0, 5 ) = -Ke( 0, 1 );
  Ke( 0, 6 ) = -( fac * hy * ( -1. + nu ) ) / ( 6. * hx )
               + ( fac * hx * ( -1. + 2. * nu ) ) / ( 6. * hy );
  Ke( 0, 7 ) = -( fac * ( -1. + 4. * nu ) ) / 8.;

  // --- --- Ke( 1, 2:7 )
  Ke( 1, 2 ) = Ke( 0, 7 );
  Ke( 1, 3 ) = -( fac * ( hy * hy * ( 1. - 2. * nu ) + hx * hx *( -1. + nu ) ) ) / ( 6. * hx * hy );
  Ke( 1, 4 ) = Ke( 0, 5 );
  Ke( 1, 5 ) = ( fac * hx * ( -1. + nu ) ) / ( 6. * hy ) + ( fac * hy * ( -1. + 2. * nu ) ) / ( 12. * hx );
  Ke( 1, 6 ) = Ke( 0, 3 );
  Ke( 1, 7 ) = ( fac * hy * ( 1. - 2. * nu ) ) / ( 12. * hx )
               + ( fac * hx * ( -1. + nu ) ) / ( 3. * hy );

  // --- --- Ke( 2, 3:7 )
  Ke( 2, 3 ) =  Ke( 0, 5 );
  Ke( 2, 4 ) =  Ke( 0, 6 );
  Ke( 2, 5 ) =  Ke( 1, 6 );
  Ke( 2, 6 ) = ( fac * hy * ( -1 + nu ) ) / ( 6. * hx ) + ( fac * hx * ( -1. + 2. * nu ) ) / ( 12. * hy );
  Ke( 2, 7 ) =  Ke( 0, 1 );

  // --- --- Ke( 3, 4:7 )
  Ke( 3, 4 ) = Ke( 1, 2 );
  Ke( 3, 5 ) = Ke( 1, 7 );
  Ke( 3, 6 ) = Ke( 0, 1 );
  Ke( 3, 7 ) = Ke( 1, 5 );

  // --- --- Ke( 4, 5:7 )
  Ke( 4, 5 ) = Ke( 0, 1 );
  Ke( 4, 6 ) = Ke( 0, 2 );
  Ke( 4, 7 ) = Ke( 0, 3 );

  // --- --- Ke( 5, 6:7 )
  Ke( 5, 6 ) = Ke( 0, 7 );
  Ke( 5, 7 ) = Ke( 1, 3 );

  // --- --- Ke( 6, 7 )
  Ke( 6, 7 ) = Ke( 0, 5 );

  // --- Fill lower triangular part
  for( localIndex i = 1; i < 8; ++i )
  {
    for( localIndex j = 0; j < i; ++j )
    {
      Ke( i, j ) = Ke( j, i );
    }
  }
}
// END_RST_NARRATIVE

//@}

#endif //GEOSX_LINEARALGEBRA_UNITTESTS_TESTLINEARALGEBRAUTILS_HPP
