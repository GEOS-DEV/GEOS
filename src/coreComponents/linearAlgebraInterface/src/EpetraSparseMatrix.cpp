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

/**
 * @file TrilinosSparseMatrix.cpp
 */


#include "EpetraSparseMatrix.hpp"

namespace geosx
{
// -----------------------------------------------------------------------------------
// ------------------------------- Constructors --------------------------------------
// -----------------------------------------------------------------------------------
// Create an empty matrix (meant to be used for declaration)
EpetraSparseMatrix::EpetraSparseMatrix()
{}

// Copy constructor
EpetraSparseMatrix::EpetraSparseMatrix( EpetraSparseMatrix const &in_mat )
{
  // Check if the vector to be copied is not empty
  if( in_mat.getPointer() != nullptr )
  {
    // Create a unique pointer to an Epetra_CrsMatrix. The data from the input matrix is
    // copied to a new memory location.
    m_matrix = std::unique_ptr<Epetra_CrsMatrix>( new Epetra_CrsMatrix( *in_mat.getPointer() ) );
  }
}

// -----------------------------------------------------------------------------------
// ------------------------------ Create/Finalize ------------------------------------
// -----------------------------------------------------------------------------------
// Allocate matrix (prepare to be filled with data).

// Create a matrix from number of elements
void EpetraSparseMatrix::create( MPI_Comm const comm,
                                 trilinosTypes::gid const size,
                                 trilinosTypes::lid const nMaxEntriesPerRow )
{
  // Create an Epetra_Map of size size.
  Epetra_Map map = Epetra_Map( size, 0, Epetra_MpiComm( comm ) );
  // Create a unique pointer to an Epetra_CrsMatrix defined from the Epetra_Map.
  // Specifies the maximum number of non-zeros per row.
  m_matrix = std::unique_ptr<Epetra_CrsMatrix>( new Epetra_CrsMatrix( Copy, map, nMaxEntriesPerRow, false ) );
}

// Create a rectangular matrix from number of elements (TODO this has not been tested
// and from previous experience it may not be the correct way of getting a rectangular
// matrix).
void EpetraSparseMatrix::create( MPI_Comm const comm,
                                 trilinosTypes::gid const in_m_nRowGlobal,
                                 trilinosTypes::gid const in_m_nColGlobal,
                                 trilinosTypes::lid const nMaxEntriesPerRow )
{
  // Create an Epetra_Map of size in_m_nRowGlobal.
  Epetra_Map rowMap = Epetra_Map( in_m_nRowGlobal, 0, Epetra_MpiComm( comm ));
  // Create an Epetra_Map of size in_m_nColGlobal.
  Epetra_Map colMap = Epetra_Map( in_m_nColGlobal, 0, Epetra_MpiComm( comm ));
  // Create a unique pointer to an Epetra_CrsMatrix defined from the two Epetra_Maps.
  // Specifies the maximum number of non-zeros per row.
  m_matrix = std::unique_ptr<Epetra_CrsMatrix>( new Epetra_CrsMatrix( Copy, rowMap, colMap, nMaxEntriesPerRow, false ) );
}

// Create a matrix from Epetra_Map
void EpetraSparseMatrix::create( Epetra_Map const &input_map,
                                 trilinosTypes::lid const nMaxEntriesPerRow )
{
  // Create a unique pointer to an Epetra_CrsMatrix defined from the input Epetra_Map.
  // Specifies the maximum number of non-zeros per row.
  m_matrix = std::unique_ptr<Epetra_CrsMatrix>( new Epetra_CrsMatrix( Copy, input_map, nMaxEntriesPerRow, false ) );
}

// Create a matrix from two Epetra_Maps
void EpetraSparseMatrix::create( Epetra_Map const &input_row_map,
                                 Epetra_Map const &input_col_map,
                                 trilinosTypes::lid const nMaxEntriesPerRow )
{
  // Create a unique pointer to an Epetra_CrsMatrix defined from the two input Epetra_Map.
  // (TODO this has not been tested and from previous experience it may not be the correct \
  // way of getting a rectangular matrix).
  // Specifies the maximum number of non-zeros per row.
  m_matrix = std::unique_ptr<Epetra_CrsMatrix>( new Epetra_CrsMatrix( Copy, input_row_map, input_col_map, nMaxEntriesPerRow, false ) );
}

// Create a matrix from an Epetra_CrsGraph.
void EpetraSparseMatrix::create( Epetra_CrsGraph const &graph )
{
  // Create a unique pointer to an Epetra_CrsMatrix defined from the input Epetra_CrsGraph.
  // An Epetra_CrsGraph can be used as a sparsity pattern.
  m_matrix = std::unique_ptr<Epetra_CrsMatrix>( new Epetra_CrsMatrix( Copy, graph ) );
}

// Create a matrix from an Epetra_CrsMatrix.
void EpetraSparseMatrix::create( Epetra_CrsMatrix &in_matrix )
{
  // Create a unique pointer to an Epetra_CrsMatrix defined from the input Epetra_CrsMatrix.
  // Similar to the copy constructor but does no copy the data, just points to the input matrix.
  m_matrix = std::unique_ptr<Epetra_CrsMatrix>( &in_matrix );
}

// Reinitialize. Keeps the map and graph but sets all values to 0.
void EpetraSparseMatrix::zero()
{
  // Set all existing values of the matrix to 0.
  m_matrix->PutScalar( 0 );
}

// Empty open function (implemented for HYPRE compatibility).
void EpetraSparseMatrix::open()
{}

// Assemble the matrix when filled
void EpetraSparseMatrix::close()
{
  // Fix the sparsity pattern, make the data contiguous in memory and optimize storage.
  m_matrix->FillComplete();
  // Switch boolean to true when done assembling.
  assembled = true;
}

// -----------------------------------------------------------------------------------
// ---------------------------------- Add/Set ----------------------------------------
// -----------------------------------------------------------------------------------

// Add single value at row iRow and column iCol
void EpetraSparseMatrix::add( trilinosTypes::gid const iRow,
                              trilinosTypes::gid const iCol,
                              real64 const value )
{
  // Add the value to the element (iRow,iCol).
  m_matrix->SumIntoGlobalValues( iRow, 1, &value, &iCol );
}

// Add values at row iRow and columns cols (size nCols)
void EpetraSparseMatrix::add( trilinosTypes::gid const iRow,
                              trilinosTypes::lid const nCols,
                              real64 const *values,
                              trilinosTypes::gid const *cols )
{

#if 1
  // Add the values to the elements (iRow,[nCols]).
  m_matrix->SumIntoGlobalValues( iRow, nCols, values, cols );
#else
//  template<typename int_type>
//  int Epetra_CrsMatrix::TSumIntoGlobalValues(int_type Row,
//              int NumEntries,
//              const double * srcValues,
//              const int_type *Indices)
//  {
  int j;
  int ierr = 0;
  int Loc = 0;


  int locRow = Graph_.LRID( Row ); // Normalize row range

  if( locRow < 0 || locRow >= NumMyRows_ )
  {
    EPETRA_CHK_ERR( -1 ); // Not in Row range
  }

  if( StaticGraph() && !Graph_.HaveColMap())
  {
    EPETRA_CHK_ERR( -1 );
  }

  double * RowValues = Values( locRow );

  if( !StaticGraph())
  {
    for( j=0 ; j<NumEntries ; j++ )
    {
      int_type Index = Indices[j];
      if( Graph_.FindtrilinosTypes::gidLoc( locRow, Index, j, Loc ))
//  #ifdef EPETRA_HAVE_OMP
//  #ifdef EPETRA_HAVE_OMP_NONASSOCIATIVE
//  #pragma omp atomic
//  #endif
//  #endif
//          RowValues[Loc] += srcValues[j];
        RAJA::atomic::atomicAdd<ATOMIC_POL2>( &RowValues [Loc], srcValues[j] );
      else
        ierr = 2;   // Value Excluded
    }
  }
  else
  {
    const Epetra_BlockMap& colmap = Graph_.ColMap();
    int NumColIndices = Graph_.NumMyIndices( locRow );
    const int* ColIndices = Graph_.Indices( locRow );

    if( Graph_.Sorted())
    {
      int insertPoint;
      for( j=0 ; j<NumEntries ; j++ )
      {
        int Index = colmap.LID( Indices[j] );

        // Check whether the next added element is the subsequent element in
        // the graph indices, then we can skip the binary search
        if( Loc < NumColIndices && Index == ColIndices[Loc] )
//  #ifdef EPETRA_HAVE_OMP
//  #ifdef EPETRA_HAVE_OMP_NONASSOCIATIVE
//  #pragma omp atomic
//  #endif
//  #endif
//            RowValues[Loc] += srcValues[j];
          RAJA::atomic::atomicAdd<ATOMIC_POL2>( &RowValues [Loc], srcValues[j] );
        else
        {
          Loc = Epetra_Util_binary_search( Index, ColIndices, NumColIndices, insertPoint );
          if( Loc > -1 )
//  #ifdef EPETRA_HAVE_OMP
//  #ifdef EPETRA_HAVE_OMP_NONASSOCIATIVE
//  #pragma omp atomic
//  #endif
//  #endif
//              RowValues[Loc] += srcValues[j];
            RAJA::atomic::atomicAdd<ATOMIC_POL2>( &RowValues [Loc], srcValues[j] );
          else
            ierr = 2;   // Value Excluded
        }
        ++Loc;
      }
    }
    else
      for( j=0 ; j<NumEntries ; j++ )
      {
        int Index = colmap.LID( Indices[j] );
        if( Graph_.FindMyIndexLoc( NumColIndices, ColIndices, Index, j, Loc ))
//  #ifdef EPETRA_HAVE_OMP
//  #ifdef EPETRA_HAVE_OMP_NONASSOCIATIVE
//  #pragma omp atomic
//  #endif
//  #endif
//            RowValues[Loc] += srcValues[j];
          RAJA::atomic::atomicAdd<ATOMIC_POL2>( &RowValues [Loc], srcValues[j] );
        else
          ierr = 2;   // Value Excluded
      }
  }

  NormOne_ = -1.0;   // Reset Norm so it will be recomputed.
  NormInf_ = -1.0;   // Reset Norm so it will be recomputed.
  NormFrob_ = -1.0;

  EPETRA_CHK_ERR( ierr );

  return(0);
//  }
#endif
}

// Set single value at row iRow and column iCol
void EpetraSparseMatrix::set( trilinosTypes::gid const iRow,
                              trilinosTypes::gid const iCol,
                              real64 const value )
{
  // Set the value of the element (iRow,iCol) to value.
  m_matrix->ReplaceGlobalValues( iRow, 1, &value, &iCol );
}

// Set values at row iRow and columns cols (size nCols)
void EpetraSparseMatrix::set( trilinosTypes::gid const iRow,
                              trilinosTypes::lid const nCols,
                              real64 const *values,
                              trilinosTypes::gid const *cols )
{
  // Add the value of elements (iRow,[nCols]) to values.
  m_matrix->ReplaceGlobalValues( iRow, nCols, values, cols );
}

// Set values at row iRow and columns cols (size nCols)
// TODO remove the possibility to dynamically construct the sparsity pattern.
void EpetraSparseMatrix::insert( trilinosTypes::gid const iRow,
                                 trilinosTypes::lid const nCols,
                                 real64 const *values,
                                 trilinosTypes::gid const *cols )
{
  // Insert the value to the element (iRow,iCol). This routine should not be used unless
  // there is a guarantee that each element will only be accessed once. If that is not the
  // case, it will allocate (way) too much memory and sum all the values allocated to the
  // same element.
  m_matrix->InsertGlobalValues( iRow, nCols, values, cols );
}

// -----------------------------------------------------------------------------------
// ------------------------------ Linear Algebra -------------------------------------
// -----------------------------------------------------------------------------------

// Matrix/vector multiplication with src. Result sent to dst.
void EpetraSparseMatrix::multiply( EpetraVector const &src,
                                   EpetraVector &dst ) const
{
  // Perform the matrix-vector product A*src = dst.
  m_matrix->Multiply( false, *src.getPointer(), *dst.getPointer() );
}

// Compute res = b - Ax (residual form).
void EpetraSparseMatrix::residual( EpetraVector const &x,
                                   EpetraVector const &b,
                                   EpetraVector &r ) const
{
  // Compute the matrix-vector product Ax = r.
  m_matrix->Multiply( false, *x.getPointer(), *r.getPointer() );
  // Update r as r = b - r.
  r.axpby( 1.0, b, -1.0 );
}

// Compute gemv <tt>y = alpha*A*x + beta*y</tt>.
void EpetraSparseMatrix::gemv( real64 const alpha,
                               EpetraVector const &x,
                               real64 const beta,
                               EpetraVector &y,
                               bool useTranspose )
{
  // Declare Ax = y.
  EpetraVector Ax( y );
  // Compute Ax = A*x.
  m_matrix->Multiply( useTranspose, *x.getPointer(), *Ax.getPointer() );
  // Compute y = beta*y + alpha*Ax
  y.axpby( alpha, Ax, beta );
}

// Multiply all elements by scalingFactor.
void EpetraSparseMatrix::scale( real64 const scalingFactor )
{
  // Scale every element of the matrix.
  m_matrix->Scale( scalingFactor );
}

// Pre-multiplies (left) with diagonal matrix consisting of the values in vec.
void EpetraSparseMatrix::leftScale( EpetraVector const &vec )
{
  // Pre-multiply all elements of the matrix with the input vector.
  m_matrix->LeftScale( *vec.getPointer() );
}

// Post-multiplies (right) with diagonal matrix consisting of the values in vec.
void EpetraSparseMatrix::rightScale( EpetraVector const &vec )
{
  // Post-multiply all elements of the matrix with the input vector.
  m_matrix->RightScale( *vec.getPointer() );
}

// Pre-multiplies (left) with diagonal matrix consisting of the values in vecLeft and
// Post-multiplies (right) with diagonal matrix consisting of the values in vecRight.
void EpetraSparseMatrix::leftRightScale( EpetraVector const &vecLeft,
                                         EpetraVector const &vecRight )
{
  // Post-multiply all elements of the matrix with the rigth input vector.
  m_matrix->RightScale( *vecRight.getPointer() );
  // Pre-multiply all elements of the matrix with the left input vector.
  m_matrix->LeftScale( *vecLeft.getPointer() );
}

// Clear row and multiply diagonal term by factor.
void EpetraSparseMatrix::clearRow( trilinosTypes::gid const row,
                                   real64 const factor )
{
  // Get the local index corresponding to the global index input. Returns -1
  // if the global index is not owned by this processor.
  int local_row = m_matrix->LRID( row );

  // If this processor owns the input row:
  if( local_row >= 0 )
  {
    double *values = nullptr;
    int *col_indices = nullptr;
    int num_entries;

    // Extract the row and get the number of entries, the values and the column indices.
    m_matrix->ExtractMyRowView( local_row, num_entries, values, col_indices );

    // TODO this if may be unnecessary?
    if( values != nullptr && col_indices != nullptr && num_entries > 0 )
    {
      int* diag_find = std::find( col_indices, col_indices+num_entries-1, local_row );
      long int diag_index = (diag_find - col_indices);

      // Set all non-diagonal entries to 0
      for( int j=0 ; j<num_entries ; ++j )
      {
        if( diag_index != j )
        {
          values[j] = 0.;
        }
      }
      // Scale the diagonal (existing) value with factor.
      values[diag_index] *= factor;
    }
  }
}

void EpetraSparseMatrix::getRow( trilinosTypes::gid GlobalRow,
                                 trilinosTypes::lid &NumEntries,
                                 real64* Values,
                                 trilinosTypes::gid* Indices ) const
{
  // Extract the global row and output the number of elements, values and column indices.
  m_matrix->ExtractGlobalRowView( GlobalRow, NumEntries, Values, Indices );
}

void EpetraSparseMatrix::getLocalRow( trilinosTypes::lid localRow,
                                      trilinosTypes::lid & NumEntries,
                                      real64 * & Values,
                                      trilinosTypes::lid * & Indices ) const
{
  // Extract the local row and output the number of elements, values and column indices.
  m_matrix->ExtractMyRowView( localRow, NumEntries, Values, Indices );
}

void EpetraSparseMatrix::getRow( trilinosTypes::gid GlobalRow,
                                 trilinosTypes::lid &NumEntries,
                                 std::vector<real64> &vecValues,
                                 std::vector<trilinosTypes::gid> &vecIndices ) const
{
  // Same as the other getRow function but outputs the values and column indices to
  // standard vectors.
  real64* Values;
  trilinosTypes::gid* Indices;
  m_matrix->ExtractGlobalRowView( GlobalRow, NumEntries, Values, Indices );
  if( m_matrix->MyGRID( GlobalRow ))
  {
    vecIndices.assign( Indices, Indices+NumEntries );
    vecValues.assign( Values, Values+NumEntries );
  }
}

void EpetraSparseMatrix::getLocalRow( trilinosTypes::lid localRow,
                                      trilinosTypes::lid &NumEntries,
                                      std::vector<real64> &vecValues,
                                      std::vector<trilinosTypes::lid> &vecIndices ) const
{
  // Same as the other getLocalRow function but outputs the values and column indices to
  // standard vectors.
  real64* Values;
  trilinosTypes::lid* Indices;
  m_matrix->ExtractMyRowView( localRow, NumEntries, Values, Indices );
  vecIndices.assign( Indices, Indices+NumEntries );
  vecValues.assign( Values, Values+NumEntries );
}

// -----------------------------------------------------------------------------------
// --------------------------------- Accessors ---------------------------------------
// -----------------------------------------------------------------------------------

// Accessor for the pointer to the matrix
Epetra_CrsMatrix * EpetraSparseMatrix::getPointer() const
{
  return m_matrix.get();
}

// Accessor for the number of global rows
trilinosTypes::gid EpetraSparseMatrix::globalRows() const
{
  return m_matrix->RowMap().NumGlobalElements64();
}

// Accessor for the number of global columns
trilinosTypes::gid EpetraSparseMatrix::globalCols() const
{
  return m_matrix->ColMap().NumGlobalElements64();
}

// Accessor for the number of global columns
trilinosTypes::gid EpetraSparseMatrix::uniqueCols() const
{
  return m_matrix->DomainMap().NumGlobalElements64();
}

// Accessor for the index of the first global row
trilinosTypes::gid EpetraSparseMatrix::ilower() const
{
  return m_matrix->RowMap().MyGlobalElements64()[0];
}

// Accessor for the index of the last global row
trilinosTypes::gid EpetraSparseMatrix::iupper() const
{
  return m_matrix->RowMap().MyGlobalElements64()[0] + m_matrix->RowMap().NumMyElements();
}

// Accessor for the number of global rows
Epetra_Map const & EpetraSparseMatrix::RowMap() const
{
  return m_matrix->RowMap();
}

// Accessor for the number of global columns
Epetra_Map const & EpetraSparseMatrix::ColMap() const
{
  return m_matrix->ColMap();
}

// Accessor for the number of global columns
Epetra_Map const & EpetraSparseMatrix::DomainMap() const
{
  return m_matrix->DomainMap();
}

// Accessor for the number of local rows
int EpetraSparseMatrix::myRows() const
{
  return m_matrix->NumMyRows();
}

// Accessor for the number of local columns
int EpetraSparseMatrix::myCols() const
{
  return m_matrix->NumMyCols();
}

// Accessor for the number of local columns
trilinosTypes::lid EpetraSparseMatrix::rowMapLID( trilinosTypes::gid const GID ) const
{
  return m_matrix->RowMap().LID( GID );
}

// Wrapper to print the trilinos output of the matrix
void EpetraSparseMatrix::print() const
{
  if( m_matrix.get() != nullptr )
    std::cout << *m_matrix.get() << std::endl;
}

// Returns the infinity norm of the matrix.
real64 EpetraSparseMatrix::normInf() const
{
  return m_matrix->NormInf();
}

// Returns the one norm of the matrix.
real64 EpetraSparseMatrix::norm1() const
{
  return m_matrix->NormOne();
}

// Returns the Frobenius norm of the matrix.
real64 EpetraSparseMatrix::normFrobenius() const
{
  return m_matrix->NormFrobenius();
}

// Boolean indicator. True = matrix assembled and ready to be used.
bool EpetraSparseMatrix::isAssembled() const
{
  return assembled;
}

}
