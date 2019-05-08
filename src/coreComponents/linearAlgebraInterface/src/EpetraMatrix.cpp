/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
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
 * @file EpetraMatrix.cpp
 */

// BEGIN_RST_NARRATIVE EpetraMatrix.rst
// ==============================
// Epetra-based Matrix Object
// ==============================
// This class contains the ParallelMatrix wrappers based on Epetra_Crs Objects.
// The class contains a unique pointer to an Epetra_CrsMatrix as well as constructors,
// functions and accessors for Epetra objects.

// Include the corresponding header file.
#include "EpetraMatrix.hpp"

// Put everything under the geosx namespace.
namespace geosx
{

// ----------------------------
// Constructors
// ----------------------------

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Create an empty matrix (meant to be used for declaration)
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
EpetraMatrix::EpetraMatrix()
{}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Copy constructor
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
EpetraMatrix::EpetraMatrix( EpetraMatrix const &src )
{
  GEOS_ERROR_IF( src.unwrappedPointer() == nullptr, "Input matrix looks empty" );
  //TODO GEOS_ERROR_IF( !src.isClosed(), "Input matrix hasn't been properly closed before copy");

  m_matrix  = std::unique_ptr<Epetra_FECrsMatrix>( new Epetra_FECrsMatrix( *src.unwrappedPointer() ) );
  m_src_map = std::unique_ptr<Epetra_Map>( new Epetra_Map( m_matrix->DomainMap()));
  m_dst_map = std::unique_ptr<Epetra_Map>( new Epetra_Map( m_matrix->RangeMap()));
}

// -----------------------------
// Create
// -----------------------------
// Allocate matrix (prepare to be filled with data).

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Create a matrix from an Epetra_FECrsGraph.
// """""""""""""""""""""""""""""""""""""""""""""""
void EpetraMatrix::create( Epetra_FECrsGraph const &graph )
{
  m_matrix  = std::unique_ptr<Epetra_FECrsMatrix>( new Epetra_FECrsMatrix( Copy, graph ) );
  m_src_map = std::unique_ptr<Epetra_Map>( new Epetra_Map( m_matrix->DomainMap()));
  m_dst_map = std::unique_ptr<Epetra_Map>( new Epetra_Map( m_matrix->RangeMap()));
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Create a matrix from number of elements
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
void EpetraMatrix::createWithGlobalSize( globalIndex const globalSize,
                                         localIndex const maxEntriesPerRow,
                                         MPI_Comm const & comm )
{
  createWithGlobalSize( globalSize, globalSize, maxEntriesPerRow, comm ); // just call general version
}

void EpetraMatrix::createWithGlobalSize( globalIndex const globalRows,
                                         globalIndex const globalCols,
                                         localIndex const maxEntriesPerRow,
                                         MPI_Comm const & comm )
{
  m_dst_map = std::unique_ptr<Epetra_Map>( new Epetra_Map( globalRows, 0, Epetra_MpiComm( comm ) ));
  m_src_map = std::unique_ptr<Epetra_Map>( new Epetra_Map( globalCols, 0, Epetra_MpiComm( comm ) ));
  m_matrix = std::unique_ptr<Epetra_FECrsMatrix>( new Epetra_FECrsMatrix( Copy, *m_dst_map, integer_conversion<int, localIndex>( maxEntriesPerRow ), false ) );
}

void EpetraMatrix::createWithLocalSize( localIndex const localSize,
                                        localIndex const maxEntriesPerRow,
                                        MPI_Comm const & comm )
{
  createWithLocalSize( localSize, localSize, maxEntriesPerRow, comm ); // just call general version
}

void EpetraMatrix::createWithLocalSize( localIndex const localRows,
                                        localIndex const localCols,
                                        localIndex const maxEntriesPerRow,
                                        MPI_Comm const & comm )
{
  m_dst_map = std::unique_ptr<Epetra_Map>( new Epetra_Map( integer_conversion<globalIndex>( -1 ), integer_conversion<int, localIndex>( localRows ), 0, Epetra_MpiComm( comm ) ));
  m_src_map = std::unique_ptr<Epetra_Map>( new Epetra_Map( integer_conversion<globalIndex>( -1 ), integer_conversion<int, localIndex>( localCols ), 0, Epetra_MpiComm( comm ) ));
  m_matrix = std::unique_ptr<Epetra_FECrsMatrix>( new Epetra_FECrsMatrix( Copy, *m_dst_map, integer_conversion<int, localIndex>( maxEntriesPerRow ), false ) );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Reinitialize.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Keeps the map and graph but sets all values to 0.
void EpetraMatrix::zero()
{
  m_matrix->PutScalar( 0 );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Open
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Empty open function (implemented for HYPRE compatibility).
void EpetraMatrix::open()
{}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Close
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Fix the sparsity pattern, make the data contiguous in memory and optimize storage.
void EpetraMatrix::close()
{
  m_matrix->GlobalAssemble( *m_src_map, *m_dst_map );
  assembled = true;
}

// -------------------------
// Add/Set
// -------------------------

// 1x1
void EpetraMatrix::add( globalIndex const rowIndex,
                        globalIndex const colIndex,
                        real64 const value )
{
  m_matrix->SumIntoGlobalValues( rowIndex, 1, &value, &colIndex );
}

void EpetraMatrix::set( globalIndex const rowIndex,
                        globalIndex const colIndex,
                        real64 const value )
{
  m_matrix->ReplaceGlobalValues( rowIndex, 1, &value, &colIndex );
}

void EpetraMatrix::insert( globalIndex const rowIndex,
                           globalIndex const colIndex,
                           real64 const value )
{
  m_matrix->InsertGlobalValues( rowIndex, 1, &value, &colIndex );
}

// 1xN c-style
void EpetraMatrix::add( globalIndex const rowIndex,
                        globalIndex const * colIndices,
                        real64 const * values,
                        localIndex size )
{
  m_matrix->SumIntoGlobalValues( rowIndex, integer_conversion<int, localIndex>( size ), values, colIndices );
}

void EpetraMatrix::set( globalIndex const rowIndex,
                        globalIndex const * colIndices,
                        real64 const * values,
                        localIndex size )
{
  m_matrix->ReplaceGlobalValues( rowIndex, integer_conversion<int, localIndex>( size ), values, colIndices );
}

void EpetraMatrix::insert( globalIndex const rowIndex,
                           globalIndex const * colIndices,
                           real64 const * values,
                           localIndex size )
{
  m_matrix->InsertGlobalValues( rowIndex, integer_conversion<int, localIndex>( size ), values, colIndices );
}

// 1xN array1d style
void EpetraMatrix::add( globalIndex const rowIndex,
                        array1d<globalIndex> const &colIndices,
                        array1d<real64> const &values )
{
  m_matrix->SumIntoGlobalValues( rowIndex, integer_conversion<int, localIndex>( colIndices.size()), values.data(), colIndices.data() );
}

void EpetraMatrix::set( globalIndex const rowIndex,
                        array1d<globalIndex> const &colIndices,
                        array1d<real64> const &values )
{
  m_matrix->ReplaceGlobalValues( rowIndex, integer_conversion<int, localIndex>( colIndices.size()), values.data(), colIndices.data() );
}

void EpetraMatrix::insert( globalIndex const rowIndex,
                           array1d<globalIndex> const &colIndices,
                           array1d<real64> const &values )
{
  m_matrix->InsertGlobalValues( rowIndex, integer_conversion<int, localIndex>( colIndices.size()), values.data(), colIndices.data() );
}

// MxN array2d style
void EpetraMatrix::add( array1d<globalIndex> const & rowIndices,
                        array1d<globalIndex> const & colIndices,
                        array2d<real64> const & values )
{
  m_matrix->SumIntoGlobalValues( integer_conversion<int, localIndex>( rowIndices.size()), rowIndices.data(),
                                 integer_conversion<int, localIndex>( colIndices.size()), colIndices.data(),
                                 values.data(), Epetra_FECrsMatrix::ROW_MAJOR );
}

void EpetraMatrix::set( array1d<globalIndex> const & rowIndices,
                        array1d<globalIndex> const & colIndices,
                        array2d<real64> const & values )
{
  m_matrix->ReplaceGlobalValues( integer_conversion<int, localIndex>( rowIndices.size()), rowIndices.data(),
                                 integer_conversion<int, localIndex>( colIndices.size()), colIndices.data(),
                                 values.data(), Epetra_FECrsMatrix::ROW_MAJOR );
}

void EpetraMatrix::insert( array1d<globalIndex> const & rowIndices,
                           array1d<globalIndex> const & colIndices,
                           array2d<real64> const & values )
{
  m_matrix->InsertGlobalValues( integer_conversion<int, localIndex>( rowIndices.size()), rowIndices.data(),
                                integer_conversion<int, localIndex>( colIndices.size()), colIndices.data(),
                                values.data(), Epetra_FECrsMatrix::ROW_MAJOR );
}

// -------------------------
// Linear Algebra
// -------------------------

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Matrix/vector multiplication
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Perform the matrix-vector product A*src = dst.
void EpetraMatrix::multiply( EpetraVector const &src,
                             EpetraVector &dst ) const
{
  m_matrix->Multiply( false, *src.unwrappedPointer(), *dst.unwrappedPointer() );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Matrix/matrix multiplication
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Perform the matrix-matrix product A*src = dst.
void EpetraMatrix::multiply( EpetraMatrix const &src,
                             EpetraMatrix &dst ) const
{
  int err = EpetraExt::MatrixMatrix::Multiply(*m_matrix,
                                              false,//don't use transpose
                                              *src.unwrappedPointer(),
                                              false,//don't use tranpose,
                                              *dst.unwrappedPointer());

  GEOS_ERROR_IF(err != 0,"Error thrown in matrix/matrix multiply routine");
}


// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Compute residual.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Compute res = b - Ax (residual form).
void EpetraMatrix::residual( EpetraVector const &x,
                             EpetraVector const &b,
                             EpetraVector &r ) const
{
  m_matrix->Multiply( false, *x.unwrappedPointer(), *r.unwrappedPointer() );
  r.axpby( 1.0, b, -1.0 );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Generalized matrix/vector product.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Compute gemv <tt>y = alpha*A*x + beta*y</tt>.
void EpetraMatrix::gemv( real64 const alpha,
                         EpetraVector const &x,
                         real64 const beta,
                         EpetraVector &y,
                         bool useTranspose )
{
  EpetraVector Ax( y );
  m_matrix->Multiply( useTranspose, *x.unwrappedPointer(), *Ax.unwrappedPointer() );
  y.axpby( alpha, Ax, beta );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Scale.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Multiply all elements by scalingFactor.
void EpetraMatrix::scale( real64 const scalingFactor )
{
  m_matrix->Scale( scalingFactor );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Left and right scaling
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
void EpetraMatrix::leftScale( EpetraVector const &vec )
{
  m_matrix->LeftScale( *(*vec.unwrappedPointer())( 0 ) );
}

void EpetraMatrix::rightScale( EpetraVector const &vec )
{
  m_matrix->RightScale( *(*vec.unwrappedPointer())( 0 ) );
}

void EpetraMatrix::leftRightScale( EpetraVector const &vecLeft,
                                   EpetraVector const &vecRight )
{
  leftScale( vecLeft );
  rightScale( vecRight );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get global row copy
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// The challenge here is that columns are stored with local, not global,
// indices, so we need to do conversions back and forth
void EpetraMatrix::getRowCopy( globalIndex globalRow,
                               array1d<globalIndex> & colIndices,
                               array1d<real64> & values ) const
{
  int n_entries = m_matrix->NumGlobalEntries( globalRow );

  localIndex length = integer_conversion<localIndex, int>( n_entries );

  values.resize( length );
  colIndices.resize( length );

  array1d<int> local_indices ( length );

  int localRow = m_matrix->LRID( globalRow );
  int err = m_matrix->ExtractMyRowCopy( localRow, n_entries, n_entries, values.data(), local_indices.data() );
  GEOS_ERROR_IF( err!=0,
                 "getRowCopy failed. This often happens if the requested global row is not local to this processor, or if close() hasn't been called." );

  for( localIndex i=0 ; i<length ; ++i )
    colIndices[i] = m_matrix->GCID64( local_indices[i] );
}


// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Clear row
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Clear the row.  By default the diagonal value will be set
// to zero, but the user can pass a desired diagValue.
void EpetraMatrix::clearRow( globalIndex const globalRow,
                             real64 const diagValue )
{
  double *values = nullptr;
  int length;

  int err = m_matrix->ExtractGlobalRowView( globalRow, length, values );
  GEOS_ERROR_IF( err != 0,
                 "getRowView failed. This often happens if the requested global row is not local to this processor, or if close() hasn't been called." );

  for( int j=0 ; j<length ; ++j )
    values[j] = 0.0;

  set( globalRow, globalRow, diagValue );
}

// ---------------------------------------------------------
//  Accessors
// ---------------------------------------------------------

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get pointer.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Accessor for the pointer to the raw Epetra matrix
Epetra_FECrsMatrix * EpetraMatrix::unwrappedPointer() const
{
  return m_matrix.get();
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get number of global rows.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Accessor for the number of global rows
globalIndex EpetraMatrix::globalRows() const
{
  return m_matrix->NumGlobalRows64();
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get number of global columns.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Accessor for the number of global columns
globalIndex EpetraMatrix::globalCols() const
{
  return m_matrix->NumGlobalCols64();
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get the lower index owned by processor.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Accessor for the index of the first global row
globalIndex EpetraMatrix::ilower() const
{
  return m_matrix->RowMap().MyGlobalElements64()[0];
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get the upper index owned by processor.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Accessor for the index of the last global row
globalIndex EpetraMatrix::iupper() const
{
  return m_matrix->RowMap().MyGlobalElements64()[0] + m_matrix->RowMap().NumMyElements();
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Print to terminal.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Wrapper to print the trilinos output of the matrix
void EpetraMatrix::print() const
{
  if( m_matrix.get() != nullptr )
  {
    std::cout << *m_matrix << std::endl;
  }
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Write to matlab-compatible file
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Note: EpetraExt also supports a MatrixMarket format as well
void EpetraMatrix::write( string const & filename ) const
{
  EpetraExt::RowMatrixToMatlabFile( filename.c_str(), *m_matrix );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Inf-norm.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Returns the infinity norm of the matrix.
real64 EpetraMatrix::normInf() const
{
  return m_matrix->NormInf();
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// 1-norm.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Returns the one norm of the matrix.
real64 EpetraMatrix::norm1() const
{
  return m_matrix->NormOne();
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Frobenius-norm.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Returns the Frobenius norm of the matrix.
real64 EpetraMatrix::normFrobenius() const
{
  return m_matrix->NormFrobenius();
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Is-assembled.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Boolean indicator. True = matrix assembled and ready to be used.
bool EpetraMatrix::isAssembled() const
{
  return assembled;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// MatrixMatrixMultiply
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Perform the matrix-matrix product A*src = dst.
void EpetraMatrix::MatrixMatrixMultiply( bool const transA,
                                         EpetraMatrix const &B,
                                         bool const transB,
                                         EpetraMatrix &C,
                                         bool const call_FillComplete ) const
{
  int
  err = EpetraExt::MatrixMatrix::Multiply( *m_matrix,
                                           transA,
                                           *B.unwrappedPointer(),
                                           transB,
                                           *C.unwrappedPointer(),
                                           call_FillComplete );

  GEOS_ERROR_IF( err != 0, "Error thrown in matrix/matrix multiply routine" );

  // Using "call_FillComplete_on_result = false" with rectangular matrices because in this
  // case the function does not work. After the multiplication is performed, call close().
  if( !call_FillComplete )
  {
    C.close();
  }
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// getLocalRowID
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Map a global row index to local row index
localIndex EpetraMatrix::getLocalRowID( globalIndex const index ) const
{
  return m_matrix->LRID( index );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// getGlobalRowID
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Map a local row index to global row index
localIndex EpetraMatrix::getGlobalRowID( localIndex const index ) const
{
  return m_matrix->GRID64( integer_conversion<int>( index ) );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// getGlobalRowID
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Map a local row index to global row index
localIndex EpetraMatrix::getGlobalRowID( globalIndex const index ) const
{
  return m_matrix->GRID64( integer_conversion<int>( index ) );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// numMyCols
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Return the local number of columns on each processor
// NOTE: direct use of NumMyCols() counts also for overlays. To avoid those, DomainMap() is needed
localIndex EpetraMatrix::numMyCols( ) const
{
  return m_matrix->DomainMap().NumMyElements();
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// printParallelMatrix
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Print the given parallel matrix in Matrix Market format (MTX file)
void EpetraMatrix::printParallelMatrix( string const & fileName ) const
{
  // Ensure the ".mtx" extension
  string name( fileName );
  if( fileName.substr( fileName.find_last_of( "." ) + 1 ) != "mtx" )
  {
    name = fileName.substr( 0, fileName.find_last_of( "." ) ) + ".mtx";
  }

  EpetraExt::RowMatrixToMatrixMarketFile( name.c_str(), *m_matrix );
}

} // end geosx namespace

// END_RST_NARRATIVE

/* TODO: We should make a decision about thread safety in another
 * pull request.  Either we make Epetra threadsafe or we move to
 * Tpetra as an alternative.
 */

/* SCRATCH CODE - possible template for threadsafe assembly */
/* DELETE WHEN NO LONGER NEEDED */

/*
   // """""""""""""""""""""""""""""""""""""""""""""""""""""""""""
   // Ignore the following. Thread-safe prototype, commented out.
   // """""""""""""""""""""""""""""""""""""""""""""""""""""""""""
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
      if( Graph_.FindglobalIndexLoc( locRow, Index, j, Loc ))
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

 */
