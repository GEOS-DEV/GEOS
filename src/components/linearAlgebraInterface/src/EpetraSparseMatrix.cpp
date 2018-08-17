/**
 * @file TrilinosSparseMatrix.cpp
 */


#include "EpetraSparseMatrix.hpp"

namespace geosx
{
// -----------------------------------------------------------------------------------
// ------------------------------- Constructors --------------------------------------
// -----------------------------------------------------------------------------------
// Create an empty matrix
EpetraSparseMatrix::EpetraSparseMatrix()
{}

// -----------------------------------------------------------------------------------
// ------------------------------ Create/Finalize ------------------------------------
// -----------------------------------------------------------------------------------
// Allocate matrix (prepare to be filled with data).

// Create a matrix from number of elements
void EpetraSparseMatrix::create( MPI_Comm const comm,
                                 globalIndex const in_m_nRowGlobal,
                                 integer const     nMaxEntriesPerRow )
{
  Epetra_Map map = Epetra_Map( in_m_nRowGlobal, 0, Epetra_MpiComm( comm ));
  matrix = std::unique_ptr<Epetra_CrsMatrix>( new Epetra_CrsMatrix( Copy, map, nMaxEntriesPerRow, false ));
}

// Create a matrix from number of elements
void EpetraSparseMatrix::create( MPI_Comm const comm,
                                 globalIndex const in_m_nRowGlobal,
                                 globalIndex const in_m_nColGlobal,
                                 integer const nMaxEntriesPerRow )
{
  Epetra_Map rowMap = Epetra_Map( in_m_nRowGlobal, 0, Epetra_MpiComm( comm ));
  Epetra_Map colMap = Epetra_Map( in_m_nColGlobal, 0, Epetra_MpiComm( comm ));
  matrix = std::unique_ptr<Epetra_CrsMatrix>( new Epetra_CrsMatrix( Copy, rowMap, colMap, nMaxEntriesPerRow, false ));
}

//// Create a matrix from number of elements
//void EpetraSparseMatrix::create( MPI_Comm const comm,
//                                 integer const in_m_nRowGlobal,
//                                 std::vector<integer> const nMaxEntriesPerRow )
//{
//  Epetra_Map map = Epetra_Map( in_m_nRowGlobal, 0, Epetra_MpiComm( comm ));
//  matrix = std::unique_ptr<Epetra_CrsMatrix>( new Epetra_CrsMatrix( Copy, map, (int *)const_cast<int *>(&(nMaxEntriesPerRow[0])), false ));
//}
//
//// Create a matrix from number of elements
//void EpetraSparseMatrix::create( MPI_Comm const comm,
//                                 globalIndex const in_m_nRowGlobal,
//                                 globalIndex const in_m_nColGlobal,
//                                 std::vector<integer> const nMaxEntriesPerRow )
//{
//  Epetra_Map rowMap = Epetra_Map( in_m_nRowGlobal, 0, Epetra_MpiComm( comm ));
//  Epetra_Map colMap = Epetra_Map( in_m_nColGlobal, 0, Epetra_MpiComm( comm ));
//  matrix = std::unique_ptr<Epetra_CrsMatrix>( new Epetra_CrsMatrix( Copy, rowMap, colMap, (int *)const_cast<int *>(&(nMaxEntriesPerRow[0])), false ));
//}

// Create a matrix from Epetra_Map
void EpetraSparseMatrix::create( Epetra_Map const &input_map,
                                 integer const nMaxEntriesPerRow )
{
  matrix = std::unique_ptr<Epetra_CrsMatrix>( new Epetra_CrsMatrix( Copy, input_map, nMaxEntriesPerRow, false ));
}

// Create a matrix from two Epetra_Maps
void EpetraSparseMatrix::create( Epetra_Map const &input_row_map,
                                 Epetra_Map const &input_col_map,
                                 integer const nMaxEntriesPerRow )
{
  matrix = std::unique_ptr<Epetra_CrsMatrix>( new Epetra_CrsMatrix( Copy, input_row_map, input_col_map, nMaxEntriesPerRow, false ));
}

// Create a matrix from an Epetra_CrsMatrix.
void EpetraSparseMatrix::create( Epetra_CrsMatrix &in_matrix )
{
  matrix = std::unique_ptr<Epetra_CrsMatrix>( &in_matrix );
}

// Reinitialize. Keeps the map and graph but sets all values to 0.
void EpetraSparseMatrix::zero()
{
  matrix->PutScalar( 0 );
}

// Empty open function (implemented fo HYPRE compatibility).
void EpetraSparseMatrix::open()
{}

// Assemble the matrix when filled
void EpetraSparseMatrix::close()
{
  matrix->FillComplete();
  assembled = true;
}

// -----------------------------------------------------------------------------------
// ---------------------------------- Add/Set ----------------------------------------
// -----------------------------------------------------------------------------------

// Add single value at row iRow and column iCol
void EpetraSparseMatrix::add( globalIndex const iRow,
                              globalIndex const iCol,
                              real64 const value )
{
  matrix->SumIntoGlobalValues( iRow, 1, &value, &iCol );
}

// Add values at row iRow and columns cols (size nCols)
void EpetraSparseMatrix::add( globalIndex const iRow,
                              integer const nCols,
                              real64 const *values,
                              globalIndex const *cols )
{

#if 1
  matrix->SumIntoGlobalValues( iRow, nCols, values, cols );
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
      if( Graph_.FindGlobalIndexLoc( locRow, Index, j, Loc ))
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
void EpetraSparseMatrix::set( globalIndex const iRow,
                              globalIndex const iCol,
                              real64 const value )
{
  matrix->ReplaceGlobalValues( iRow, 1, &value, &iCol );
}

// Set values at row iRow and columns cols (size nCols)
void EpetraSparseMatrix::set( globalIndex const iRow,
                              integer const nCols,
                              real64 const *values,
                              globalIndex const *cols )
{
  matrix->ReplaceGlobalValues( iRow, nCols, values, cols );
}

// Set values at row iRow and columns cols (size nCols)
void EpetraSparseMatrix::insert( globalIndex const iRow,
                                 integer const nCols,
                                 real64 const *values,
                                 globalIndex const *cols )
{
  matrix->InsertGlobalValues( iRow, nCols, values, cols );
}

// -----------------------------------------------------------------------------------
// ------------------------------ Linear Algebra -------------------------------------
// -----------------------------------------------------------------------------------

// Matrix/vector multiplication with src. Result sent to dst.
void EpetraSparseMatrix::multiply( EpetraVector const &src,
                                   EpetraVector &dst )
{
  matrix->Multiply( false, *src.getPointer(), *dst.getPointer());
}

// Compute res = b - Ax (residual form).
void EpetraSparseMatrix::residual( EpetraVector const &x,
                                   EpetraVector const &b,
                                   EpetraVector &res )
{
  matrix->Multiply( false, *x.getPointer(), *res.getPointer());
  res.update( -1.0, b, 1.0 );
}

// Compute gaxpy r = alpha*A*x + beta*b.
void EpetraSparseMatrix::gaxpy( real64 alpha,
                                EpetraVector const &x,
                                real64 beta,
                                EpetraVector const &b,
                                EpetraVector &res,
                                bool useTranspose)
{
  matrix->Multiply( useTranspose, *x.getPointer(), *res.getPointer());
  res.update(alpha,b,beta);
}

// Multiply all elements by scalingFactor.
void EpetraSparseMatrix::scale( real64 scalingFactor )
{
  matrix->Scale( scalingFactor );
}

// Pre-multiplies (left) with diagonal matrix consisting of the values in vec.
void EpetraSparseMatrix::leftScale( EpetraVector const &vec )
{
  matrix->LeftScale( *vec.getPointer());
}

// Post-multiplies (right) with diagonal matrix consisting of the values in vec.
void EpetraSparseMatrix::rightScale( EpetraVector const &vec )
{
  matrix->RightScale( *vec.getPointer());
}

// Pre-multiplies (left) with diagonal matrix consisting of the values in vecLeft and
// Post-multiplies (right) with diagonal matrix consisting of the values in vecRight.
void EpetraSparseMatrix::leftRightScale( EpetraVector const &vecLeft,
                                         EpetraVector const &vecRight )
{
  matrix->RightScale(*vecRight.getPointer());
  matrix->LeftScale(*vecLeft.getPointer());
}

void EpetraSparseMatrix::getRow( int GlobalRow,
                                 int &NumEntries,
                                 std::vector<real64> &vecValues,
                                 std::vector<int>    &vecIndices )
{
  real64* Values;
  int* Indices;
  matrix->ExtractGlobalRowView( GlobalRow, NumEntries, Values, Indices );
  if( matrix->MyGRID( GlobalRow ))
  {
    vecIndices.assign( Indices, Indices+NumEntries );
    vecValues.assign( Values, Values+NumEntries );
  }
}

void EpetraSparseMatrix::getLocalRow( localIndex GlobalRow,
                                      integer &NumEntries,
                                      std::vector<real64> &vecValues,
                                      std::vector<localIndex>    &vecIndices )
{
  real64* Values;
  localIndex* Indices;
  matrix->ExtractMyRowView( GlobalRow, NumEntries, Values, Indices );
  vecIndices.assign( Indices, Indices+NumEntries );
  vecValues.assign( Values, Values+NumEntries );
}

// -----------------------------------------------------------------------------------
// --------------------------------- Accessors ---------------------------------------
// -----------------------------------------------------------------------------------

// Accessor for the pointer to the matrix
Epetra_CrsMatrix * EpetraSparseMatrix::getPointer() const
{
  return matrix.get();
}

// Accessor for the number of global rows
globalIndex EpetraSparseMatrix::globalRows() const
{
  return matrix->RowMap().NumGlobalElements64();
}

// Accessor for the number of global columns
globalIndex EpetraSparseMatrix::globalCols() const
{
  return matrix->ColMap().NumGlobalElements64();
}

// Accessor for the number of global columns
globalIndex EpetraSparseMatrix::uniqueCols() const
{
  return matrix->DomainMap().NumGlobalElements64();
}

// Accessor for the number of global rows
Epetra_Map const & EpetraSparseMatrix::RowMap() const
{
  return matrix->RowMap();
}

// Accessor for the number of global columns
Epetra_Map const & EpetraSparseMatrix::ColMap() const
{
  return matrix->ColMap();
}

// Accessor for the number of global columns
Epetra_Map const & EpetraSparseMatrix::DomainMap() const
{
  return matrix->DomainMap();
}

// Accessor for the number of local rows
int EpetraSparseMatrix::myRows() const
{
  return matrix->NumMyRows();
}

// Accessor for the number of local columns
int EpetraSparseMatrix::myCols() const
{
  return matrix->NumMyCols();
}

// Wrapper to print the trilinos output of the matrix
void EpetraSparseMatrix::print() const
{
  std::cout << *matrix.get() << std::endl;
}

// Returns the infinity norm of the matrix.
real64 EpetraSparseMatrix::normInf() const
{
  return matrix->NormInf();
}

// Returns the one norm of the matrix.
real64 EpetraSparseMatrix::norm1() const
{
  return matrix->NormOne();
}

// Returns the Frobenius norm of the matrix.
real64 EpetraSparseMatrix::normFrobenius() const
{
  return matrix->NormFrobenius();
}

// Boolean indicator. True = matrix assembled and ready to be used.
bool EpetraSparseMatrix::isAssembled() const
{
  return assembled;
}

}
