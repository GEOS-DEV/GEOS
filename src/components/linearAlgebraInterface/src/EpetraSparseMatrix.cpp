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
{
}

// -----------------------------------------------------------------------------------
// ------------------------------ Create/Finalize ------------------------------------
// -----------------------------------------------------------------------------------
// Allocate matrix (prepare to be filled with data).

// Create a matrix from number of elements
void EpetraSparseMatrix::create(const MPI_Comm    comm,
                                const globalIndex in_m_nRowGlobal,
                                const integer     nMaxEntriesPerRow)
{
  Epetra_Map map = Epetra_Map(in_m_nRowGlobal,0,Epetra_MpiComm(comm));
  matrix = std::unique_ptr<Epetra_CrsMatrix>(new Epetra_CrsMatrix(Copy, map, nMaxEntriesPerRow, false));
}

// Create a matrix from number of elements
void EpetraSparseMatrix::create(const MPI_Comm    comm,
                                const globalIndex in_m_nRowGlobal,
                                const globalIndex in_m_nColGlobal,
                                const integer nMaxEntriesPerRow)
{
  Epetra_Map rowMap = Epetra_Map(in_m_nRowGlobal,0,Epetra_MpiComm(comm));
  Epetra_Map colMap = Epetra_Map(in_m_nColGlobal,0,Epetra_MpiComm(comm));
  matrix = std::unique_ptr<Epetra_CrsMatrix>(new Epetra_CrsMatrix(Copy, rowMap, colMap, nMaxEntriesPerRow, false));
}

// Create a matrix from number of elements
void EpetraSparseMatrix::create(const MPI_Comm             comm,
                                const integer              in_m_nRowGlobal,
                                const std::vector<integer> nMaxEntriesPerRow)
{
  Epetra_Map map = Epetra_Map(in_m_nRowGlobal,0,Epetra_MpiComm(comm));
  matrix = std::unique_ptr<Epetra_CrsMatrix>(new Epetra_CrsMatrix(Copy, map, (int *)const_cast<int *>(&(nMaxEntriesPerRow[0])), false));
}

// Create a matrix from number of elements
void EpetraSparseMatrix::create(const MPI_Comm             comm,
                                const globalIndex          in_m_nRowGlobal,
                                const globalIndex          in_m_nColGlobal,
                                const std::vector<integer> nMaxEntriesPerRow)
{
  Epetra_Map rowMap = Epetra_Map(in_m_nRowGlobal,0,Epetra_MpiComm(comm));
  Epetra_Map colMap = Epetra_Map(in_m_nColGlobal,0,Epetra_MpiComm(comm));
  matrix = std::unique_ptr<Epetra_CrsMatrix>(new Epetra_CrsMatrix(Copy, rowMap, colMap, (int *)const_cast<int *>(&(nMaxEntriesPerRow[0])), false));
}

// Create a matrix from Epetra_Map
void EpetraSparseMatrix::create(const Epetra_Map &input_map,
                                const integer     nMaxEntriesPerRow)
{
  matrix = std::unique_ptr<Epetra_CrsMatrix>(new Epetra_CrsMatrix(Copy, input_map, nMaxEntriesPerRow, false));
}

// Create a matrix from two Epetra_Maps
void EpetraSparseMatrix::create(const Epetra_Map &input_row_map,
                                const Epetra_Map &input_col_map,
                                const integer         nMaxEntriesPerRow)
{
  matrix = std::unique_ptr<Epetra_CrsMatrix>(new Epetra_CrsMatrix(Copy, input_row_map, input_col_map, nMaxEntriesPerRow, false));
}

// Create a matrix from an Epetra_CrsMatrix.
void EpetraSparseMatrix::create(Epetra_CrsMatrix &in_matrix)
{
  matrix = std::unique_ptr<Epetra_CrsMatrix>(&in_matrix);
}

// Reinitialize. Keeps the map and graph but sets all values to 0.
void EpetraSparseMatrix::zero()
{
  matrix->PutScalar(0);
}

// Empty open function (implemented fo HYPRE compatibility).
void EpetraSparseMatrix::open()
{
}

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
void EpetraSparseMatrix::add(const globalIndex    iRow,
                             const globalIndex    iCol,
                             const real64 value)
{
  matrix->SumIntoGlobalValues(iRow,1,&value,&iCol);
}

// Add values at row iRow and columns cols (size nCols)
void EpetraSparseMatrix::add(const globalIndex  iRow,
                             const integer      nCols,
                             const real64      *values,
                             const globalIndex *cols)
{
  matrix->SumIntoGlobalValues(iRow,nCols,values,cols);
}

// Set single value at row iRow and column iCol
void EpetraSparseMatrix::set(const globalIndex    iRow,
                             const globalIndex    iCol,
                             const real64         value)
{
  matrix->ReplaceGlobalValues(iRow,1,&value,&iCol);
}

// Set values at row iRow and columns cols (size nCols)
void EpetraSparseMatrix::set(const globalIndex     iRow,
                             const integer         nCols,
                             const real64         *values,
                             const globalIndex    *cols)
{
  matrix->ReplaceGlobalValues(iRow,nCols,values,cols);
}

// -----------------------------------------------------------------------------------
// ------------------------------ Linear Algebra -------------------------------------
// -----------------------------------------------------------------------------------

// Apply operator. Matrix/vector multiplication with src. Result sent to dst.
void EpetraSparseMatrix::apply(      EpetraVector &dst,
                               const EpetraVector &src)
{
  matrix->Multiply(false,*src.getPointer(),*dst.getPointer());
}


void EpetraSparseMatrix::getRow(int 				 GlobalRow,
                                int 				&NumEntries,
                                std::vector<real64> &vecValues,
                                std::vector<int>    &vecIndices)
{
  real64* Values;
  int* Indices;
  matrix->ExtractGlobalRowView(GlobalRow, NumEntries, Values, Indices);
  if (matrix->MyGRID(GlobalRow))
  {
    vecIndices.assign(Indices,Indices+NumEntries);
    vecValues.assign(Values,Values+NumEntries);
  }
}

// -----------------------------------------------------------------------------------
// --------------------------------- Accessors ---------------------------------------
// -----------------------------------------------------------------------------------

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
  std::cout << *matrix << std::endl;
}

// Boolean indicator. True = matrix assembled and ready to be used.
bool EpetraSparseMatrix::isAssembled() const
{
  return assembled;
}

}
