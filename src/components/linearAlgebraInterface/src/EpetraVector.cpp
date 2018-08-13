/**
 * @file EpetraVector.cpp
 */

#include "EpetraVector.hpp"

namespace geosx
{
// Create an empty vector
EpetraVector::EpetraVector()
{}

// Create a vector from array
void EpetraVector::create( const globalIndex size, double *V )
{
  Epetra_Map map = Epetra_Map( size, 0, Epetra_MpiComm( MPI_COMM_WORLD ));
  vector = std::unique_ptr<Epetra_Vector>( new Epetra_Vector( View, map, V ));
}

// Create a vector from array
void EpetraVector::create( const Epetra_Map& Map, double *V )
{
  Epetra_Map map = Epetra_Map( Map );
  vector = std::unique_ptr<Epetra_Vector>( new Epetra_Vector( View, map, V ));
}

// Create a vector from vector
void EpetraVector::create( std::vector<double> &vec )
{
  globalIndex m_size = vec.size();
  Epetra_Map map = Epetra_Map( m_size, 0, Epetra_MpiComm( MPI_COMM_WORLD ));
  vector = std::unique_ptr<Epetra_Vector>( new Epetra_Vector( View, map, vec.data()));
}

// Multiply all elements by scalingFactor.
void EpetraVector::scale( real64 const scalingFactor )
{
  vector.get()->Scale( scalingFactor );
}
// Dot product with the vector vec.
void EpetraVector::dot( EpetraVector const &vec,
                        real64 *dst )
{
  vector.get()->Dot( *vec.getPointer(), dst );
}
// Multiply all elements by scalingFactor.
void EpetraVector::update( real64 const alpha,
                           EpetraVector const &vec,
                           real64 const beta )
{
  vector.get()->Update( alpha, *vec.getPointer(), beta );
}

// Return the global size of the vector (total number of elements).
globalIndex EpetraVector::globalSize() const
{
  return vector.get()->GlobalLength64();
}

// Return the local size of the vector (total number of local elements).
localIndex EpetraVector::localSize() const
{
  return vector.get()->MyLength();
}

// Print vector to the terminal in Trilinos format.
void EpetraVector::print() const
{
  std::cout << *vector.get() << std::endl;
}

// Accessors
// Get const pointer
const Epetra_Vector* EpetraVector::getPointer() const
{
  return vector.get();
}

// Get non-const pointer
Epetra_Vector* EpetraVector::getPointer()
{
  return vector.get();
}

}
