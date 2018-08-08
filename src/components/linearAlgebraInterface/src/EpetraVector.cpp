/**
 * @file EpetraVector.cpp
 */

#include "EpetraVector.hpp"

namespace geosx
{
// Create an empty vector
EpetraVector::EpetraVector()
{
  Epetra_Map map = Epetra_Map(0,0,Epetra_MpiComm(MPI_COMM_WORLD));
  vector = new Epetra_Vector(View, map, 0);
}

// Create a vector from array
void EpetraVector::create(const globalIndex size, double *V)
{
  Epetra_Map map = Epetra_Map(size,0,Epetra_MpiComm(MPI_COMM_WORLD));
  vector = new Epetra_Vector(View,map,V);
}

// Create a vector from array
void EpetraVector::create(const Epetra_Map& Map, double *V)
{
  Epetra_Map map = Epetra_Map(Map);
  vector = new Epetra_Vector(View,map,V);
}

// Create a vector from vector
void EpetraVector::create(std::vector<double> &vec)
{
  globalIndex m_size = vec.size();
  Epetra_Map map = Epetra_Map(m_size,0,Epetra_MpiComm(MPI_COMM_WORLD));
  vector = new Epetra_Vector(View, map, vec.data());
}


globalIndex EpetraVector::globalSize()
{
  return vector->GlobalLength64();
}

localIndex EpetraVector::localSize()
{
  return vector->MyLength();
}

void EpetraVector::print()
{
  std::cout << *vector << std::endl;
}

const Epetra_Vector* EpetraVector::getPointer() const
{
  return vector;
}

Epetra_Vector* EpetraVector::getPointer()
{
  return vector;
}

}

