#include "PETScVector.hpp"

PETScVector::PETScVector()
{
	// Do nothing
}

PETScVector::PETScVector(PETScVector& vec)
{
	VecDuplicate(vec._vec, &_vec); 
	VecCopy(vec._vec, _vec); 
}

PETScVector::PETScVector(Vec vec)
{
	_vec = vec;
}

// PETScVector::~PETScVector()
// {
//   if (_vec)
//     VecDestroy(&_vec);
// }

void PETScVector::create(const int size, double *vals)
{

  int indices[size];
  for (int i = 0; i < size; i++){
    indices[i] = i;
  }

  VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, size, &_vec);
  VecSetValues(_vec, size, indices, vals, INSERT_VALUES);
  VecAssemblyBegin(_vec);
  VecAssemblyEnd(_vec);

}

void PETScVector::scale(double const scalingFactor)
{
	VecScale(_vec, scalingFactor);
}

void PETScVector::dot(PETScVector vec, double *dst)
{
	VecDot(_vec, vec._vec, dst);
}

void PETScVector::update(double const alpha, PETScVector vec, double const beta)
{
	VecScale(_vec, beta);
	VecAXPY(_vec, alpha, vec._vec);
}

void PETScVector::norm1(double &result) const
{
	VecNorm(_vec, NORM_1, &result);
}

void PETScVector::norm2(double &result) const
{
	VecNorm(_vec, NORM_2, &result);
}

void PETScVector::normInf(double &result) const
{
	VecNorm(_vec, NORM_INFINITY, &result);
}

int PETScVector::globalSize() const
{
	int size;
	VecGetSize(_vec, &size);
	return size;
}

int PETScVector::localSize() const
{
	int size;
	VecGetLocalSize(_vec, &size);
	return size;
}

Vec PETScVector::getVec()
{
	return _vec;
}

void PETScVector::print() const
{
	VecView(_vec, PETSC_VIEWER_STDOUT_WORLD);
}
