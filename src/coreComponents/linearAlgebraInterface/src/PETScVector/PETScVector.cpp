#include "PETScVector.hpp"

#define CHECK_ERROR(NAME) do { if (ierr != 0) petsc_error(ierr, __FILE__, NAME); } while(0)

PETScVector::PETScVector()
{

}

PETScVector::PETScVector(PETScVector& vec)
{
	VecDuplicate(vec.m_vector, &m_vector); 
	VecCopy(vec.m_vector, m_vector); 
}

// PETScVector::~PETScVector()
// {
//   if (m_vector)
//     VecDestroy(&m_vector);
// }

void PETScVector::create(const int size, double *values)
{
  Vec x;

  int indices[size];
  for (int i = 0; i < size; i++){
    indices[i] = i;
  }

  VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, size, &x);
  VecSetValues(x, size, indices, values, INSERT_VALUES);
  VecAssemblyBegin(x);
  VecAssemblyEnd(x);

  m_vector = x;

}

void PETScVector::scale(double const scalingFactor)
{
	VecScale(m_vector, scalingFactor);
}

void PETScVector::dot(PETScVector vec, double *dst)
{
	VecDot(m_vector, vec.m_vector, dst);
}

void PETScVector::update(double const alpha, PETScVector vec, double const beta)
{
	VecScale(m_vector, beta);
	VecAXPY(m_vector, alpha, vec.m_vector);
}

void PETScVector::norm1(double &result) const
{
	VecNorm(m_vector, NORM_1, &result);
}

void PETScVector::norm2(double &result) const
{
	VecNorm(m_vector, NORM_2, &result);
}

void PETScVector::normInf(double &result) const
{
	VecNorm(m_vector, NORM_INFINITY, &result);
}

double PETScVector::globalSize() const
{
	int size;
	VecGetSize(m_vector, &size);
	return size;
}

int PETScVector::localSize() const
{
	int size;
	VecGetLocalSize(m_vector, &size);
	return size;
}

void PETScVector::print() const
{
	VecView(m_vector, PETSC_VIEWER_STDOUT_WORLD);
}
