#include "PetscVector.hpp"
#include <time.h>

PetscVector::PetscVector()
{
	// Do nothing
}

PetscVector::PetscVector(PetscVector const & vec)
{
	VecDuplicate(vec._vec, &_vec); 
	VecCopy(vec._vec, _vec); 
}

PetscVector::PetscVector(Vec vec)
{
	_vec = vec;
}

// PetscVector::~PetscVector()
// {
//   if (_vec) VecDestroy(&_vec);
// }

void PetscVector::create(const int size, double *vals)
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

void PetscVector::create(std::vector<double> &vec)
{
  int size = vec.size();

  int indices[size];
  for (int i = 0; i < size; i++){
    indices[i] = i;
  }

  VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, size, &_vec);
  VecSetValues(_vec, size, indices, vec.data(), INSERT_VALUES);
  VecAssemblyBegin(_vec);
  VecAssemblyEnd(_vec);
}

/* set value of vector element */
void PetscVector::set(int element, double value)
{
	VecSetValue(_vec, element, value, INSERT_VALUES);
	VecAssemblyBegin(_vec);
  	VecAssemblyEnd(_vec);
}

/* add value to vector element */
void PetscVector::add(int element, double value)
{
	VecSetValue(_vec, element, value, ADD_VALUES);
	VecAssemblyBegin(_vec);
  	VecAssemblyEnd(_vec);
}

// /* add values to vector elements */
// void PetscVector::set(array1d<int> elements, array1d<double> values){
	
// }

void PetscVector::scale(double const scalingFactor)
{
	VecScale(_vec, scalingFactor);
}

void PetscVector::dot(PetscVector const vec, double *dst)
{
	VecDot(_vec, vec._vec, dst);
}

/* x Vector to add ? */
void PetscVector::copy(PetscVector const &x)
{
	VecSet(_vec, 0);
	VecAXPY(_vec, 1.0, x._vec);
}

/* alpha*x + this */
void PetscVector::axpy(double const alpha, PetscVector const &x)
{
	VecAXPY(_vec, alpha, x._vec);
}

void PetscVector::axpby(double const alpha, PetscVector &x, double const beta)
{
	VecScale(_vec, beta);
	VecAXPY(_vec, alpha, x._vec);
}

void PetscVector::norm1(double &result) const
{
	VecNorm(_vec, NORM_1, &result);
}

void PetscVector::norm2(double &result) const
{
	VecNorm(_vec, NORM_2, &result);
}

void PetscVector::normInf(double &result) const
{
	VecNorm(_vec, NORM_INFINITY, &result);
}

int PetscVector::globalSize() const
{
	int size;
	VecGetSize(_vec, &size);
	return size;
}

int PetscVector::localSize() const
{
	int size;
	VecGetLocalSize(_vec, &size);
	return size;
}

double PetscVector::getElement(int i) const
{

	double value[1];
	int index[1] = {i};
	VecGetValues(_vec, 1, index, value);
  	return value[0];
}

const Vec* PetscVector::getPointer() const
{
	return &(_vec);
}

const Vec* PetscVector::getPointer()
{
	return &(_vec);
}

Vec PetscVector::getConstVec() const
{
	return _vec;
}

Vec PetscVector::getVec()
{
	return _vec;
}

void PetscVector::print() const
{
	VecView(_vec, PETSC_VIEWER_STDOUT_WORLD);
}

void PetscVector::zero()
{
	VecZeroEntries(_vec);
}

void PetscVector::rand()
{
	PetscRandom ran;
	PetscRandomCreate(PETSC_COMM_WORLD, &ran);

	time_t seconds;
  	seconds = time (NULL);
	PetscRandomSetSeed(ran, seconds);
	PetscRandomSeed(ran);

	VecSetRandom(_vec, ran);
	PetscRandomDestroy(&ran);
}

void PetscVector::createWithLocalSize(int const localSize, MPI_Comm const & comm)
{
  VecCreateMPI(comm, localSize, PETSC_DETERMINE, &_vec);
}

void PetscVector::createWithGlobalSize(int const globalSize, MPI_Comm const & comm)
{
	VecCreateMPI(comm, PETSC_DECIDE, globalSize, &_vec);	
}

void PetscVector::set(int const *globalIndices, double const *values, int size)
{
	VecSetValues(_vec, size, globalIndices, values, INSERT_VALUES);
	VecAssemblyBegin(_vec);
  	VecAssemblyEnd(_vec);
}  

void PetscVector::add(int const *globalIndices, double const *values, int size)
{
	VecSetValues(_vec, size, globalIndices, values, ADD_VALUES);
	VecAssemblyBegin(_vec);
  	VecAssemblyEnd(_vec);
}  

void PetscVector::set(double const value)
{
	VecSet(_vec, value);
	VecAssemblyBegin(_vec);
  	VecAssemblyEnd(_vec);
}

double PetscVector::get(int globalRow) const
{
	double value[1];
	int index[1] = {globalRow};
	VecGetValues(_vec, 1, index, value);
  	return value[0];
}

void PetscVector::write(std::string const & filename)
{
	// need a char[] for PETSc
	char filename_char[filename.length() + 1]; 
    strcpy(filename_char, filename.c_str()); 

	// set up PETSc viewer
	PetscViewer viewer;
	PetscViewerCreate(PETSC_COMM_WORLD, &viewer);
	PetscViewerSetType(viewer, PETSCVIEWERBINARY);
	PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
	PetscViewerFileSetMode(viewer, FILE_MODE_WRITE);
	PetscViewerFileSetName(viewer, filename_char);
	VecView(_vec, viewer);

	printf("Hannah to do: write()\n");
}

void PetscVector::open()
{
	// do nothing
}

void PetscVector::close()
{
	VecAssemblyBegin(_vec);
  	VecAssemblyEnd(_vec);	
}