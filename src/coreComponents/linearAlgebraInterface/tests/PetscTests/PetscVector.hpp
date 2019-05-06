#include <petscvec.h>
#include <vector>
#include <string>

class PetscVector
{

 public:

  /* Constructor */
  PetscVector();

  /* Copy constructor */
  PetscVector(PetscVector const & vec);

  /* Construct from PETSc vector */
  PetscVector(Vec vec);

  /* Destructor */
  // virtual ~PetscVector() = default;

  /* Create from array */
  void create(int const size, double *values);

  /* Create from std::vector */
  void create(std::vector<double> &vec);

  /* set value of vector element */
  void set(int element, double value);

  // /* set values of vector elements */
  // void set(array1d<int> elements, array1d<double> values);

  /* add value to element */
  void add(int element, double value);

  // /* add values to vector elements */
  // void add(array1d<int> elements, array1d<double> values);
  
  /* Multiply all elements by scalingFactor */
  void scale(double const scalingFactor);

  /* Doc product with the vector vec */
  void dot(PetscVector const vec, double *dst);

  /* x Vector to add ? */
  void copy(PetscVector const &x);

  /* alpha*x + this */
  void axpy(double const alpha, PetscVector const &x);

  /* Update vector as this = alpha*x + beta*this */
  void axpby(double const alpha, PetscVector &x, double const beta);

  /* 1-norm of the vector */
  void norm1(double &result) const;

  /* 2-norm of the vector */
  void norm2(double &result) const;

  /* Infinity-norm of the vector */
  void normInf(double &result) const;

  /* Returns the global size of the vector */
  int globalSize() const;

  /* Returns the local size of the vector */
  int localSize() const;

  /* Returns element of the vector */
  double getElement(int i) const;

  /* Returns a const pointer to the underlying Petsc Vec */
  const Vec* getPointer() const;

  /* Returns a non-const pointer to the underlying Petsc Vec */
  const Vec* getPointer();

  /* Returns vector */
  Vec getConstVec() const;

  /* Returns vector */
  Vec getVec();

  /* Print the vector in Petsc format to the terminal */
  void print() const;

  /* zero entries in vector */
  void zero();

  /* Set vector to random numbers */
  void rand();

  /* create a vector based on the number of local elements */
  void createWithLocalSize(int const localSize, MPI_Comm const & comm = MPI_COMM_WORLD);

  /* create a vector based on the global number of elements */
  void createWithGlobalSize(int const globalSize, MPI_Comm const & comm = MPI_COMM_WORLD);

  /* set vector values at given elements */
  void set(int const *globalIndices, double const *values, int size);

  /* add to vector values at given elements */
  void add(int const *globalIndices, double const *values, int size);

  /* set vector elements to a value */
  void set(double const value);

  /* get the value of an element */
  double get(int globalRow) const;

  /* write vector to Matlab file */
  void write(std::string const & filename);

  /* open vector (not needed for PETSc) */
  void open();

  /* close vector and assemble */
  void close();

 protected:
  
  // Underlying Petsc Vec type
  Vec _vec;
};


