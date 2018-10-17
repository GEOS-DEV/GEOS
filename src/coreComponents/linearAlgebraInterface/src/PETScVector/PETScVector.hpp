# include <petscvec.h>

class PETScVector
{

 public:

  /* Constructor */
  PETScVector();

  /* Copy constructor */
  PETScVector(PETScVector& vec);

  /* Destructor */
  // virtual ~PETScVector();

  /* Create from array */
  void create(const int size, double *values);

  /* Create from Epetra_Map and array of values */
  // void create(const Epetra_Map &map, real64 *values);

  /* Create from std::vector */
  //void create(std::vector<real64> &vec);

  /* Multiply all elements by scalingFactor */
  void scale(double const scalingFactor);

  /* Doc product with the vector vec */
  void dot(PETScVector vec, double *dst);

  /* Update vector as this = alpha*vec + beta*this */
  void update(double const alpha, PETScVector vec, double const beta);

  /* 1-norm of the vector */
  void norm1(double &result) const;

  /* 2-norm of the vector */
  void norm2(double &result) const;

  /* Infinity-norm of the vector */
  void normInf(double &result) const;

  /* Returns the global size of the vector */
  double globalSize() const;

  /* Returns the local size of the vector */
  int localSize() const;

  /* Returns a const pointer to the underlying Petsc Vec */
  // const Vec* getPointer() const;

  /* Returns a non-const pointer to the underlying Petsc Vec */
  // Vec* getPointer();

  /* Print the vector in Petsc format to the terminal */
  void print() const;

 protected:
  
  // Underlying Petsc Vec type
  Vec m_vector;
};


