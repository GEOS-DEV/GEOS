// #include "PetscVector.hpp"
#include "PetscSparseMatrix.hpp"

// to compile: make all (check makefile)
// to run: mpiexec -n 2 ./testPetscOperations

void test_PetscVector(int rank)
{
  // TESTING PetscVector

  // make a PetscVector
  PetscVector vec1;

  // test create from array
  if (rank == 0) printf("create a vector:\n");
  double values1[5] = {2, 2, 3, 1, 4};
  vec1.create(5, values1);

  // test create from vector
  std::vector<double> vec (4, 100);
  printf("size: %lu\n", vec.size());

  PetscVector vec4;
  vec4.create(vec);
  
  // test print
  vec1.print();
  vec4.print();

  // test copy
  if (rank == 0) printf("copy a vector:\n");
  PetscVector vec2;
  vec2 = PetscVector(vec1);
  vec2.print();

  // test set
  if (rank == 0) printf("set value in a vector:\n");
  vec2.set(1, 10);
  vec2.print();

  // test add
  if (rank == 0) printf("add to value in a vector:\n");
  vec2.add(1, 2);
  vec2.print();

  // test scale
  if (rank == 0) printf("scale a vector:\n");
  vec1.scale(.5);
  vec1.print();

  // test dot
  PetscVector vec3;
  double values2[5] = {1, 2, 3, 1, 2};
  vec3.create(5, values2);

  // test copy
  if (rank == 0) printf("copy a vector:\n");
  vec2.copy(vec3);
  vec2.print();

  // test dot
  double dotproduct;
  vec1.dot(vec3, &dotproduct);
  if(rank == 0) printf("dot product is: %f\n", dotproduct);

  // test axpby
  if (rank == 0) printf("axpy a vector:\n");
  vec1.print();
  vec1.axpy(2, vec2);
  vec1.print();
  
  // test axpby
  if (rank == 0) printf("axpby a vector:\n");
  vec1.axpby(2, vec3, 3);
  vec1.print();

  // test getVec
  if (rank == 0) printf("get a vector:\n");
  VecView(vec1.getVec(), PETSC_VIEWER_STDOUT_WORLD);

  // test getPointer
  if (rank == 0) printf("get a pointer:\n");
  printf("%p\n", vec1.getPointer());

  // test norms
  double norm1, norm2, normInf;
  vec3.norm1(norm1);
  vec3.norm2(norm2);
  vec3.normInf(normInf);

  if(rank == 0){
    printf("1-norm is: %f\n", norm1);
    printf("2-norm is: %f\n", norm2);
    printf("infinity-norm is: %f\n", normInf);
  }
  
  // test globalSize, localsize
  if(rank == 0) printf("global size is: %d\n", vec1.globalSize());
  printf("local size of vec1 is: %d\n", vec1.localSize());

  return;
}

void test_PetscSparseMatrix(int rank)
{
  // TESTING PetscSparseMatrix
  printf("\n\n\n");

  // make vector
  PetscVector vec3;
  double values2[5] = {1, 2, 3, 1, 2};
  vec3.create(5, values2);

  // test create
  if (rank == 0) printf("create a square matrix:\n");
  PetscSparseMatrix mat1;
  // 5 rows with 3 nonzeros per row
  mat1.create(PETSC_COMM_WORLD, 5, 3);
  mat1.print();

  // test add
  if (rank == 0) printf("add value:\n");
  mat1.add(3, 3, 1);
  mat1.print();

  // test set
  if (rank == 0) printf("set value:\n");
  mat1.add(4, 3, -1);
  mat1.print();

  // test constructor
  if (rank == 0) printf("copy a matrix:\n");
  PetscSparseMatrix mat2(mat1);
  mat2.print();

  // test create
  if (rank == 0) printf("create a rectangular matrix:\n");
  PetscSparseMatrix mat3;
  // 4 rows, 3 columns with 3 nonzeros per row
  mat3.create(PETSC_COMM_WORLD, 5, 3, 3);
  mat3.print();

  // test create
  if (rank == 0) printf("copy a matrix:\n");
  PetscSparseMatrix mat4;
  mat4.create(mat3);
  mat4.print();

  // test zero
  mat4.set(0, 1, .5);
  mat4.set(2, 2, 3);
  mat4.print();
  mat4.zero();
  if (rank == 0) printf("zero a matrix:\n");
  mat4.print();

  // test set
  if (rank == 0) printf("set values to a matrix:\n");
  if (rank == 0) printf("before:\n");
  mat1.print();

  int row = 2;
  const int ncols = 2;
  int cols[ncols] = {0, 2};
  double values[ncols] = {3, -1};

  mat1.set(row, ncols, values, cols);
  if (rank == 0) printf("after:\n");
  mat1.print();

  // test add
  if (rank == 0) printf("add values to a matrix:\n");
  if (rank == 0) printf("before:\n");
  mat1.print();

  int row2 = 3;
  const int ncols2 = 3;
  int cols2[3] = {0, 3, 4};
  double values3[3] = {1, .5, -.1};

  mat1.add(row2, ncols2, values3, cols2);
  if (rank == 0) printf("after:\n");
  mat1.print();

  // test multiply
  if (rank == 0) printf("multiply a matrix and vector:\n");
  PetscVector vec4(vec3);
  mat1.multiply(vec3, vec4);
  vec4.print();

  // make new vectors and stuff
  PetscVector vec5;
  double values5[3] = {1, 0, 2};
  vec5.create(3, values5);
  vec5.print();

  PetscVector vec6;
  double values6[3] = {2, 4.5, 2};
  vec6.create(3, values6);
  vec6.print();

  PetscVector vec7;
  vec7 = PetscVector(vec6);
  vec7.print();

  PetscSparseMatrix mat5;
  mat5.create(PETSC_COMM_WORLD, 3, 3, 3);
  int cols_[3] = {0, 1, 2};
  double row1[3] = {2, 1, 0};
  double row2_[3] = {.5, -1, 2};
  double row3[3] = {3, 2, 1};
  mat5.set(0, 3, row1, cols_);
  mat5.set(1, 3, row2_, cols_);
  mat5.set(2, 3, row3, cols_);
  mat5.print();

  // test residual
  if (rank == 0) printf("compute residual:\n");
  mat5.residual(vec5, vec6, vec7);
  vec7.print();

  // test scale
  if (rank == 0) printf("scale a matrix:\n");
  mat2.scale(.5);
  mat2.print();

  // test leftScale
  if (rank == 0) printf("left scale a matrix:\n");
  PetscSparseMatrix mat6(mat5);
  mat6.leftScale(vec5);
  mat6.print();

  // test rightScale
  if (rank == 0) printf("right scale a matrix:\n");
  PetscSparseMatrix mat7(mat5);
  mat7.rightScale(vec5);
  mat7.print();

  // test leftRightScale
  if (rank == 0) printf("left and right scale a matrix:\n");
  PetscSparseMatrix mat8(mat5);
  mat8.leftRightScale(vec5, vec7);
  mat8.print();

  // test gemv
  if (rank == 0) printf("compute gemv:\n");
  mat6.gemv(.5, vec5, 2, vec7, false);
  vec7.print();

  // test clearRow
  if (rank == 0) printf("clear a matrix row:\n");
  mat1.clearRow(3, 2);
  mat1.print();

  // test globalRows, globalCols
  if(rank == 0){
    printf("global rows: %d\n", mat4.globalRows());
    printf("global columns: %d\n", mat4.globalCols());
  }

  // test ilower and iupper
  mat4.print();
  printf("rank %d first row: %d\n", rank, mat4.ilower());
  printf("rank %d last row: %d\n", rank, mat4.iupper());

  // test myRows and myCols
  printf("rank %d number of rows: %d\n", rank, mat4.myRows());
  printf("rank %d number of columns: %d\n", rank, mat4.myCols());

  // test rowMapGID
  printf("rank %d LID of GID 2: %d\n", rank, mat4.rowMapLID(2));
  printf("rank %d LID of GID 4: %d\n", rank, mat4.rowMapLID(4));


  // test norms
  double matnorm1, matnormInf, matnormFrob;
  matnorm1 = mat1.norm1();
  matnormInf = mat1.normInf();
  matnormFrob = mat1.normFrobenius();

  mat1.print();
  if(rank == 0){
    printf("1-norm of mat1 is: %f\n", matnorm1);
    printf("infinity-norm of mat1 is: %f\n", matnormInf);
    printf("Frobenius of mat1 is: %f\n", matnormFrob);
  }

  // test getPointer
  if (rank == 0) printf("get a pointer:\n");
  printf("%p\n", mat1.getPointer());

  // test getMat
  MatView(mat1.getMat(), PETSC_VIEWER_STDOUT_WORLD);

  // test getRow
  int row_ = 2;
  int numEntries_;
  int* indices_;
  double* values_;
  if(mat1.ilower() <= row_ && row_ <= mat1.iupper())
  {
    mat1.getRow(row_, numEntries_, values_, indices_);
    printf("rank %d number of entires in row %d: %d\n", rank, row_, numEntries_);
  }

  // test getLocalRow
  row_ = 1;
  const double* values__;
  const int* indices__;
  mat1.print();
  mat1.getLocalRow(row_, numEntries_, values__, indices__);
  printf("rank %d number of entires in row %d: %d\n", rank, row_, numEntries_);
  if (rank == 0) {
    row_ = 2;
    mat1.getLocalRow(row_, numEntries_, values__, indices__);
    printf("rank %d number of entires in row %d: %d\n", rank, row_, numEntries_);
  }
  
  return;

}

int main()
{
  PetscInitializeNoArguments();

  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  // run tests
  test_PetscVector(rank);
  test_PetscSparseMatrix(rank);

  PetscFinalize();
  return 0;

}