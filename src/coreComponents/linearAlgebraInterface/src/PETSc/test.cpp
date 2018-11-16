// #include "PETScSparseMatrix.hpp"
#include "PETScSolver.hpp"
// #include "../hannah-PETScVector/PETScVector.hpp"

// to compile: make all
// to run: mpiexec -n 2 ./test

int main()
{
  PetscInitializeNoArguments();

  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  // TESTING PETScVector

  // make a PETScVector
  PETScVector vec1;

  // test create from array
  if (rank == 0) printf("create a vector:\n");
  double values1[5] = {2, 2, 3, 1, 4};
  vec1.create(5, values1);

  // test print
  vec1.print();

  // test copy
  if (rank == 0) printf("copy a vector:\n");
  PETScVector vec2;
  vec2 = PETScVector(vec1);
  vec2.print();

  // test scale
  if (rank == 0) printf("scale a vector:\n");
  vec1.scale(.5);
  vec1.print();

  // test dot
  PETScVector vec3;
  double values2[5] = {1, 2, 3, 1, 2};
  vec3.create(5, values2);

  double dotproduct;
  vec1.dot(vec3, &dotproduct);
  if(rank == 0) printf("dot product is: %f\n", dotproduct);
  
  // test update
  if (rank == 0) printf("update a vector:\n");
  vec1.update(2, vec3, 3);
  vec1.print();

  // test getVec
  if (rank == 0) printf("get a vector:\n");
  VecView(vec1.getVec(), PETSC_VIEWER_STDOUT_WORLD);

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

  // TESTING PETScSparseMatrix

  // test create
  if (rank == 0) printf("create a square matrix:\n");
  PETScSparseMatrix mat1;
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
  PETScSparseMatrix mat2(mat1);
  mat2.print();

  // test create
  if (rank == 0) printf("create a rectangular matrix:\n");
  PETScSparseMatrix mat3;
  // 4 rows, 3 columns with 3 nonzeros per row
  mat3.create(PETSC_COMM_WORLD, 5, 3, 3);
  mat3.print();

  // test create
  if (rank == 0) printf("copy a matrix:\n");
  PETScSparseMatrix mat4;
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
  PETScVector vec4(vec3);
  mat1.multiply(vec3, vec4);
  vec4.print();

  // make new vectors and stuff
  PETScVector vec5;
  double values5[3] = {1, 0, 2};
  vec5.create(3, values5);
  vec5.print();

  PETScVector vec6;
  double values6[3] = {2, 4.5, 2};
  vec6.create(3, values6);
  vec6.print();

  PETScVector vec7;
  vec7 = PETScVector(vec6);
  vec7.print();

  PETScSparseMatrix mat5;
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
  PETScSparseMatrix mat6(mat5);
  mat6.leftScale(vec5);
  mat6.print();

  // test rightScale
  if (rank == 0) printf("right scale a matrix:\n");
  PETScSparseMatrix mat7(mat5);
  mat7.rightScale(vec5);
  mat7.print();

  // test leftRightScale
  if (rank == 0) printf("left and right scale a matrix:\n");
  PETScSparseMatrix mat8(mat5);
  mat8.leftRightScale(vec5, vec7);
  mat8.print();

  // test gaxpy
  if (rank == 0) printf("compute gaxpy:\n");
  PETScVector vec8(vec7);
  mat6.gaxpy(.5, vec5, 2, vec7, vec8, false);
  vec8.print();

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

  // test getMat
  MatView(mat1.getMat(), PETSC_VIEWER_STDOUT_WORLD);

  // TESTING PETScSolver

  PETScSolver solver;
  PETScVector vec9;
  double values9[3] = {2, 12.5, 9};
  vec9.create(3, values9);
  PETScVector vec10(vec9);
  // solver.solve(mat7, vec9, vec10, 100, 10);
  solver.dsolve(mat7, vec9, vec10);
  mat7.print();
  vec10.print();

  PetscFinalize();
  return 0;




}