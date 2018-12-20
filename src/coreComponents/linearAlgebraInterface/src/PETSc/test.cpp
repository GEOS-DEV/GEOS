// #include "PETScVector.hpp"
// #include "PETScSparseMatrix.hpp"
#include "PETScSolver.hpp"


// to compile: make all
// to run: mpiexec -n 2 ./test

void test_PETScVector(int rank)
{
  // TESTING PETScVector

  // make a PETScVector
  PETScVector vec1;

  // test create from array
  if (rank == 0) printf("create a vector:\n");
  double values1[5] = {2, 2, 3, 1, 4};
  vec1.create(5, values1);

  // test create from vector
  std::vector<double> vec (4, 100);
  printf("size: %lu\n", vec.size());

  PETScVector vec4;
  vec4.create(vec);
  
  // test print
  vec1.print();
  vec4.print();

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

void test_PETScSparseMatrix(int rank)
{
  // TESTING PETScSparseMatrix

  // make vector
  PETScVector vec3;
  double values2[5] = {1, 2, 3, 1, 2};
  vec3.create(5, values2);

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

  // test getPointer
  if (rank == 0) printf("get a pointer:\n");
  printf("%p\n", mat1.getPointer());

  // test getMat
  MatView(mat1.getMat(), PETSC_VIEWER_STDOUT_WORLD);

  return;

}

void PETSc_KSP_example(int rank)
{
  // PETSc example
  Vec x, b, u;
  Mat A;
  KSP ksp;
  double norm;
  int i, j, Ii, J, Istart, Iend, m=3, n=3, its;
  double v;

  // create matrix
  MatCreate(PETSC_COMM_WORLD, &A);
  MatSetType(A, MATMPIAIJ);
  MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, m*n, m*n);
  MatMPIAIJSetPreallocation(A, 5, NULL, 5, NULL);
  MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  MatSetUp(A);

  // populate matrix
  MatGetOwnershipRange(A, &Istart, &Iend);

  for (Ii=Istart; Ii<Iend; Ii++) 
  {
    v = -1.0; i = Ii/n; j = Ii - i*n;
    if (i>0)   {J = Ii - n; MatSetValues(A,1,&Ii,1,&J,&v,ADD_VALUES);}
    if (i<m-1) {J = Ii + n; MatSetValues(A,1,&Ii,1,&J,&v,ADD_VALUES);}
    if (j>0)   {J = Ii - 1; MatSetValues(A,1,&Ii,1,&J,&v,ADD_VALUES);}
    if (j<n-1) {J = Ii + 1; MatSetValues(A,1,&Ii,1,&J,&v,ADD_VALUES);}
    v = 4.0; MatSetValues(A,1,&Ii,1,&Ii,&v,ADD_VALUES);
  }

  // assemble matrix
  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

  // create vectors
  VecCreate(PETSC_COMM_WORLD, &u);
  VecSetSizes(u, PETSC_DECIDE, m*n);
  VecSetFromOptions(u);
  VecDuplicate(u, &b);
  VecDuplicate(b, &x);

  VecSet(u, 1.0);
  MatMult(A, u, b);

  // create linear solver
  KSPCreate(PETSC_COMM_WORLD, &ksp);
  KSPSetOperators(ksp, A, A);

  // solve the system
  KSPSolve(ksp, b, x);

  // check error
  VecAXPY(x, -1.0, u);
  VecNorm(x, NORM_2, &norm);
  KSPGetIterationNumber(ksp, &its);

  // print info
  MatView(A, PETSC_VIEWER_STDOUT_WORLD);
  VecView(u, PETSC_VIEWER_STDOUT_WORLD);
  VecView(x, PETSC_VIEWER_STDOUT_WORLD);
  VecView(b, PETSC_VIEWER_STDOUT_WORLD);

  printf("rank: %d owns %d - %d rows\n", rank, Istart, Iend);
  if (rank == 0) printf("Norm of error: %f, number of iterations: %d\n", norm, its);

  return;
}

void test_KSP_solver(int rank)
{

  // solve a problem with PETScVector and PETScSparseMatrix

  // set up variables
  PETScVector x, b, u;
  PETScSparseMatrix A;
  KSP ksp;
  double norm;
  int i, j, Istart, Iend, m=3, n=3, its;

  // set up matrix
  A.create(PETSC_COMM_WORLD, m*n, 5);

  // populate matrix
  // need a process to fill its own rows
  for(i = 0; i < m*n; i++){
    for(j = 0; j < m*n; j++){
      if(i == j) {A.set(i, j, 4.0);}
      if(i == j + 1 && i < m*n) {A.set(i, j, -1.0);}
      if(i == j - 1 && j < m*n) {A.set(i, j, -1.0);}
      if(i == j + 3 && i < m*n) {A.set(i, j, -1.0);}
      if(i == j - 3 && j < m*n) {A.set(i, j, -1.0);}
    }
  }
  
  // set vectors
  double ones[m*n];
  double zeros[m*n];
  for(i = 0; i < m*n; i++){
    ones[i] = 1.0;
    zeros[i] = 0.0;
  }
  u.create(m*n, ones);

  // set x and b to something
  x.create(m*n, zeros);
  b.create(m*n, zeros);

  A.multiply(u, b);
  
  // create linear solver and solve system
  KSPCreate(PETSC_COMM_WORLD, &ksp);
  KSPSetOperators(ksp, A.getMat(), A.getMat());
  KSPSolve(ksp, b.getVec(), x.getVec());

  // check error
  x.update(1.0, u, -1.0);
  x.norm2(norm);
  KSPGetIterationNumber(ksp, &its);

  // print info
  A.print();
  u.print();
  x.print(); // results worse than PETSc?
  b.print();
  
  Istart = A.ilower(); Iend = A.iupper();

  printf("rank: %d owns %d - %d rows\n", rank, Istart, Iend);
  if (rank == 0) printf("Norm of error: %f, number of iterations: %d\n", norm, its);

  return;
}

void test_PETScSolver(int rank)
{

  // TESTING PETScSolver

  // set up variables
  PETScSolver solver;
  PETScVector x, b, u;
  PETScSparseMatrix A;
  double norm;
  int i, j, m=3, n=3, its;

  // set up matrix
  A.create(PETSC_COMM_WORLD, m*n, 5);

  // populate matrix
  // need a process to fill its own rows
  for(i = 0; i < m*n; i++){
    for(j = 0; j < m*n; j++){
      if(i == j) {A.set(i, j, 4.0);}
      if(i == j + 1 && i < m*n) {A.set(i, j, -1.0);}
      if(i == j - 1 && j < m*n) {A.set(i, j, -1.0);}
      if(i == j + 3 && i < m*n) {A.set(i, j, -1.0);}
      if(i == j - 3 && j < m*n) {A.set(i, j, -1.0);}
    }
  }
  
  // set vectors
  double ones[m*n];
  double zeros[m*n];
  for(i = 0; i < m*n; i++){
    ones[i] = 1.0;
    zeros[i] = 0.0;
  }
  u.create(m*n, ones);

  // set x and b to something
  x.create(m*n, zeros);
  b.create(m*n, zeros);

  x.print();

  A.multiply(u, b);

  solver.solve(A, b, x, 100, 0.000001);
  // solver.dsolve(A, b, x);
  x.print();

  // check error
  x.update(1.0, u, -1.0);
  x.norm2(norm);

  // print info
  A.print();
  u.print();
  x.print(); // results worse than PETSc?
  b.print();
  if (rank == 0) printf("Norm of error: %f\n", norm);

  return;
}

int main()
{
  PetscInitializeNoArguments();

  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  // run tests
  // test_PETScVector(rank);
  test_PETScSparseMatrix(rank);
  // PETSc_KSP_example(rank);
  // test_KSP_solver(rank);
  // test_PETScSolver(rank);



  

  PetscFinalize();
  return 0;

}