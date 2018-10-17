#include "PETScVector.hpp"

int main()
{
  PetscInitializeNoArguments();

  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  // make a PETScVector
  PETScVector vec1;

  // test create from array
  double values1[5] = {2, 2, 3, 1, 4};
  vec1.create(5, values1);

  // test print
  printf("vec1 is:");
  vec1.print();

  // test copy
  PETScVector vec2;
  vec2 = PETScVector(vec1);
  printf("vec2 is:");
  vec2.print();

  // test scale
  vec1.scale(.5);
  vec1.print();

  // test dot
  PETScVector vec3;
  double values2[5] = {1, 2, 3, 1, 2};
  vec3.create(5, values2);

  double dotproduct;
  vec1.dot(vec3, &dotproduct);
  if(rank == 0){
    printf("dot product is: %f\n", dotproduct);
  }
  
  // test update
  vec1.update(2, vec3, 3);
  vec1.print();

  // test norms
  double norm1, norm2, normInf;
  vec3.norm1(norm1);
  vec3.norm2(norm2);
  vec3.normInf(normInf);

  if(rank == 0){
    printf("1-norm of vec3 is: %f\n", norm1);
    printf("2-norm of vec3 is: %f\n", norm2);
    printf("infinity-norm of vec3 is: %f\n", normInf);
  }
  
  // test globalSize, localsize
  if(rank == 0){
    printf("global size of vec1 is: %f\n", vec1.globalSize());
  }
  printf("local size of vec1 is: %d\n", vec1.localSize());

  PetscFinalize();
  return 0;
}
