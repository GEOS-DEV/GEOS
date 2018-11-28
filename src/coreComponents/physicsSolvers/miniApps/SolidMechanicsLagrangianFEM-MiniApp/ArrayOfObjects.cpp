#include <iostream>
#include <cstring>

#include "MemoryManager_impl.hpp"
#include "RAJA/RAJA.hpp"
#include "MatrixMath_impl.hpp"
#include "Layout.hpp"
#include "MeshGen_impl.hpp"
#include "ConstitutiveUpdate_impl.hpp"
#include "ShapeFun_impl.hpp"
#include "omp.h"

#include "../../src/SolidMechanicsLagrangianFEMKernels_impl.hpp"

//
//Driver for the GEOSX proxy app
//Assumes an Array of Objects data layout
//

int main(int argc, char* const argv[])
{

  if(argc != 3)
    {
      std::cout<<"usage ./main NoElem Iter"<<std::endl;
      exit(-1);
    }
#if 1
  size_t dataAllocated   = 0.0; 
  localIndex Kx      = atoi(argv[1]); 
  localIndex Niter   =  atoi(argv[2]);

  std::cout<<"GEOSX mini-app: Array of Objects data structures"<<std::endl;
  printParameters();

  real64 dt = 0.125; //step size
  const localIndex nx = Kx+1; //Number of nodes in a cartesian dimension
  const localIndex NoElem = Kx*Kx*Kx; //Total number of elements
  const localIndex numNodes = nx*nx*nx; //Total number of nodes    

  //
  //Generate an element list
  //
#if defined(USE_GEOSX_ARRAY)
  LvArray::Array<localIndex, 1, localIndex> _elementList(NoElem);
  LvArray::ArrayView<localIndex, 1, localIndex>  & elementList = _elementList;
#else
  localIndex * const elementList = memoryManager::allocate<localIndex>(NoElem, dataAllocated);
#endif

  for(localIndex i=0; i<NoElem; ++i) elementList[i] = i;

  //Allocate space for constitutive map
  //
#if defined(USE_GEOSX_ARRAY)
  LvArray::Array<localIndex, 2, localIndex> _constitutiveMap(NoElem, inumQuadraturePoints);
  LvArray::ArrayView<localIndex,2,localIndex> & constitutiveMap = _constitutiveMap;
#else
  localIndex *  constitutiveMap = memoryManager::allocate<localIndex>(inumQuadraturePoints*NoElem, dataAllocated); 
#endif

  //
  //Generate space for an element to node list
  //
#if defined(USE_GEOSX_ARRAY)
  LvArray::Array<localIndex, 2, localIndex> _elemsToNodes(NoElem, inumNodesPerElement);
  LvArray::ArrayView<localIndex,2,localIndex> & elemsToNodes = _elemsToNodes;
#else
  localIndex  *  elemsToNodes = memoryManager::allocate<localIndex>(inumNodesPerElement*NoElem, dataAllocated);
#endif
 
  //
  //Allocate space for a list of vertices, generate a mesh, and populate the constitutive map
  //
#if defined(USE_GEOSX_ARRAY)
  LvArray::Array<real64, 2, localIndex> _VX(numNodes, local_dim);
  LvArray::ArrayView<real64,2, localIndex> & VX = _VX;
#else
  geosxData VX = memoryManager::allocate<real64>(numNodes*local_dim, dataAllocated);
#endif

  meshGen(VX, elemsToNodes,constitutiveMap,Kx);

  //
  //Precompute evaluate the quadrature points at every basis function
  //with respect to the parent basis function
  //

  P_Wrapper P; 
  generateP(P, inumNodesPerElement, inumQuadraturePoints);

  //
  //Allocate space for shape function derivatives and compute their values
  //
#if defined(USE_GEOSX_ARRAY)
  LvArray::Array<localIndex, 4, localIndex> _dNdX(NoElem, inumQuadraturePoints, inumNodesPerElement, local_dim);
  LvArray::ArrayView<localIndex,4, localIndex> & dNdX = _dNdX;
#else
  geosxData dNdX = memoryManager::allocate<real64>(inumNodesPerElement*inumQuadraturePoints*NoElem*local_dim, dataAllocated);
#endif
  
  make_dNdX(dNdX, VX, elemsToNodes, NoElem, inumQuadraturePoints, inumNodesPerElement);


  ///
  //Allocate space for nodal degrees of freedom as an Array of Objects
  ///
#if defined(USE_GEOSX_ARRAY)
  LvArray::Array<real64, 2, localIndex> _u(numNodes, local_dim);
  LvArray::Array<real64, 2, localIndex> _uhat(numNodes, local_dim);
  LvArray::Array<real64, 2, localIndex> _acc(numNodes, local_dim);

  LvArray::ArrayView<real64, 2, localIndex> & u = _u;
  LvArray::ArrayView<real64, 2, localIndex> & uhat = _uhat;
  LvArray::ArrayView<real64, 2, localIndex> & acc = _acc;
  std::memset(u.data(),0, numNodes*local_dim*sizeof(real64));
  std::memset(uhat.data(),0, numNodes*local_dim*sizeof(real64));
  std::memset(acc.data(),0, numNodes*local_dim*sizeof(real64));
#else
  geosxData u = memoryManager::allocate<real64>(numNodes*local_dim, dataAllocated);
  geosxData uhat = memoryManager::allocate<real64>(numNodes*local_dim, dataAllocated);
  geosxData acc = memoryManager::allocate<real64>(numNodes*local_dim, dataAllocated);
  std::memset(u,0, numNodes*local_dim*sizeof(real64));
  std::memset(uhat,0, numNodes*local_dim*sizeof(real64));
  std::memset(acc,0, numNodes*local_dim*sizeof(real64));
#endif

  //
  //In the case of three kernel launches, we need to allocate
  //extra memory for intermediate results
  //
#if defined(USE_GEOSX_ARRAY) && defined(THREE_KERNEL_UPDATE)
  LvArray::Array<real64, 3, localIndex> _Dadt(NoElem, inumQuadraturePoints, localMatSz);
  LvArray::Array<real64, 3, localIndex> _Rot(NoElem, inumQuadraturePoints, localMatSz);
  LvArray::Array<real64, 3, localIndex> _detF(NoElem, inumQuadraturePoints);
  LvArray::Array<real64, 3, localIndex> _inverseF(NoElem, inumQuadraturePoints, localMatSz);
  
  LvArray::ArrayView<real64, 3, localIndex> & Dadt     = _Dadt;
  LvArray::ArrayView<real64, 3, localIndex> & Rot      = _Rot;
  LvArray::ArrayView<real64, 3, localIndex> & detF     = _detF;
  LvArray::ArrayView<real64, 3, localIndex> & inverseF = _inverseF;

  std::memset(Dadt.data(), 0, localMatSz*inumQuadraturePoints*NoElem*sizeof(real64));
  std::memset(Rot.data(), 0, localMatSz*inumQuadraturePoints*NoElem*sizeof(real64));
  std::memset(detF.data(), 0, inumQuadraturePoints*NoElem*sizeof(real64));
  std::memset(inverseF.data(), 0, localMatSz*inumQuadraturePoints*NoElem*sizeof(real64));

#elif !defined(USE_GEOSX_ARRAY) && defined(THREE_KERNEL_UPDATE)
  geosxData Dadt = memoryManager::allocate<real64>(localMatSz*inumQuadraturePoints*NoElem, dataAllocated);
  geosxData Rot  = memoryManager::allocate<real64>(localMatSz*inumQuadraturePoints*NoElem, dataAllocated);
  geosxData detF  = memoryManager::allocate<real64>(inumQuadraturePoints*NoElem, dataAllocated);
  geosxData inverseF  = memoryManager::allocate<real64>(localMatSz*inumQuadraturePoints*NoElem, dataAllocated);
  
  std::memset(Dadt, 0, localMatSz*inumQuadraturePoints*NoElem*sizeof(real64));
  std::memset(Rot, 0, localMatSz*inumQuadraturePoints*NoElem*sizeof(real64));
  std::memset(detF, 0, inumQuadraturePoints*NoElem*sizeof(real64));
  std::memset(inverseF, 0, localMatSz*inumQuadraturePoints*NoElem*sizeof(real64));
#endif

  //
  //Allocate and set physical parameters
  //
  localIndex noSymEnt = 6;
  real64 bulkModulus = 10, shearModulus=20;

#if defined(USE_GEOSX_ARRAY)
  
  LvArray::Array<real64,2,localIndex> _detJ(NoElem, inumQuadraturePoints);
  LvArray::Array<real64,3,localIndex> _devStressData(NoElem, inumQuadraturePoints, noSymEnt);
  LvArray::Array<real64,1,localIndex> _meanStress(NoElem*inumQuadraturePoints); //Reformulate to 1D

  LvArray::ArrayView<real64,2,localIndex> & detJ          = _detJ;
  LvArray::ArrayView<real64,3,localIndex> & devStressData = _devStressData;  
  LvArray::ArrayView<real64,1,localIndex> & meanStress    = _meanStress;

  std::memset(detJ.data(), 1 ,inumQuadraturePoints * NoElem * sizeof(real64));
  std::memset(meanStress.data(), 1 ,inumQuadraturePoints * NoElem * sizeof(real64));
  std::memset(devStressData.data(), 1 , noSymEnt * inumQuadraturePoints * NoElem * sizeof(real64));
#else
  geosxData detJ            = memoryManager::allocate<real64>(inumQuadraturePoints*NoElem, dataAllocated);
  geosxData meanStress      = memoryManager::allocate<real64>(inumQuadraturePoints*NoElem, dataAllocated);
  geosxData devStressData   = memoryManager::allocate<real64>(noSymEnt*inumQuadraturePoints*NoElem, dataAllocated);

  std::memset(detJ, 1 ,inumQuadraturePoints * NoElem * sizeof(real64));
  std::memset(meanStress, 1 ,inumQuadraturePoints * NoElem * sizeof(real64));
  std::memset(devStressData, 1 , noSymEnt * inumQuadraturePoints * NoElem * sizeof(real64));
#endif
  

#if 1//DELETE
  //
  //Set up function pointer for constitutive relationship
  //
#if defined(USE_CUDA)
  std::cout<<"Using CUDA"<<std::endl;
  constUpdate myUpdate;
  cudaMemcpyFromSymbol(&myUpdate,deviceUpdate,sizeof(deviceUpdate));
#else
  constUpdate myUpdate = UpdateStatePoint;
#endif  

  //-----set up timer-----
  double start, end, diff;
  double myMin = 100000000;

  std::cout<<"Allocated "<<dataAllocated*1e-9<<" GBs of data"<<std::endl;  
  std::cout<<"Launching kernel. . . "<<std::endl;  
  //-----------------[Kernel Launch]---------------------------------    
  for(localIndex it=0; it<Niter; ++it)
    {
      start = omp_get_wtime();


#if !defined(COMPUTE_SHAPE_FUN) && !defined(THREE_KERNEL_UPDATE)
      //Standard Monolithic Kernel
      SolidMechanicsLagrangianFEMKernels::ArrayOfObjectsKernel<kernelPol>(NoElem,elementList, dt, elemsToNodes,
                                                                          u, uhat, dNdX, constitutiveMap, devStressData, meanStress,shearModulus,
                                                                          bulkModulus, detJ, acc, myUpdate, nx, nx, nx);
#elif defined(COMPUTE_SHAPE_FUN) && !defined(THREE_KERNEL_UPDATE)
      SolidMechanicsLagrangianFEMKernels::ArrayOfObjectsKernel_Shape<kernelPol>(NoElem,elementList, dt, elemsToNodes,
                                                                                u, uhat, VX, P, constitutiveMap, devStressData, meanStress,shearModulus,
                                                                                bulkModulus, detJ, acc, myUpdate, nx, nx, nx);            
#elif defined(THREE_KERNEL_UPDATE)

      //Kinematic step
      SolidMechanicsLagrangianFEMKernels::ArrayOfObjects_KinematicKernel<kernelPol>(NoElem,elementList, dt, elemsToNodes,
                                                                                    u, uhat, dNdX, constitutiveMap, devStressData, meanStress,shearModulus,
                                                                                    bulkModulus, detJ, acc, Dadt, Rot, detF, inverseF,
                                                                                    nx, nx, nx);
      //Constitutive update
      SolidMechanicsLagrangianFEMKernels::ConstitutiveUpdateKernel<kernelPol>(NoElem, elementList, Dadt, Rot, constitutiveMap, devStressData,
                                                                              meanStress, shearModulus, bulkModulus);      
      //Integration step
      SolidMechanicsLagrangianFEMKernels::ArrayOfObjects_IntegrationKernel<kernelPol>(NoElem,elementList, dt, elemsToNodes,
                                                                                      u, uhat, dNdX, constitutiveMap, devStressData, meanStress,shearModulus,
                                                                                      bulkModulus, detJ, acc, Dadt, Rot, detF, inverseF,
                                                                                      nx, nx, nx);
#endif

#if defined (USE_CUDA)
      cudaDeviceSynchronize();
#endif      
      end = omp_get_wtime();
      diff = end-start; 
      if( diff < myMin) myMin = diff;
    }

#if defined(USE_CUDA)
    std::cout<<"Computing on the GPU"<<std::endl;
#else
    std::cout<<"Computing on the CPU"<<std::endl;
#endif
  
  std::cout<<"Run time = "<<myMin<<" sec"<<std::endl;
  std::cout<<"No of Elements:= "<<NoElem<<std::endl;
#endif//DELETE

  //
  //Free data
  //
#if !defined(USE_GEOSX_ARRAY)
  memoryManager::deallocate(elementList);
  memoryManager::deallocate(constitutiveMap);
  memoryManager::deallocate(elemsToNodes);
  memoryManager::deallocate(VX);

#if !defined(COMPUTE_SHAPE_FUN)  
  memoryManager::deallocate(dNdX);
#endif

  memoryManager::deallocate(u);
  memoryManager::deallocate(uhat);
  memoryManager::deallocate(acc);
  
  memoryManager::deallocate(detJ);
  memoryManager::deallocate(meanStress);
  memoryManager::deallocate(devStressData);
  
#if defined(THREE_KERNEL_UPDATE)
  memoryManager::deallocate(Dadt);
  memoryManager::deallocate(Rot);
  memoryManager::deallocate(detF);
  memoryManager::deallocate(inverseF);
#endif  
#endif


#endif //DELETE
  return 0;
}
