#include "RAJA/RAJA.hpp"
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

#include "../../src/SolidMechanicsLagrangianFEMKernels_brief_impl.hpp"

//
//Driver for the GEOSX proxy app
//Assumes an Objects of Array data layout
//

int main(int argc, char* const argv[])
{

#if defined(USE_GEOSX_ARRAY)
  std::cout<<"Using GEOSX ARRAY"<<std::endl;
#elif defined(USE_RAJA_VIEW)
  std::cout<<"Using RAJA VIEW"<<std::endl;
#else
  std::cout<<"Using RAW PTRS"<<std::endl;
#endif

  if(argc != 3)
    {
      std::cout<<"usage ./main NoElem Iter"<<std::endl;
      exit(-1);
    }

  size_t dataAllocated = 0.0;
  localIndex Kx      = atoi(argv[1]); 
  localIndex Niter   =  atoi(argv[2]);

  std::cout<<"GEOSX mini-app: Object of Arrays data structures"<<std::endl;

  real64 dt = 0.125; //step size
  const localIndex nx = Kx+1; //Number of nodes in a cartesian dimension
  const localIndex NoElem = Kx*Kx*Kx; //Total number of elements
  const localIndex numNodes = nx*nx*nx; //Total number of nodes

  //
  //Generate an element list
  //
  localIndex * const elementList = memoryManager::allocate<localIndex>(NoElem, dataAllocated);
  for(localIndex i=0; i<NoElem; ++i) elementList[i] = i;
  

  //
  //Allocate space for constitutive map
  //
  //Allocate space for constitutive map
  //
#if defined(USE_GEOSX_ARRAY)
  LvArray::Array<localIndex, 2, localIndex> _constitutiveMap(NoElem, NUMQUADPTS);
  LvArray::ArrayView<localIndex,2,localIndex> & constitutiveMap = _constitutiveMap;

#elif defined(USE_RAJA_VIEW)

  localIndex *  _constitutiveMap = memoryManager::allocate<localIndex>(NUMQUADPTS*NoElem, dataAllocated); 
  RAJA::View<localIndex, RAJA::Layout<2,localIndex,1>> constitutiveMap(_constitutiveMap, NoElem, NUMQUADPTS);
#else

  localIndex *  constitutiveMap = memoryManager::allocate<localIndex>(NUMQUADPTS*NoElem, dataAllocated); 
#endif

  //
  //Generate space for an element to node list
  //
#if defined(USE_GEOSX_ARRAY)
  LvArray::Array<localIndex, 2, localIndex> _elemsToNodes(NoElem, NODESPERELEM);
  LvArray::ArrayView<localIndex,2,localIndex> & elemsToNodes = _elemsToNodes;

#elif defined(USE_RAJA_VIEW)
  localIndex  *  elemsToNodes = memoryManager::allocate<localIndex>(NODESPERELEM*NoElem, dataAllocated);
#else
  localIndex  *  elemsToNodes = memoryManager::allocate<localIndex>(NODESPERELEM*NoElem, dataAllocated);
#endif

  //
  //Allocate space for a list of vertices, generate a mesh, and populate the constitutive map
  //
#if defined(USE_GEOSX_ARRAY)
  LvArray::Array<real64, 2, localIndex> _VX(numNodes, LOCAL_DIM);
  LvArray::ArrayView<real64,2, localIndex> & VX = _VX;
#elif defined(USE_RAJA_VIEW)

  geosxData _VX = memoryManager::allocate<real64>(numNodes*LOCAL_DIM, dataAllocated);
  RAJA::View<real64, RAJA::Layout<2, localIndex, 1> > VX(_VX, numNodes, LOCAL_DIM);
#else
  geosxData VX = memoryManager::allocate<real64>(numNodes*LOCAL_DIM, dataAllocated);
#endif
  meshGen(VX, elemsToNodes,constitutiveMap,Kx);

  //
  //Precompute evaluate the quadrature points at every basis function
  //with respect to the parent basis function
  //
  P_Wrapper P; 
  generateP(P, NODESPERELEM, NUMQUADPTS);
  

  //
  //Allocate space for shape function derivatives and compute their values
  //
#if defined(USE_GEOSX_ARRAY)
  LvArray::Array<real64, 3, localIndex> _dNdX_x(NoElem, NUMQUADPTS, NODESPERELEM);
  LvArray::Array<real64, 3, localIndex> _dNdX_y(NoElem, NUMQUADPTS, NODESPERELEM);
  LvArray::Array<real64, 3, localIndex> _dNdX_z(NoElem, NUMQUADPTS, NODESPERELEM);
  LvArray::ArrayView<real64, 3, localIndex> & dNdX_x = _dNdX_x;
  LvArray::ArrayView<real64, 3, localIndex> & dNdX_y = _dNdX_y;
  LvArray::ArrayView<real64, 3, localIndex> & dNdX_z = _dNdX_z;

#elif defined(USE_RAJA_VIEW)
  geosxData _dNdX_x = memoryManager::allocate<real64>(NODESPERELEM*NUMQUADPTS*NoElem*LOCAL_DIM, dataAllocated);
  geosxData _dNdX_y = memoryManager::allocate<real64>(NODESPERELEM*NUMQUADPTS*NoElem*LOCAL_DIM, dataAllocated);
  geosxData _dNdX_z = memoryManager::allocate<real64>(NODESPERELEM*NUMQUADPTS*NoElem*LOCAL_DIM, dataAllocated);
  RAJA::View<real64, RAJA::Layout<3, localIndex, 2> > dNdX_x(_dNdX_x, NoElem, NUMQUADPTS, NODESPERELEM);
  RAJA::View<real64, RAJA::Layout<3, localIndex, 2> > dNdX_y(_dNdX_y, NoElem, NUMQUADPTS, NODESPERELEM);
  RAJA::View<real64, RAJA::Layout<3, localIndex, 2> > dNdX_z(_dNdX_z, NoElem, NUMQUADPTS, NODESPERELEM);
#else
  geosxData dNdX_x = memoryManager::allocate<real64>(NODESPERELEM*NUMQUADPTS*NoElem, dataAllocated);
  geosxData dNdX_y = memoryManager::allocate<real64>(NODESPERELEM*NUMQUADPTS*NoElem, dataAllocated);
  geosxData dNdX_z = memoryManager::allocate<real64>(NODESPERELEM*NUMQUADPTS*NoElem, dataAllocated);
#endif

  make_dNdX(dNdX_x, dNdX_y, dNdX_z, VX, elemsToNodes, NoElem, NUMQUADPTS, NODESPERELEM);

  //
  //Allocate space for nodal degrees of freedom as Object of Arrays
  //
#if defined(USE_GEOSX_ARRAY)
  LvArray::Array<real64, 1, localIndex> _u_x(numNodes);
  LvArray::Array<real64, 1, localIndex> _u_y(numNodes);
  LvArray::Array<real64, 1, localIndex> _u_z(numNodes);

  LvArray::Array<real64, 1, localIndex> _uhat_x(numNodes);
  LvArray::Array<real64, 1, localIndex> _uhat_y(numNodes);
  LvArray::Array<real64, 1, localIndex> _uhat_z(numNodes);

  LvArray::Array<real64, 1, localIndex> _acc_x(numNodes);
  LvArray::Array<real64, 1, localIndex> _acc_y(numNodes);
  LvArray::Array<real64, 1, localIndex> _acc_z(numNodes); 

  LvArray::ArrayView<real64, 1, localIndex> & u_x = _u_x;
  LvArray::ArrayView<real64, 1, localIndex> & u_y = _u_y;
  LvArray::ArrayView<real64, 1, localIndex> & u_z = _u_z;

  LvArray::ArrayView<real64, 1, localIndex> & uhat_x = _uhat_x;
  LvArray::ArrayView<real64, 1, localIndex> & uhat_y = _uhat_y;
  LvArray::ArrayView<real64, 1, localIndex> & uhat_z = _uhat_z;

  LvArray::ArrayView<real64, 1, localIndex> & acc_x = _acc_x;
  LvArray::ArrayView<real64, 1, localIndex> & acc_y = _acc_y;
  LvArray::ArrayView<real64, 1, localIndex> & acc_z = _acc_z;
  std::memset(u_x.data(),0, numNodes*sizeof(real64));
  std::memset(u_y.data(),0, numNodes*sizeof(real64));
  std::memset(u_z.data(),0, numNodes*sizeof(real64));
  std::memset(uhat_x.data(),0, numNodes*sizeof(real64));
  std::memset(uhat_y.data(),0, numNodes*sizeof(real64));
  std::memset(uhat_z.data(),0, numNodes*sizeof(real64));
  std::memset(acc_x.data(),0, numNodes*sizeof(real64));
  std::memset(acc_y.data(),0, numNodes*sizeof(real64));
  std::memset(acc_z.data(),0, numNodes*sizeof(real64));

#elif defined(USE_RAJA_VIEW)

  geosxData _u_x = memoryManager::allocate<real64>(numNodes, dataAllocated);
  geosxData _u_y = memoryManager::allocate<real64>(numNodes, dataAllocated);
  geosxData _u_z = memoryManager::allocate<real64>(numNodes, dataAllocated);

  geosxData _uhat_x = memoryManager::allocate<real64>(numNodes, dataAllocated);
  geosxData _uhat_y = memoryManager::allocate<real64>(numNodes, dataAllocated);
  geosxData _uhat_z = memoryManager::allocate<real64>(numNodes, dataAllocated);

  geosxData _acc_x = memoryManager::allocate<real64>(numNodes, dataAllocated);
  geosxData _acc_y = memoryManager::allocate<real64>(numNodes, dataAllocated);
  geosxData _acc_z = memoryManager::allocate<real64>(numNodes, dataAllocated);
  
  RAJA::View<real64, RAJA::Layout<1, localIndex, 0>> u_x(_u_x, numNodes);
  RAJA::View<real64, RAJA::Layout<1, localIndex, 0>> u_y(_u_y, numNodes);
  RAJA::View<real64, RAJA::Layout<1, localIndex, 0>> u_z(_u_z, numNodes);

  RAJA::View<real64, RAJA::Layout<1, localIndex, 0>> uhat_x(_uhat_x, numNodes);
  RAJA::View<real64, RAJA::Layout<1, localIndex, 0>> uhat_y(_uhat_y, numNodes);
  RAJA::View<real64, RAJA::Layout<1, localIndex, 0>> uhat_z(_uhat_z, numNodes);

  RAJA::View<real64, RAJA::Layout<1, localIndex, 0>> acc_x(_acc_x, numNodes);
  RAJA::View<real64, RAJA::Layout<1, localIndex, 0>> acc_y(_acc_y, numNodes);
  RAJA::View<real64, RAJA::Layout<1, localIndex, 0>> acc_z(_acc_z, numNodes);

  std::memset(_u_x,0, numNodes*sizeof(real64));
  std::memset(_u_y,0, numNodes*sizeof(real64));
  std::memset(_u_z,0, numNodes*sizeof(real64));

  std::memset(_uhat_x,0, numNodes*sizeof(real64));
  std::memset(_uhat_y,0, numNodes*sizeof(real64));
  std::memset(_uhat_z,0, numNodes*sizeof(real64));

  std::memset(_acc_x,0, numNodes*sizeof(real64));
  std::memset(_acc_y,0, numNodes*sizeof(real64));
  std::memset(_acc_z,0, numNodes*sizeof(real64));
 
#else
  geosxData u_x = memoryManager::allocate<real64>(numNodes, dataAllocated);
  geosxData u_y = memoryManager::allocate<real64>(numNodes, dataAllocated);
  geosxData u_z = memoryManager::allocate<real64>(numNodes, dataAllocated);
  
  geosxData uhat_x = memoryManager::allocate<real64>(numNodes, dataAllocated);
  geosxData uhat_y = memoryManager::allocate<real64>(numNodes, dataAllocated);
  geosxData uhat_z = memoryManager::allocate<real64>(numNodes, dataAllocated);
  
  geosxData acc_x = memoryManager::allocate<real64>(numNodes, dataAllocated);
  geosxData acc_y = memoryManager::allocate<real64>(numNodes, dataAllocated);
  geosxData acc_z = memoryManager::allocate<real64>(numNodes, dataAllocated);
  
  std::memset(u_x,0,numNodes*sizeof(real64));
  std::memset(u_y,0,numNodes*sizeof(real64));
  std::memset(u_z,0,numNodes*sizeof(real64));
  
  std::memset(uhat_x,0,numNodes*sizeof(real64));
  std::memset(uhat_y,0,numNodes*sizeof(real64));
  std::memset(uhat_z,0,numNodes*sizeof(real64));
  
  std::memset(acc_x,0,numNodes*sizeof(real64));
  std::memset(acc_y,0,numNodes*sizeof(real64));
  std::memset(acc_z,0,numNodes*sizeof(real64));
#endif


  //
  //Allocate and set physical parameters
  //
  localIndex noSymEnt = 6;
  real64 bulkModulus = 10, shearModulus=20;

#if defined(USE_GEOSX_ARRAY)
  
  LvArray::Array<real64,2,localIndex> _detJ(NoElem, NUMQUADPTS);
  LvArray::Array<real64,3,localIndex> _devStressData(NoElem, NUMQUADPTS, noSymEnt);
  LvArray::Array<real64,1,localIndex> _meanStress(NoElem*NUMQUADPTS); //Reformulate to 1D

  LvArray::ArrayView<real64,2,localIndex> & detJ          = _detJ;
  LvArray::ArrayView<real64,3,localIndex> & devStressData = _devStressData;  
  LvArray::ArrayView<real64,1,localIndex> & meanStress    = _meanStress;

  std::memset(detJ.data(), 1 ,NUMQUADPTS * NoElem * sizeof(real64));
  std::memset(meanStress.data(), 1 ,NUMQUADPTS * NoElem * sizeof(real64));
  std::memset(devStressData.data(), 1 , noSymEnt * NUMQUADPTS * NoElem * sizeof(real64));

#elif defined(USE_RAJA_VIEW)
  
  geosxData _detJ            = memoryManager::allocate<real64>(NUMQUADPTS*NoElem, dataAllocated);
  geosxData _meanStress      = memoryManager::allocate<real64>(NUMQUADPTS*NoElem, dataAllocated);
  geosxData _devStressData   = memoryManager::allocate<real64>(noSymEnt*NUMQUADPTS*NoElem, dataAllocated);

  RAJA::View<real64, RAJA::Layout<2,localIndex,1>> detJ(_detJ, NoElem, NUMQUADPTS);
  RAJA::View<real64, RAJA::Layout<3,localIndex,2>> devStressData(_devStressData, NoElem, NUMQUADPTS, noSymEnt);
  RAJA::View<real64, RAJA::Layout<1,localIndex,0>> meanStress(_meanStress, NoElem*NUMQUADPTS);

  std::memset(_detJ, 1 ,NUMQUADPTS * NoElem * sizeof(real64));
  std::memset(_meanStress, 1 ,NUMQUADPTS * NoElem * sizeof(real64));
  std::memset(_devStressData, 1 , noSymEnt * NUMQUADPTS * NoElem * sizeof(real64));  

#else
  geosxData detJ            = memoryManager::allocate<real64>(NUMQUADPTS*NoElem, dataAllocated);
  geosxData meanStress      = memoryManager::allocate<real64>(NUMQUADPTS*NoElem, dataAllocated);
  geosxData devStressData   = memoryManager::allocate<real64>(noSymEnt*NUMQUADPTS*NoElem, dataAllocated);

  std::memset(detJ, 1 ,NUMQUADPTS * NoElem * sizeof(real64));
  std::memset(meanStress, 1 ,NUMQUADPTS * NoElem * sizeof(real64));
  std::memset(devStressData, 1 , noSymEnt * NUMQUADPTS * NoElem * sizeof(real64));
#endif  

  //
  //Set up function pointer to constitutive relationship
  //
#if defined(USE_GPU)
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
  for(localIndex it=0; it<Niter; ++it)
    {
      start = omp_get_wtime();


      //Monolithic Kernel: Assumes precomputed shape function derivatives
    SolidMechanicsLagrangianFEMKernels::ObjectOfArraysKernel<kernelPol>(NoElem, elementList, dt, elemsToNodes,
                                        u_x, u_y, u_z, uhat_x, uhat_y, uhat_z,
                                        dNdX_x, dNdX_y, dNdX_z, constitutiveMap, devStressData,
                                        meanStress,shearModulus, bulkModulus, detJ, acc_x, acc_y, acc_z, myUpdate, nx, nx, nx);

#if defined(USE_GPU)
      cudaDeviceSynchronize();
#endif      
      end = omp_get_wtime();
      diff = end-start; 
      if( diff < myMin) myMin = diff;
    }

  
#if defined(USE_GPU)
    std::cout<<"Computing on the GPU"<<std::endl;
#else
    std::cout<<"Computing on the CPU"<<std::endl;
#endif 

  std::cout<<"Run time = "<<myMin<<" sec"<<std::endl;
  std::cout<<"No of Elements:= "<<NoElem<<std::endl;

  //  
    //Free the data
    //
  /*
    memoryManager::deallocate(elementList);
    memoryManager::deallocate(constitutiveMap);
    memoryManager::deallocate(elemsToNodes);
    memoryManager::deallocate(VX);

#if !defined(COMPUTE_SHAPE_FUN)    
    memoryManager::deallocate(dNdX_x);
    memoryManager::deallocate(dNdX_y);
    memoryManager::deallocate(dNdX_z);
#endif    

    memoryManager::deallocate(u_x);
    memoryManager::deallocate(u_y);
    memoryManager::deallocate(u_z);
    memoryManager::deallocate(uhat_x);
    memoryManager::deallocate(uhat_y);
    memoryManager::deallocate(uhat_z);
    memoryManager::deallocate(acc_x);
    memoryManager::deallocate(acc_y);
    memoryManager::deallocate(acc_z);
    
    memoryManager::deallocate(detJ);
    memoryManager::deallocate(meanStress);
    memoryManager::deallocate(devStressData);
        
#if defined(THREE_KERNEL_UPDATE)
    memoryManager::deallocate(Dadt);
    memoryManager::deallocate(Rot);
    memoryManager::deallocate(detF);
    memoryManager::deallocate(inverseF);
#endif  
  */
   
  return 0;
}
