/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

#ifndef __SOLID_MECHANICS_LAGRANGIAN_BRIEF_FEM_KERNELS_HPP__
#define __SOLID_MECHANICS_LAGRANGIAN_BRIEF_FEM_KERNELS_HPP__
#include "../miniApps/SolidMechanicsLagrangianFEM-MiniApp/MatrixMath_impl.hpp"
#include "../miniApps/SolidMechanicsLagrangianFEM-MiniApp/Layout.hpp"
#include "../miniApps/SolidMechanicsLagrangianFEM-MiniApp/ShapeFun_impl.hpp"
#include "../miniApps/SolidMechanicsLagrangianFEM-MiniApp/ConstitutiveUpdate_impl.hpp"
#include "../../rajaInterface/GEOS_RAJA_Interface.hpp"

/*
  Collection of Solid Mechanics Kernels

  Here we explore storing nodal degrees of freedom and shape
  function derivatives in two different formats. 
  
  Array of Objects vs Object of Arrays 

  Array of objects stores nodal degrees of freedom in a single 
  array via and x,y,z format:

        u = [x_0, y_0 , z_0 , x_1 , y_1 , z_1, .... ]

  Object of Arrays stores nodal degrees of freedom by cartsian 
  dimension:

        u_x = [x_0, x_1, x_2, . . . . ]
        u_y = [y_0, y_1, y_2, . . . . ]
        u_z = [z_0, z_1, z_2, . . . . ]

   The two main kernels in this header are the
   ObjectsOfArray and ArrayOfObjects kernel. 
   They store nodal and shape function derivatives 
   in the follower manner:

                           Nodal Dofs  | Quad Dofs
                          _________________________
    ObjectOfArraysKernel  |   ObjofArr | ObjofArr | 
                          -------------|-----------
    ArraysOfObjectsKernel |   ArrofObj | ArrofObj |
                          -------------------------

   We also include breaking up the kernel into three steps.
   
   1. Kinematic step
   2. Constitutive update step
   3. Integration step. 

   The consequence of breaking up the monolithic kernel
   into three kernels is the extra storage needed to store
   intermediate computations. 

   Lastly, An underlying assumption of these kernels 
   is that data access is done throught the paranthesis operator
   for multi-dimensional arrays. Changing data layouts
   for the fastest running index is accomplished by
   defining/undefing macros in Layout.hpp. The element index
   may either be the fast or slowest runnin index.
                          
*/

using namespace LvArray;
namespace SolidMechanicsLagrangianFEMKernels{

//using namespace SolidMechanicsLagrangianFEMKernels;

///
///Solid mechanics update kernel with nodal degrees of freedom and
///shape function derivatives stored in object of arrays format. 
///Computations are carried out in a monothilic kernel.
///

template<typename Pol>
RAJA_INLINE
void ObjectOfArraysKernel(localIndex noElem, geosxIndex elemList, real64 dt,
                          const localIndex * elemsToNodes, geosxData iu_x,
                          geosxData iu_y, geosxData iu_z,
                          geosxData iuhat_x,
                          geosxData iuhat_y, geosxData iuhat_z,
                          geosxData idNdX_x,
                          geosxData idNdX_y,geosxData idNdX_z,
                          localIndex const * iconstitutiveMap, geosxData idevStressData,
                          geosxData imeanStress, real64 ishearModulus, real64 ibulkModulus,
                          const real64 * idetJ, geosxData iacc_x,
                          geosxData iacc_y, geosxData iacc_z, constUpdate updateState_ptr, 
                          localIndex nx=2, localIndex ny=2, localIndex nz=2)
{

  /*
  geosx::forall_in_set<Pol>(elemList, noElem, GEOSX_LAMBDA (globalIndex k) {
      
      real64 uhat_local_x[NODESPERELEM];
      real64 uhat_local_y[NODESPERELEM];
      real64 uhat_local_z[NODESPERELEM];

      real64 u_local_x[NODESPERELEM];
      real64 u_local_y[NODESPERELEM];
      real64 u_local_z[NODESPERELEM];
      
      real64 f_local_x[NODESPERELEM] = {0};
      real64 f_local_y[NODESPERELEM] = {0};
      real64 f_local_z[NODESPERELEM] = {0};


#if defined(STRUCTURED_GRID)
       localIndex nodeList[NODESPERELEM];       
       structuredElemToNodes(nodeList,k,nx,ny,nz);
#else
       const localIndex *nodeList = (&elemsToNodes[NODESPERELEM*k]);
#endif       
      
       //Copy Global To Local
       GlobalToLocal(nodeList, k, 
                     u_local_x, u_local_y, u_local_z,
                     uhat_local_x, uhat_local_y, uhat_local_z,
                     iu_x, iu_y, iu_z, iuhat_x, iuhat_y, iuhat_z);
              
      //Compute Quadrature
      for(localIndex q=0; q<NUMQUADPTS; ++q)
        {


          real64 dUhatdX[LOCAL_DIM][LOCAL_DIM] = {{0.0}};
          real64 dUdX[LOCAL_DIM][LOCAL_DIM] = {{0.0}};

          //Calculate gradient
          CalculateGradient(dUdX, u_local_x, u_local_y, u_local_z,
                            idNdX_x, idNdX_y, idNdX_z, k, q, noElem);
          
          CalculateGradient(dUhatdX, uhat_local_x, uhat_local_y, uhat_local_z,
                            idNdX_x, idNdX_y, idNdX_z, k, q, noElem);                            

          real64 F[LOCAL_DIM][LOCAL_DIM];
          real64 Finv[LOCAL_DIM][LOCAL_DIM];
          real64 L[LOCAL_DIM][LOCAL_DIM];

          {            
            real64 dvdX[LOCAL_DIM][LOCAL_DIM];

            
            for(localIndex ty=0; ty<LOCAL_DIM; ++ty){
              for(localIndex tx=0; tx<LOCAL_DIM; ++tx){
                dvdX[ty][tx] = dUhatdX[ty][tx]*(1.0/dt);
              }
            }
            
            //calculate du/dX
            for(localIndex row=0; row<LOCAL_DIM; ++row){
              for(localIndex col=0; col<LOCAL_DIM; ++col){
                F[row][col] = dUhatdX[row][col];
              }
            }
            
            for(localIndex row=0; row<LOCAL_DIM; ++row){
              for(localIndex col=0; col<LOCAL_DIM; ++col){
                      F[row][col] *= 0.5;
              }
            }
            
            for(localIndex row=0; row<LOCAL_DIM; ++row){
              for(localIndex col=0; col<LOCAL_DIM; ++col){
                F[row][col] += dUdX[row][col];
              }
            }
            
            for(localIndex tx=0; tx<LOCAL_DIM; ++tx){
              F[tx][tx] += 1.0;
            }

            
            Finverse(F, Finv);
            

            AijBjk(dvdX,Finv,L);
          }

          //Calculate gradient (end of step)
          for(localIndex ty=0; ty<LOCAL_DIM; ++ty){
            for(localIndex tx=0; tx<LOCAL_DIM; ++tx){
              F[ty][tx] = dUhatdX[ty][tx] + dUdX[ty][tx];
            }
          }

          for(int tx=0; tx<LOCAL_DIM; ++tx){ F[tx][tx] += 1.0;};
          
          real64 detF = det<real64>(F);
          Finverse<real64> (F, Finv);

          real64 Rot[LOCAL_DIM][LOCAL_DIM];
          real64 Dadt[LOCAL_DIM][LOCAL_DIM];

          //-------------[Hughues Winget]--------------
          HughesWinget(Rot, Dadt, L, dt);
          //-------------------------------------------
          
          //-------------[Constitutive update]-------------
          localIndex m = iconstitutiveMap(k,q);          
          //UpdateStatePoint(Dadt,Rot,m, q,k, idevStressData, imeanStress,ishearModulus, ibulkModulus, noElem);
          updateState_ptr(Dadt,Rot,m, q,k, idevStressData, imeanStress, ishearModulus, ibulkModulus, noElem);
          //-------------------------------------------          

          real64 TotalStress[LOCAL_DIM][LOCAL_DIM];

          TotalStress[0][0] = idevStressData(k,q,0);
          TotalStress[1][0] = idevStressData(k,q,1);
          TotalStress[1][1] = idevStressData(k,q,2);
          TotalStress[2][0] = idevStressData(k,q,3);
          TotalStress[2][1] = idevStressData(k,q,4);
          TotalStress[2][2] = idevStressData(k,q,5);

          TotalStress[0][1] = TotalStress[1][0];
          TotalStress[0][2] = TotalStress[2][0];
          TotalStress[1][2] = TotalStress[2][1];

          for(localIndex i=0; i<LOCAL_DIM; ++i)
            {
              TotalStress[i][i] += imeanStress[m];
            }

          Integrate(f_local_x, f_local_y, f_local_z,
                    idetJ(k,q), detF, Finv, TotalStress, idNdX_x, idNdX_y, idNdX_z, k, q, noElem);
        }//end of quadrature

      //Atomic policy
      AddLocalToGlobal<atomicPol>(nodeList,f_local_x, f_local_y, f_local_z, iacc_x, iacc_y, iacc_z);

    });
  */

}
      
///
///Solid mechanics update kernel with nodal degrees of freedom and
///shape function derivatives stored as an object of structs format.
///All computations are done in a monolithic kernel.
///      

#if defined(USE_GEOSX_ARRAY)
template<typename Pol>
RAJA_INLINE void ArrayOfObjectsKernel(localIndex noElem,  ArrayView<localIndex, 1, localIndex> elemList,
                                      real64 dt, ArrayView<localIndex,2,localIndex> elemsToNodes, 
                                      ArrayView<real64, 2, localIndex> iu,
                                      ArrayView<real64, 2, localIndex> iuhat,
                                      ArrayView<real64, 4, localIndex> idNdX,
                                      ArrayView<localIndex,2,localIndex> iconstitutiveMap,
                                      ArrayView<real64,3,localIndex> idevStressData,
                                      ArrayView<real64,1,localIndex> imeanStress, real64 ishearModulus, real64 ibulkModulus,
                                      ArrayView<real64,2,localIndex> idetJ,
                                      ArrayView<real64, 2, localIndex> iacc, constUpdate updateState_ptr,
                                      localIndex nx=2, localIndex ny=2, localIndex nz=2)

#elif defined(USE_RAJA_VIEW)
template<typename Pol>
RAJA_INLINE void ArrayOfObjectsKernel(localIndex noElem, geosxIndex elemList,
                                      real64 dt, const localIndex * const elemsToNodes,
                                      RAJA::View<real64, RAJA::Layout<2, localIndex,1> > iu,
                                      RAJA::View<real64, RAJA::Layout<2, localIndex,1> > iuhat,
                                      RAJA::View<real64, RAJA::Layout<4, localIndex,3> > idNdX,
                                      RAJA::View<localIndex,RAJA::Layout<2,localIndex,1> > iconstitutiveMap,
                                      RAJA::View<real64,RAJA::Layout<3,localIndex,2> > idevStressData,
                                      RAJA::View<real64,RAJA::Layout<1,localIndex,0> > imeanStress, real64 ishearModulus, real64 ibulkModulus,
                                      RAJA::View<real64,RAJA::Layout<2,localIndex,1> > idetJ,
                                      RAJA::View<real64, RAJA::Layout<2, localIndex,1 > > iacc, constUpdate updateState_ptr,
                                      localIndex nx=2, localIndex ny=2, localIndex nz=2)

#else //Using Raw pointers
template<typename Pol>
RAJA_INLINE void ArrayOfObjectsKernel(localIndex noElem, geosxIndex elemList, real64 dt,
                                      const localIndex * const elemsToNodes, geosxData iu,
                                      geosxData iuhat, geosxData idNdX,
                                      localIndex const * iconstitutiveMap, geosxData idevStressData,
                                      geosxData imeanStress, real64 ishearModulus, real64 ibulkModulus,
                                      const real64 * idetJ, geosxData iacc, constUpdate updateState_ptr,
                                      localIndex nx=2, localIndex ny=2, localIndex nz=2)
#endif
{

  //geosx::forall_in_set<Pol>(elemList.data(), noElem, GEOSX_LAMBDA (globalIndex k) {
  geosx::forall_in_range<Pol>(0, noElem, GEOSX_LAMBDA (globalIndex k) {

      real64 uhat_local[LOCAL_DIM*NODESPERELEM];
      real64 u_local[LOCAL_DIM*NODESPERELEM];
      real64 f_local[LOCAL_DIM*NODESPERELEM] = {0};

#if defined(STRUCTURED_GRID)
      localIndex nodeList[NODESPERELEM];       
      structuredElemToNodes(nodeList,k,nx,ny,nz);
#else

#if defined(USE_GEOSX_ARRAY)
      const localIndex *nodeList = (&elemsToNodes.data()[NODESPERELEM*k]);
#else      
      const localIndex *nodeList = (&elemsToNodes[NODESPERELEM*k]);
#endif//GEOSX_ARRAY_CHECK

#endif//STRUCTURED_GRID_CHECK


       //Copy Global to Local
       GlobalToLocal(nodeList, k, 
                     u_local,  uhat_local, iu, iuhat);

      //Compute Quadrature
      for(localIndex q=0; q<NUMQUADPTS; ++q)
        {

          real64 dUdX[LOCAL_DIM][LOCAL_DIM] = {{0.0}};          
          real64 dUhatdX[LOCAL_DIM][LOCAL_DIM] = {{0.0}};

          //Calculate Gradient
          CalculateGradient(dUdX, u_local, idNdX, k, q, noElem);
          CalculateGradient(dUhatdX, uhat_local, idNdX, k, q, noElem);


          real64 F[LOCAL_DIM][LOCAL_DIM];
          real64 Finv[LOCAL_DIM][LOCAL_DIM];
          real64 L[LOCAL_DIM][LOCAL_DIM];
          
          //Compute L
          {            
            real64 dvdX[LOCAL_DIM][LOCAL_DIM];
            
            for(localIndex ty=0; ty<LOCAL_DIM; ++ty){
              for(localIndex tx=0; tx<LOCAL_DIM; ++tx){
                dvdX[ty][tx] = dUhatdX[ty][tx]*(1.0/dt);
              }
            }
            
            //calculate du/dX
            for(localIndex row=0; row<LOCAL_DIM; ++row){
              for(localIndex col=0; col<LOCAL_DIM; ++col){
                F[row][col] = dUhatdX[row][col];
              }
            }
            
            for(localIndex row=0; row<LOCAL_DIM; ++row){
              for(localIndex col=0; col<LOCAL_DIM; ++col){
                      F[row][col] *= 0.5;
              }
            }
            
            for(localIndex row=0; row<LOCAL_DIM; ++row){
              for(localIndex col=0; col<LOCAL_DIM; ++col){
                F[row][col] += dUdX[row][col];
              }
            }
            
            for(localIndex tx=0; tx<LOCAL_DIM; ++tx){
              F[tx][tx] += 1.0;
            }
            
            Finverse(F, Finv);
            

            AijBjk(dvdX,Finv,L);
          }

          //Calculate gradient (end of step)
          for(localIndex ty=0; ty<LOCAL_DIM; ++ty){
            for(localIndex tx=0; tx<LOCAL_DIM; ++tx){
              F[ty][tx] = dUhatdX[ty][tx] + dUdX[ty][tx];
            }
          }

          for(int tx=0; tx<LOCAL_DIM; ++tx){ F[tx][tx] += 1.0;};
          
          real64 detF = det<real64>(F);
          Finverse<real64> (F, Finv);


          real64 Rot[LOCAL_DIM][LOCAL_DIM];
          real64 Dadt[LOCAL_DIM][LOCAL_DIM];

          //-------------[Hughues Winget]--------------
          HughesWinget(Rot, Dadt, L, dt);
          //-------------------------------------------
          

          //-------------[Constitutive update]-------------
          localIndex m = iconstitutiveMap(k,q);
          //UpdateStatePoint(Dadt,Rot,m, q,k, idevStressData,
          //imeanStress,ishearModulus, ibulkModulus,noElem);
          
          updateState_ptr(Dadt,Rot, m, q, k, idevStressData,
                          imeanStress, ishearModulus, ibulkModulus, noElem);

          //-------------------------------------------

          real64 TotalStress[LOCAL_DIM][LOCAL_DIM];          
          TotalStress[0][0] = idevStressData(k,q,0);
          TotalStress[1][0] = idevStressData(k,q,1);
          TotalStress[1][1] = idevStressData(k,q,2);
          TotalStress[2][0] = idevStressData(k,q,3);
          TotalStress[2][1] = idevStressData(k,q,4);
          TotalStress[2][2] = idevStressData(k,q,5);

          TotalStress[0][1] = TotalStress[1][0];
          TotalStress[0][2] = TotalStress[2][0];
          TotalStress[1][2] = TotalStress[2][1];

          for(localIndex i=0; i<LOCAL_DIM; ++i)
            {

              TotalStress[i][i] += imeanStress(m);
            }

          Integrate(f_local, idetJ(k,q), detF, Finv, TotalStress, idNdX, k, q, noElem); 

        }//end of quadrature

       AddLocalToGlobal<atomicPol>(nodeList, f_local, iacc);

     });
}


}//namespace

#endif
