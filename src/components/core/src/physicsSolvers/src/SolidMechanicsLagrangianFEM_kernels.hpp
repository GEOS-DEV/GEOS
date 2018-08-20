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

#ifndef __SOLID_MECHANICS_LAGRANGIAN_FEM_KERNELS_HPP__
#define __SOLID_MECHANICS_LAGRANGIAN_FEM_KERNELS_HPP__
#include "RAJA/RAJA.hpp"
#include "SolidMechanicsLagrangianFEM-MiniApp/Matrix_Math.hpp"
#include "SolidMechanicsLagrangianFEM-MiniApp/Layout.hpp"
#include "SolidMechanicsLagrangianFEM-MiniApp/Constitutive_Update.hpp"

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



namespace SolidMechanicsLagrangianFEMKernels{

using namespace SolidMechanicsLagrangianFEMKernels;

///
///Solid mechanics update kernel with nodal degrees of freedom and
///shape function derivatives stored in object of arrays format. 
///Computations are carried out in a monothilic kernel.
///
template<typename Pol>
RAJA_INLINE
void ObjectOfArraysKernel(Index_type NoElem, RAJA::TypedListSegment<Index_type> elemList, real64 dt,
                          const Index_type * elemsToNodes, geosxData iu_x,
                          geosxData iu_y, geosxData iu_z,
                          geosxData iuhat_x,
                          geosxData iuhat_y, geosxData iuhat_z,
                          geosxData idNdX_x,
                          geosxData idNdX_y,geosxData idNdX_z,
                          Index_type const * iconstitutiveMap, geosxData idevStressData,
                          geosxData imeanStress, real64 ishearModulus, real64 ibulkModulus,
                          const real64 * idetJ, geosxData iacc_x,
                          geosxData iacc_y, geosxData iacc_z, constUpdate updateState_ptr)
{


#if defined(CUDA)
  RAJA::forall<Pol>
    (elemList,[=] RAJA_DEVICE (Index_type k){
#else
  RAJA::forall<Pol>
    (elemList,[=](Index_type k){
#endif

      real64 uhat_local_x[inumNodesPerElement];
      real64 uhat_local_y[inumNodesPerElement];
      real64 uhat_local_z[inumNodesPerElement];

      real64 u_local_x[inumNodesPerElement];
      real64 u_local_y[inumNodesPerElement];
      real64 u_local_z[inumNodesPerElement];
      
      real64 f_local_x[inumNodesPerElement]; for(int i=0; i<inumNodesPerElement; ++i) f_local_x[i] = 0;
      real64 f_local_y[inumNodesPerElement]; for(int i=0; i<inumNodesPerElement; ++i) f_local_y[i] = 0;
      real64 f_local_z[inumNodesPerElement]; for(int i=0; i<inumNodesPerElement; ++i) f_local_z[i] = 0;

      //Copy Global to Local
      for(Index_type a=0; a<inumNodesPerElement; ++a)
        {
          Index_type id = elemsToNodes[a + inumNodesPerElement*k];
          u_local_x[a] = iu_x[id];
          u_local_y[a] = iu_y[id];
          u_local_z[a] = iu_z[id];
          
          uhat_local_x[a] = iuhat_x[id];
          uhat_local_y[a] = iuhat_y[id];
          uhat_local_z[a] = iuhat_z[id];
        }

      
      //Compute Quadrature
      for(Index_type q=0; q<inumQuadraturePoints; ++q)
        {


          real64 dUhatdX[local_dim][local_dim];
          real64 dUdX[local_dim][local_dim];

          for(Index_type ty=0; ty<local_dim; ++ty){
            for(Index_type tx=0; tx<local_dim; ++tx){
              dUhatdX[ty][tx] = 0.0;
              dUdX[ty][tx] = 0.0;
            }
          }


          for(Index_type a=0; a<inumNodesPerElement; ++a)
            {
              dUdX[0][0] += u_local_x[a]*idNdX_x(k,q,a);
              dUdX[0][1] += u_local_x[a]*idNdX_y(k,q,a);
              dUdX[0][2] += u_local_x[a]*idNdX_z(k,q,a);
              
              dUdX[1][0] += u_local_y[a]*idNdX_x(k,q,a);
              dUdX[1][1] += u_local_y[a]*idNdX_y(k,q,a);
              dUdX[1][2] += u_local_y[a]*idNdX_z(k,q,a);
              
              dUdX[2][0] += u_local_z[a]*idNdX_x(k,q,a);
              dUdX[2][1] += u_local_z[a]*idNdX_y(k,q,a);
              dUdX[2][2] += u_local_z[a]*idNdX_z(k,q,a);
            }

          for(Index_type a=0; a<inumNodesPerElement; ++a)
            {
              dUhatdX[0][0] += uhat_local_x[a]*idNdX_x(k,q,a);
              dUhatdX[0][1] += uhat_local_x[a]*idNdX_y(k,q,a);
              dUhatdX[0][2] += uhat_local_x[a]*idNdX_z(k,q,a);
              
              dUhatdX[1][0] += uhat_local_y[a]*idNdX_x(k,q,a);
              dUhatdX[1][1] += uhat_local_y[a]*idNdX_y(k,q,a);
              dUhatdX[1][2] += uhat_local_y[a]*idNdX_z(k,q,a);
              
              dUhatdX[2][0] += uhat_local_z[a]*idNdX_x(k,q,a);
              dUhatdX[2][1] += uhat_local_z[a]*idNdX_y(k,q,a);
              dUhatdX[2][2] += uhat_local_z[a]*idNdX_z(k,q,a);
            }
          

          real64 F[local_dim][local_dim];
          real64 Finv[local_dim][local_dim];
          real64 L[local_dim][local_dim];

          {            
            real64 dvdX[local_dim][local_dim];

            
            for(Index_type ty=0; ty<local_dim; ++ty){
              for(Index_type tx=0; tx<local_dim; ++tx){
                dvdX[ty][tx] = dUhatdX[ty][tx]*(1.0/dt);
              }
            }
            
            //calculate du/dX
            for(Index_type row=0; row<local_dim; ++row){
              for(Index_type col=0; col<local_dim; ++col){
                F[row][col] = dUhatdX[row][col];
              }
            }
            
            for(Index_type row=0; row<local_dim; ++row){
              for(Index_type col=0; col<local_dim; ++col){
                      F[row][col] *= 0.5;
              }
            }
            
            for(Index_type row=0; row<local_dim; ++row){
              for(Index_type col=0; col<local_dim; ++col){
                F[row][col] += dUdX[row][col];
              }
            }
            
            for(Index_type tx=0; tx<local_dim; ++tx){
              F[tx][tx] += 1.0;
            }

            
            Finverse(F, Finv);
            

            AijBjk(dvdX,Finv,L);
          }

          //Calculate gradient (end of step)
          for(Index_type ty=0; ty<local_dim; ++ty){
            for(Index_type tx=0; tx<local_dim; ++tx){
              F[ty][tx] = dUhatdX[ty][tx] + dUdX[ty][tx];
            }
          }

          for(int tx=0; tx<local_dim; ++tx){ F[tx][tx] += 1.0;};
          
          real64 detF = det<real64>(F);
          Finverse<real64> (F, Finv);

          real64 Rot[local_dim][local_dim];
          real64 Dadt[local_dim][local_dim];


          //-------------[Hughues Winget]--------------
          HughesWinget(Rot, Dadt, L, dt);
          //-------------------------------------------
          
          //-------------[Constitutive update]-------------
          Index_type m = iconstitutiveMap(k,q);          
          //UpdateStatePoint(Dadt,Rot,m, q,k, idevStressData, imeanStress,ishearModulus, ibulkModulus, NoElem);
          updateState_ptr(Dadt,Rot,m, q,k, idevStressData, imeanStress, ishearModulus, ibulkModulus, NoElem);
          //-------------------------------------------          

          real64 TotalStress[local_dim][local_dim];

          TotalStress[0][0] = idevStressData(k,q,0);
          TotalStress[1][0] = idevStressData(k,q,1);
          TotalStress[1][1] = idevStressData(k,q,2);
          TotalStress[2][0] = idevStressData(k,q,3);
          TotalStress[2][1] = idevStressData(k,q,4);
          TotalStress[2][2] = idevStressData(k,q,5);

          TotalStress[0][1] = TotalStress[1][0];
          TotalStress[0][2] = TotalStress[2][0];
          TotalStress[1][2] = TotalStress[2][1];

          for(Index_type i=0; i<local_dim; ++i)
            {
              TotalStress[i][i] += imeanStress[m];
            }

          //---------[Integrate - Function]---------------------
          real64 const integrationFactor = idetJ(k,q)*detF;
          real64 P[local_dim][local_dim];

          AijBkj(TotalStress,Finv,P);
          for(Index_type ty=0; ty<local_dim; ++ty){
            for(Index_type tx=0; tx<local_dim; ++tx){
              P[ty][tx] *= integrationFactor;
            }
          } 


          for(int a=0; a<inumNodesPerElement; ++a){
            //long int id = a + inumNodesPerElement*(q + inumQuadraturePoints*k);
            f_local_x[a] -= P[0][0]*idNdX_x(k,q,a) + P[0][1]*idNdX_y(k,q,a) + P[0][2]*idNdX_z(k,q,a);
            f_local_y[a] -= P[1][0]*idNdX_x(k,q,a) + P[1][1]*idNdX_y(k,q,a) + P[1][2]*idNdX_z(k,q,a);
            f_local_z[a] -= P[2][0]*idNdX_x(k,q,a) + P[2][1]*idNdX_y(k,q,a) + P[2][2]*idNdX_z(k,q,a);
          }

          
        }//end of quadrature


      //-------------[Add Local To Global]----------------
      for(Index_type a=0; a<inumNodesPerElement; ++a)
        {
          
          Index_type id = elemsToNodes[a + inumNodesPerElement*k];
          RAJA::atomic::atomicAdd<atomicPol>(&iacc_x[id],f_local_x[a]);
          RAJA::atomic::atomicAdd<atomicPol>(&iacc_y[id],f_local_y[a]);
          RAJA::atomic::atomicAdd<atomicPol>(&iacc_z[id],f_local_z[a]);
        }          
      //---------------------------------------------------

    });

} 

///
///Solid mechanics update kernel with nodal degrees of freedom and
///shape function derivatives stored in object of arrays format. 
///Computations are carried out in a monothilic kernel. Kernel
///assumes a structured mesh.
///      
template<typename Pol>
RAJA_INLINE
void ObjectOfArraysKernel_Structured(Index_type NoElem, RAJA::TypedListSegment<Index_type> elemList, real64 dt,
                                     Index_type nx, Index_type ny, Index_type nz,
                                     geosxData iu_x,
                                     geosxData iu_y, geosxData iu_z,
                                     geosxData iuhat_x,
                                     geosxData iuhat_y, geosxData iuhat_z,
                                     geosxData idNdX_x,
                                     geosxData idNdX_y,geosxData idNdX_z,
                                     Index_type const * iconstitutiveMap, geosxData idevStressData,
                                     geosxData imeanStress, real64 ishearModulus, real64 ibulkModulus,
                                     const real64 * idetJ, geosxData iacc_x,
                                     geosxData iacc_y, geosxData iacc_z, constUpdate updateState_ptr)
{


#if defined(CUDA)
  RAJA::forall<Pol>
    (elemList,[=] RAJA_DEVICE (Index_type k){
#else
  RAJA::forall<Pol>
    (elemList,[=](Index_type k){
#endif

      real64 uhat_local_x[inumNodesPerElement];
      real64 uhat_local_y[inumNodesPerElement];
      real64 uhat_local_z[inumNodesPerElement];

      real64 u_local_x[inumNodesPerElement];
      real64 u_local_y[inumNodesPerElement];
      real64 u_local_z[inumNodesPerElement];
      
      real64 f_local_x[inumNodesPerElement]; for(int i=0; i<inumNodesPerElement; ++i) f_local_x[i] = 0;
      real64 f_local_y[inumNodesPerElement]; for(int i=0; i<inumNodesPerElement; ++i) f_local_y[i] = 0;
      real64 f_local_z[inumNodesPerElement]; for(int i=0; i<inumNodesPerElement; ++i) f_local_z[i] = 0;
      
      //Gather element to node list
      Index_type nodeList[8];
      structuredElemToNodes(nodeList,k,nx,ny,nz);

      //Copy Global to Local
      for(Index_type a=0; a<inumNodesPerElement; ++a)
        {
          Index_type id = nodeList[a];
          u_local_x[a] = iu_x[id];
          u_local_y[a] = iu_y[id];
          u_local_z[a] = iu_z[id];
          
          uhat_local_x[a] = iuhat_x[id];
          uhat_local_y[a] = iuhat_y[id];
          uhat_local_z[a] = iuhat_z[id];
        }

      
      //Compute Quadrature
      for(Index_type q=0; q<inumQuadraturePoints; ++q)
        {


          real64 dUhatdX[local_dim][local_dim];
          real64 dUdX[local_dim][local_dim];

          for(Index_type ty=0; ty<local_dim; ++ty){
            for(Index_type tx=0; tx<local_dim; ++tx){
              dUhatdX[ty][tx] = 0.0;
              dUdX[ty][tx] = 0.0;
            }
          }


          for(Index_type a=0; a<inumNodesPerElement; ++a)
            {
              dUdX[0][0] += u_local_x[a]*idNdX_x(k,q,a);
              dUdX[0][1] += u_local_x[a]*idNdX_y(k,q,a);
              dUdX[0][2] += u_local_x[a]*idNdX_z(k,q,a);
              
              dUdX[1][0] += u_local_y[a]*idNdX_x(k,q,a);
              dUdX[1][1] += u_local_y[a]*idNdX_y(k,q,a);
              dUdX[1][2] += u_local_y[a]*idNdX_z(k,q,a);
              
              dUdX[2][0] += u_local_z[a]*idNdX_x(k,q,a);
              dUdX[2][1] += u_local_z[a]*idNdX_y(k,q,a);
              dUdX[2][2] += u_local_z[a]*idNdX_z(k,q,a);
            }

          for(Index_type a=0; a<inumNodesPerElement; ++a)
            {
              dUhatdX[0][0] += uhat_local_x[a]*idNdX_x(k,q,a);
              dUhatdX[0][1] += uhat_local_x[a]*idNdX_y(k,q,a);
              dUhatdX[0][2] += uhat_local_x[a]*idNdX_z(k,q,a);
              
              dUhatdX[1][0] += uhat_local_y[a]*idNdX_x(k,q,a);
              dUhatdX[1][1] += uhat_local_y[a]*idNdX_y(k,q,a);
              dUhatdX[1][2] += uhat_local_y[a]*idNdX_z(k,q,a);
              
              dUhatdX[2][0] += uhat_local_z[a]*idNdX_x(k,q,a);
              dUhatdX[2][1] += uhat_local_z[a]*idNdX_y(k,q,a);
              dUhatdX[2][2] += uhat_local_z[a]*idNdX_z(k,q,a);
            }
          

          real64 F[local_dim][local_dim];
          real64 Finv[local_dim][local_dim];
          real64 L[local_dim][local_dim];

          {            
            real64 dvdX[local_dim][local_dim];

            
            for(Index_type ty=0; ty<local_dim; ++ty){
              for(Index_type tx=0; tx<local_dim; ++tx){
                dvdX[ty][tx] = dUhatdX[ty][tx]*(1.0/dt);
              }
            }
            
            //calculate du/dX
            for(Index_type row=0; row<local_dim; ++row){
              for(Index_type col=0; col<local_dim; ++col){
                F[row][col] = dUhatdX[row][col];
              }
            }
            
            for(Index_type row=0; row<local_dim; ++row){
              for(Index_type col=0; col<local_dim; ++col){
                      F[row][col] *= 0.5;
              }
            }
            
            for(Index_type row=0; row<local_dim; ++row){
              for(Index_type col=0; col<local_dim; ++col){
                F[row][col] += dUdX[row][col];
              }
            }
            
            for(Index_type tx=0; tx<local_dim; ++tx){
              F[tx][tx] += 1.0;
            }

            
            Finverse(F, Finv);
            

            AijBjk(dvdX,Finv,L);
          }

          //Calculate gradient (end of step)
          for(Index_type ty=0; ty<local_dim; ++ty){
            for(Index_type tx=0; tx<local_dim; ++tx){
              F[ty][tx] = dUhatdX[ty][tx] + dUdX[ty][tx];
            }
          }

          for(int tx=0; tx<local_dim; ++tx){ F[tx][tx] += 1.0;};
          
          real64 detF = det<real64>(F);
          Finverse<real64> (F, Finv);

          real64 Rot[local_dim][local_dim];
          real64 Dadt[local_dim][local_dim];


          //-------------[Hughues Winget]--------------
          HughesWinget(Rot, Dadt, L, dt);
          //-------------------------------------------

          //-------------[Constitutive update]-------------
          Index_type m = iconstitutiveMap(k,q);          
          
          //UpdateStatePoint(Dadt,Rot,m, q,k, idevStressData, imeanStress,ishearModulus, ibulkModulus, NoElem);
          updateState_ptr(Dadt,Rot,m, q,k, idevStressData, imeanStress, ishearModulus, ibulkModulus, NoElem);
          //-------------------------------------------          


          real64 TotalStress[local_dim][local_dim];

          TotalStress[0][0] = idevStressData(k,q,0);
          TotalStress[1][0] = idevStressData(k,q,1);
          TotalStress[1][1] = idevStressData(k,q,2);
          TotalStress[2][0] = idevStressData(k,q,3);
          TotalStress[2][1] = idevStressData(k,q,4);
          TotalStress[2][2] = idevStressData(k,q,5);

          TotalStress[0][1] = TotalStress[1][0];
          TotalStress[0][2] = TotalStress[2][0];
          TotalStress[1][2] = TotalStress[2][1];

          for(Index_type i=0; i<local_dim; ++i)
            {
              TotalStress[i][i] += imeanStress[m];
            }

          //---------[Integrate - Function]---------------------
          real64 const integrationFactor = idetJ(k,q)*detF;
          real64 P[local_dim][local_dim];

          AijBkj(TotalStress,Finv,P);
          for(Index_type ty=0; ty<local_dim; ++ty){
            for(Index_type tx=0; tx<local_dim; ++tx){
              P[ty][tx] *= integrationFactor;
            }
          } 


          for(int a=0; a<inumNodesPerElement; ++a){
            //long int id = a + inumNodesPerElement*(q + inumQuadraturePoints*k);
            f_local_x[a] -= P[0][0]*idNdX_x(k,q,a) + P[0][1]*idNdX_y(k,q,a) + P[0][2]*idNdX_z(k,q,a);
            f_local_y[a] -= P[1][0]*idNdX_x(k,q,a) + P[1][1]*idNdX_y(k,q,a) + P[1][2]*idNdX_z(k,q,a);
            f_local_z[a] -= P[2][0]*idNdX_x(k,q,a) + P[2][1]*idNdX_y(k,q,a) + P[2][2]*idNdX_z(k,q,a);
          }

        }//end of quadrature


      //-------------[Add Local To Global]----------------
      for(Index_type a=0; a<inumNodesPerElement; ++a)
        {
          Index_type id = nodeList[a];
          RAJA::atomic::atomicAdd<atomicPol>(&iacc_x[id],f_local_x[a]);
          RAJA::atomic::atomicAdd<atomicPol>(&iacc_y[id],f_local_y[a]);
          RAJA::atomic::atomicAdd<atomicPol>(&iacc_z[id],f_local_z[a]);
        }          
      //---------------------------------------------------

    });

} 

///
///Solid mechanics update kernel with nodal degrees of freedom and
///shape function derivatives stored as an object of structs format.
///All computations are done in a monolithic kernel.
///      
template<typename Pol>
RAJA_INLINE void ArrayOfObjectsKernel(Index_type NoElem, RAJA::TypedListSegment<Index_type> elemList, real64 dt,
                                      const Index_type * elemsToNodes, geosxData iu,
                                      geosxData iuhat, geosxData idNdX,
                                      Index_type const * iconstitutiveMap, geosxData idevStressData,
                                      geosxData imeanStress, real64 ishearModulus, real64 ibulkModulus,
                                      const real64 * idetJ, geosxData iacc, constUpdate updateState_ptr)
{


#if defined(CUDA)
  RAJA::forall<Pol>
    (elemList,[=] RAJA_DEVICE (Index_type k){
#else
   RAJA::forall<Pol>
     (elemList,[=](Index_type k){
#endif
       
       real64 uhat_local[local_dim*inumNodesPerElement];
       real64 u_local[local_dim*inumNodesPerElement];
       real64 f_local[local_dim*inumNodesPerElement]; for(int i=0; i<local_dim*inumNodesPerElement; ++i) f_local[i] = 0;

              
      //Copy Global to Local
      for(Index_type a=0; a<inumNodesPerElement; ++a)
        {
          Index_type id = elemsToNodes[a + inumNodesPerElement*k];
          for(Index_type i=0; i<local_dim; ++i)
            {
              u_local[i + local_dim*a] = iu[i+local_dim*id];
              uhat_local[i + local_dim*a] = iuhat[i + local_dim*id];
            }
        }


      //Compute Quadrature
      for(Index_type q=0; q<inumQuadraturePoints; ++q)
        {
          real64 dUhatdX[local_dim][local_dim];
          real64 dUdX[local_dim][local_dim];

          for(Index_type ty=0; ty<local_dim; ++ty){
            for(Index_type tx=0; tx<local_dim; ++tx){
              dUhatdX[ty][tx] = 0.0;
              dUdX[ty][tx] = 0.0;
            }
          }

          
          for(Index_type a=0; a<inumNodesPerElement; ++a)
            {              
              dUdX[0][0] += u_local[0 + local_dim*a]*idNdX(k,q,a,0);
              dUdX[0][1] += u_local[0 + local_dim*a]*idNdX(k,q,a,1);
              dUdX[0][2] += u_local[0 + local_dim*a]*idNdX(k,q,a,2);
              
              dUdX[1][0] += u_local[1 + local_dim*a]*idNdX(k,q,a,0);
              dUdX[1][1] += u_local[1 + local_dim*a]*idNdX(k,q,a,1);
              dUdX[1][2] += u_local[1 + local_dim*a]*idNdX(k,q,a,2);
              
              dUdX[2][0] += u_local[2 + local_dim*a]*idNdX(k,q,a,0);
              dUdX[2][1] += u_local[2 + local_dim*a]*idNdX(k,q,a,1);
              dUdX[2][2] += u_local[2 + local_dim*a]*idNdX(k,q,a,2);
            }

          for(Index_type a=0; a<inumNodesPerElement; ++a)
            {              
              dUhatdX[0][0] += uhat_local[0 + local_dim*a]*idNdX(k,q,a,0);
              dUhatdX[0][1] += uhat_local[0 + local_dim*a]*idNdX(k,q,a,1);
              dUhatdX[0][2] += uhat_local[0 + local_dim*a]*idNdX(k,q,a,2);
              
              dUhatdX[1][0] += uhat_local[1 + local_dim*a]*idNdX(k,q,a,0);
              dUhatdX[1][1] += uhat_local[1 + local_dim*a]*idNdX(k,q,a,1);
              dUhatdX[1][2] += uhat_local[1 + local_dim*a]*idNdX(k,q,a,2);
              
              dUhatdX[2][0] += uhat_local[2 + local_dim*a]*idNdX(k,q,a,0);
              dUhatdX[2][1] += uhat_local[2 + local_dim*a]*idNdX(k,q,a,1);
              dUhatdX[2][2] += uhat_local[2 + local_dim*a]*idNdX(k,q,a,2);
            }
          

          real64 F[local_dim][local_dim];
          real64 Finv[local_dim][local_dim];
          real64 L[local_dim][local_dim];

          //Compute L
          {            
            real64 dvdX[local_dim][local_dim];

            
            for(Index_type ty=0; ty<local_dim; ++ty){
              for(Index_type tx=0; tx<local_dim; ++tx){
                dvdX[ty][tx] = dUhatdX[ty][tx]*(1.0/dt);
              }
            }
            
            //calculate du/dX
            for(Index_type row=0; row<local_dim; ++row){
              for(Index_type col=0; col<local_dim; ++col){
                F[row][col] = dUhatdX[row][col];
              }
            }
            
            for(Index_type row=0; row<local_dim; ++row){
              for(Index_type col=0; col<local_dim; ++col){
                      F[row][col] *= 0.5;
              }
            }
            
            for(Index_type row=0; row<local_dim; ++row){
              for(Index_type col=0; col<local_dim; ++col){
                F[row][col] += dUdX[row][col];
              }
            }
            
            for(Index_type tx=0; tx<local_dim; ++tx){
              F[tx][tx] += 1.0;
            }
            
            Finverse(F, Finv);
            

            AijBjk(dvdX,Finv,L);
          }

          //Calculate gradient (end of step)
          for(Index_type ty=0; ty<local_dim; ++ty){
            for(Index_type tx=0; tx<local_dim; ++tx){
              F[ty][tx] = dUhatdX[ty][tx] + dUdX[ty][tx];
            }
          }

          for(int tx=0; tx<local_dim; ++tx){ F[tx][tx] += 1.0;};
          
          real64 detF = det<real64>(F);
          Finverse<real64> (F, Finv);


          real64 Rot[local_dim][local_dim];
          real64 Dadt[local_dim][local_dim];

          //-------------[Hughues Winget]--------------
          HughesWinget(Rot, Dadt, L, dt);
          //-------------------------------------------
          

          //-------------[Constitutive update]-------------
          Index_type m = iconstitutiveMap(k,q);
          //UpdateStatePoint(Dadt,Rot,m, q,k, idevStressData, imeanStress,ishearModulus, ibulkModulus,NoElem);
          updateState_ptr(Dadt,Rot,m, q,k, idevStressData, imeanStress, ishearModulus, ibulkModulus, NoElem);
          //-------------------------------------------

          
          real64 TotalStress[local_dim][local_dim];          
          TotalStress[0][0] = idevStressData(k,q,0);
          TotalStress[1][0] = idevStressData(k,q,1);
          TotalStress[1][1] = idevStressData(k,q,2);
          TotalStress[2][0] = idevStressData(k,q,3);
          TotalStress[2][1] = idevStressData(k,q,4);
          TotalStress[2][2] = idevStressData(k,q,5);

          TotalStress[0][1] = TotalStress[1][0];
          TotalStress[0][2] = TotalStress[2][0];
          TotalStress[1][2] = TotalStress[2][1];

          for(Index_type i=0; i<local_dim; ++i)
            {
              TotalStress[i][i] += imeanStress[m];
            }
          

          //---------[Integrate - Function]---------------------
          real64 const integrationFactor = idetJ(k,q)*detF;
          real64 P[local_dim][local_dim];

          AijBkj(TotalStress,Finv,P);
          for(Index_type ty=0; ty<local_dim; ++ty){
            for(Index_type tx=0; tx<local_dim; ++tx){
              P[ty][tx] *= integrationFactor;
            }
          } 
          //--------------------------------------------------
          for(int a=0; a<inumNodesPerElement; ++a){
            
            f_local[0 + local_dim*a] -= P[0][0]*idNdX(k,q,a,0) + P[0][1]*idNdX(k,q,a,1) + P[0][2]*idNdX(k,q,a,2);
            f_local[1 + local_dim*a] -= P[1][0]*idNdX(k,q,a,0) + P[1][1]*idNdX(k,q,a,1) + P[1][2]*idNdX(k,q,a,2);
            f_local[2 + local_dim*a] -= P[2][0]*idNdX(k,q,a,0) + P[2][1]*idNdX(k,q,a,1) + P[2][2]*idNdX(k,q,a,2);
          }

        }//end of quadrature

      
      //-------------[Add Local To Global]----------------
      for(Index_type a=0; a<inumNodesPerElement; ++a)
        {          
          Index_type id = elemsToNodes[a + inumNodesPerElement*k];
          RAJA::atomic::atomicAdd<atomicPol>(&iacc[0 + local_dim*id],f_local[0 + local_dim*a]);
          RAJA::atomic::atomicAdd<atomicPol>(&iacc[1 + local_dim*id],f_local[1 + local_dim*a]);
          RAJA::atomic::atomicAdd<atomicPol>(&iacc[2 + local_dim*id],f_local[2 + local_dim*a]);
        }
      
      //---------------------------------------------------
     });
   
}

///
///Solid mechanics update kernel with nodal degrees of freedom and
///shape function derivatives stored as an object of structs format.
///All computations are done in a monolithic kernel. Domain is assumed
///to be a structured mesh.      
///            
template<typename Pol>
RAJA_INLINE void ArrayOfObjectsKernel_Structured(Index_type NoElem, RAJA::TypedListSegment<Index_type> elemList, real64 dt,
                                                 Index_type nx, Index_type ny, Index_type nz,
                                                 geosxData iu, geosxData iuhat, geosxData idNdX,                                    
                                                 Index_type const * iconstitutiveMap, geosxData idevStressData,
                                                 geosxData imeanStress, real64 ishearModulus, real64 ibulkModulus,
                                                 const real64 * idetJ, geosxData iacc, constUpdate updateState_ptr)
{


#if defined(CUDA)
  RAJA::forall<Pol>
    (elemList,[=] RAJA_DEVICE (Index_type k){
#else
   RAJA::forall<Pol>
     (elemList,[=](Index_type k){
#endif
       
       real64 uhat_local[local_dim*inumNodesPerElement];
       real64 u_local[local_dim*inumNodesPerElement];
       real64 f_local[local_dim*inumNodesPerElement]; for(int i=0; i<local_dim*inumNodesPerElement; ++i) f_local[i] = 0;

      //Gather element to node list
      Index_type nodeList[8];
      structuredElemToNodes(nodeList,k,nx,ny,nz);
              
      //Copy Global to Local
      for(Index_type a=0; a<inumNodesPerElement; ++a)
        {
          Index_type id = nodeList[a];
          for(Index_type i=0; i<local_dim; ++i)
            {
              u_local[i + local_dim*a] = iu[i+local_dim*id];
              uhat_local[i + local_dim*a] = iuhat[i + local_dim*id];
            }
        }


      //Compute Quadrature
      for(Index_type q=0; q<inumQuadraturePoints; ++q)
        {
          real64 dUhatdX[local_dim][local_dim];
          real64 dUdX[local_dim][local_dim];

          for(Index_type ty=0; ty<local_dim; ++ty){
            for(Index_type tx=0; tx<local_dim; ++tx){
              dUhatdX[ty][tx] = 0.0;
              dUdX[ty][tx] = 0.0;
            }
          }

          
          for(Index_type a=0; a<inumNodesPerElement; ++a)
            {              
              dUdX[0][0] += u_local[0 + local_dim*a]*idNdX(k,q,a,0);
              dUdX[0][1] += u_local[0 + local_dim*a]*idNdX(k,q,a,1);
              dUdX[0][2] += u_local[0 + local_dim*a]*idNdX(k,q,a,2);
              
              dUdX[1][0] += u_local[1 + local_dim*a]*idNdX(k,q,a,0);
              dUdX[1][1] += u_local[1 + local_dim*a]*idNdX(k,q,a,1);
              dUdX[1][2] += u_local[1 + local_dim*a]*idNdX(k,q,a,2);
              
              dUdX[2][0] += u_local[2 + local_dim*a]*idNdX(k,q,a,0);
              dUdX[2][1] += u_local[2 + local_dim*a]*idNdX(k,q,a,1);
              dUdX[2][2] += u_local[2 + local_dim*a]*idNdX(k,q,a,2);
            }

          for(Index_type a=0; a<inumNodesPerElement; ++a)
            {              
              dUhatdX[0][0] += uhat_local[0 + local_dim*a]*idNdX(k,q,a,0);
              dUhatdX[0][1] += uhat_local[0 + local_dim*a]*idNdX(k,q,a,1);
              dUhatdX[0][2] += uhat_local[0 + local_dim*a]*idNdX(k,q,a,2);
              
              dUhatdX[1][0] += uhat_local[1 + local_dim*a]*idNdX(k,q,a,0);
              dUhatdX[1][1] += uhat_local[1 + local_dim*a]*idNdX(k,q,a,1);
              dUhatdX[1][2] += uhat_local[1 + local_dim*a]*idNdX(k,q,a,2);
              
              dUhatdX[2][0] += uhat_local[2 + local_dim*a]*idNdX(k,q,a,0);
              dUhatdX[2][1] += uhat_local[2 + local_dim*a]*idNdX(k,q,a,1);
              dUhatdX[2][2] += uhat_local[2 + local_dim*a]*idNdX(k,q,a,2);
            }
          

          real64 F[local_dim][local_dim];
          real64 Finv[local_dim][local_dim];
          real64 L[local_dim][local_dim];

          //Compute L
          {            
            real64 dvdX[local_dim][local_dim];

            
            for(Index_type ty=0; ty<local_dim; ++ty){
              for(Index_type tx=0; tx<local_dim; ++tx){
                dvdX[ty][tx] = dUhatdX[ty][tx]*(1.0/dt);
              }
            }
            
            //calculate du/dX
            for(Index_type row=0; row<local_dim; ++row){
              for(Index_type col=0; col<local_dim; ++col){
                F[row][col] = dUhatdX[row][col];
              }
            }
            
            for(Index_type row=0; row<local_dim; ++row){
              for(Index_type col=0; col<local_dim; ++col){
                      F[row][col] *= 0.5;
              }
            }
            
            for(Index_type row=0; row<local_dim; ++row){
              for(Index_type col=0; col<local_dim; ++col){
                F[row][col] += dUdX[row][col];
              }
            }
            
            for(Index_type tx=0; tx<local_dim; ++tx){
              F[tx][tx] += 1.0;
            }
            
            Finverse(F, Finv);
            

            AijBjk(dvdX,Finv,L);
          }

          //Calculate gradient (end of step)
          for(Index_type ty=0; ty<local_dim; ++ty){
            for(Index_type tx=0; tx<local_dim; ++tx){
              F[ty][tx] = dUhatdX[ty][tx] + dUdX[ty][tx];
            }
          }

          for(int tx=0; tx<local_dim; ++tx){ F[tx][tx] += 1.0;};
          
          real64 detF = det<real64>(F);
          Finverse<real64> (F, Finv);


          real64 Rot[local_dim][local_dim];
          real64 Dadt[local_dim][local_dim];

          //-------------[Hughues Winget]--------------
          HughesWinget(Rot, Dadt, L, dt);
          //-------------------------------------------
          

          //-------------[Constitutive update]-------------
          Index_type m = iconstitutiveMap(k,q);
          //UpdateStatePoint(Dadt,Rot,m, q,k, idevStressData, imeanStress,ishearModulus, ibulkModulus,NoElem);
          updateState_ptr(Dadt,Rot,m, q,k, idevStressData, imeanStress, ishearModulus, ibulkModulus, NoElem);
          //-------------------------------------------

          
          real64 TotalStress[local_dim][local_dim];          
          TotalStress[0][0] = idevStressData(k,q,0);
          TotalStress[1][0] = idevStressData(k,q,1);
          TotalStress[1][1] = idevStressData(k,q,2);
          TotalStress[2][0] = idevStressData(k,q,3);
          TotalStress[2][1] = idevStressData(k,q,4);
          TotalStress[2][2] = idevStressData(k,q,5);

          TotalStress[0][1] = TotalStress[1][0];
          TotalStress[0][2] = TotalStress[2][0];
          TotalStress[1][2] = TotalStress[2][1];

          for(Index_type i=0; i<local_dim; ++i)
            {
              TotalStress[i][i] += imeanStress[m];
            }
          

          //---------[Integrate - Function]---------------------
          real64 const integrationFactor = idetJ(k,q)*detF;
          real64 P[local_dim][local_dim];

          AijBkj(TotalStress,Finv,P);
          for(Index_type ty=0; ty<local_dim; ++ty){
            for(Index_type tx=0; tx<local_dim; ++tx){
              P[ty][tx] *= integrationFactor;
            }
          } 
          //--------------------------------------------------
          for(int a=0; a<inumNodesPerElement; ++a){
            
            f_local[0 + local_dim*a] -= P[0][0]*idNdX(k,q,a,0) + P[0][1]*idNdX(k,q,a,1) + P[0][2]*idNdX(k,q,a,2);
            f_local[1 + local_dim*a] -= P[1][0]*idNdX(k,q,a,0) + P[1][1]*idNdX(k,q,a,1) + P[1][2]*idNdX(k,q,a,2);
            f_local[2 + local_dim*a] -= P[2][0]*idNdX(k,q,a,0) + P[2][1]*idNdX(k,q,a,1) + P[2][2]*idNdX(k,q,a,2);
          }

        }//end of quadrature

      
      //-------------[Add Local To Global]----------------
      for(Index_type a=0; a<inumNodesPerElement; ++a)
        {          
          Index_type id = nodeList[a];
          RAJA::atomic::atomicAdd<atomicPol>(&iacc[0 + local_dim*id],f_local[0 + local_dim*a]);
          RAJA::atomic::atomicAdd<atomicPol>(&iacc[1 + local_dim*id],f_local[1 + local_dim*a]);
          RAJA::atomic::atomicAdd<atomicPol>(&iacc[2 + local_dim*id],f_local[2 + local_dim*a]);
        }
      
      //---------------------------------------------------
     });
   
}

///
///Solid mechanics update kernel with nodal degrees of freedom and
///shape function derivatives stored as an object of structs format.
///All computations are done in a monolithic kernel. Only the kinematic step
///is taken here. 
///            
template<typename Pol>
RAJA_INLINE void ArrayOfObjects_KinematicKernel(Index_type NoElem, RAJA::TypedListSegment<Index_type> elemList, real64 dt,
                                                const Index_type * elemsToNodes,
                                                geosxData iu, geosxData iuhat, geosxData idNdX,                          
                                                Index_type const * iconstitutiveMap, geosxData idevStressData,
                                                geosxData imeanStress, real64 ishearModulus, real64 ibulkModulus,
                                                const real64 * idetJ, geosxData iacc, geosxData Dadt_ptr, geosxData Rot_ptr,
                                                geosxData detF_ptr, geosxData Finv_ptr)
{


#if defined(CUDA)
  RAJA::forall<Pol>
    (elemList,[=] RAJA_DEVICE (Index_type k){
#else
   RAJA::forall<Pol>
     (elemList,[=](Index_type k){
#endif
       
       real64 uhat_local[local_dim*inumNodesPerElement];
       real64 u_local[local_dim*inumNodesPerElement];
              
      //Copy Global to Local
      for(Index_type a=0; a<inumNodesPerElement; ++a)
        {
          Index_type id = elemsToNodes[a + inumNodesPerElement*k];
          for(Index_type i=0; i<local_dim; ++i)
            {
              u_local[i + local_dim*a] = iu[i+local_dim*id];
              uhat_local[i + local_dim*a] = iuhat[i + local_dim*id];
            }
        }
              

      //Compute Quadrature
      for(Index_type q=0; q<inumQuadraturePoints; ++q)
        {
          real64 dUhatdX[local_dim][local_dim];
          real64 dUdX[local_dim][local_dim];

          for(Index_type ty=0; ty<local_dim; ++ty){
            for(Index_type tx=0; tx<local_dim; ++tx){
              dUhatdX[ty][tx] = 0.0;
              dUdX[ty][tx] = 0.0;
            }
          }

          
          for(Index_type a=0; a<inumNodesPerElement; ++a)
            {              
              dUdX[0][0] += u_local[0 + local_dim*a]*idNdX(k,q,a,0);
              dUdX[0][1] += u_local[0 + local_dim*a]*idNdX(k,q,a,1);
              dUdX[0][2] += u_local[0 + local_dim*a]*idNdX(k,q,a,2);
              
              dUdX[1][0] += u_local[1 + local_dim*a]*idNdX(k,q,a,0);
              dUdX[1][1] += u_local[1 + local_dim*a]*idNdX(k,q,a,1);
              dUdX[1][2] += u_local[1 + local_dim*a]*idNdX(k,q,a,2);
              
              dUdX[2][0] += u_local[2 + local_dim*a]*idNdX(k,q,a,0);
              dUdX[2][1] += u_local[2 + local_dim*a]*idNdX(k,q,a,1);
              dUdX[2][2] += u_local[2 + local_dim*a]*idNdX(k,q,a,2);
            }

          for(Index_type a=0; a<inumNodesPerElement; ++a)
            {              
              dUhatdX[0][0] += uhat_local[0 + local_dim*a]*idNdX(k,q,a,0);
              dUhatdX[0][1] += uhat_local[0 + local_dim*a]*idNdX(k,q,a,1);
              dUhatdX[0][2] += uhat_local[0 + local_dim*a]*idNdX(k,q,a,2);
              
              dUhatdX[1][0] += uhat_local[1 + local_dim*a]*idNdX(k,q,a,0);
              dUhatdX[1][1] += uhat_local[1 + local_dim*a]*idNdX(k,q,a,1);
              dUhatdX[1][2] += uhat_local[1 + local_dim*a]*idNdX(k,q,a,2);
              
              dUhatdX[2][0] += uhat_local[2 + local_dim*a]*idNdX(k,q,a,0);
              dUhatdX[2][1] += uhat_local[2 + local_dim*a]*idNdX(k,q,a,1);
              dUhatdX[2][2] += uhat_local[2 + local_dim*a]*idNdX(k,q,a,2);
            }
          

          real64 F[local_dim][local_dim];
          real64 Finv[local_dim][local_dim];
          real64 L[local_dim][local_dim];

          //Compute L
          {            
            real64 dvdX[local_dim][local_dim];

            
            for(Index_type ty=0; ty<local_dim; ++ty){
              for(Index_type tx=0; tx<local_dim; ++tx){
                dvdX[ty][tx] = dUhatdX[ty][tx]*(1.0/dt);
              }
            }
            
            //calculate du/dX
            for(Index_type row=0; row<local_dim; ++row){
              for(Index_type col=0; col<local_dim; ++col){
                F[row][col] = dUhatdX[row][col];
              }
            }
            
            for(Index_type row=0; row<local_dim; ++row){
              for(Index_type col=0; col<local_dim; ++col){
                F[row][col] *= 0.5;
              }
            }
            
            for(Index_type row=0; row<local_dim; ++row){
              for(Index_type col=0; col<local_dim; ++col){
                F[row][col] += dUdX[row][col];
              }
            }
            
            for(Index_type tx=0; tx<local_dim; ++tx){
              F[tx][tx] += 1.0;
            }
            
            Finverse(F, Finv);
            

            AijBjk(dvdX,Finv,L);
          }

          //Calculate gradient (end of step)
          for(Index_type ty=0; ty<local_dim; ++ty){
            for(Index_type tx=0; tx<local_dim; ++tx){
              F[ty][tx] = dUhatdX[ty][tx] + dUdX[ty][tx];
            }
          }

          for(int tx=0; tx<local_dim; ++tx){ F[tx][tx] += 1.0;};
          
          real64 detF = det<real64>(F);
          Finverse<real64> (F, Finv);


          real64 Rot[local_dim][local_dim];
          real64 Dadt[local_dim][local_dim];

          //-------------[Hughues Winget]--------------
          HughesWinget(Rot, Dadt, L, dt);
          //-------------------------------------------

          //Write out data to global memory
          detF_ptr(k,q) = detF;
          for(Index_type r = 0; r < local_dim; ++r){
            for(Index_type c = 0; c < local_dim; ++c){
                Dadt_ptr(k,q,r,c) = Dadt[r][c];
                Rot_ptr(k,q,r,c) = Rot[r][c];
                Finv_ptr(k,q,r,c) = Finv[r][c];
            }
          }
                            
        }//end of quadrature

     });

   
}      
///
///Solid mechanics update kernel with nodal degrees of freedom and
///shape function derivatives stored as an array of objects format.
///All computations are done in a monolithic kernel. Only the kinematic step
///is taken here. Domain is assumed to be structured. 
///            
template<typename Pol>
RAJA_INLINE void ArrayOfObjects_KinematicKernel_Structured(Index_type NoElem, RAJA::TypedListSegment<Index_type> elemList, real64 dt,
                                                 Index_type nx, Index_type ny, Index_type nz,
                                                 geosxData iu, geosxData iuhat, geosxData idNdX,                          
                                                 Index_type const * iconstitutiveMap, geosxData idevStressData,
                                                 geosxData imeanStress, real64 ishearModulus, real64 ibulkModulus,
                                                 const real64 * idetJ, geosxData iacc, geosxData Dadt_ptr, geosxData Rot_ptr,
                                                 geosxData detF_ptr, geosxData Finv_ptr)
{


#if defined(CUDA)
  RAJA::forall<Pol>
    (elemList,[=] RAJA_DEVICE (Index_type k){
#else
   RAJA::forall<Pol>
     (elemList,[=](Index_type k){
#endif
       
       real64 uhat_local[local_dim*inumNodesPerElement];
       real64 u_local[local_dim*inumNodesPerElement];


      //Gather element to node list
      Index_type nodeList[8];
      structuredElemToNodes(nodeList,k,nx,ny,nz);
              
      //Copy Global to Local
      for(Index_type a=0; a<inumNodesPerElement; ++a)
        {
          Index_type id = nodeList[a];
          for(Index_type i=0; i<local_dim; ++i)
            {
              u_local[i + local_dim*a] = iu[i+local_dim*id];
              uhat_local[i + local_dim*a] = iuhat[i + local_dim*id];
            }
        }
              

      //Compute Quadrature
      for(Index_type q=0; q<inumQuadraturePoints; ++q)
        {
          real64 dUhatdX[local_dim][local_dim];
          real64 dUdX[local_dim][local_dim];

          for(Index_type ty=0; ty<local_dim; ++ty){
            for(Index_type tx=0; tx<local_dim; ++tx){
              dUhatdX[ty][tx] = 0.0;
              dUdX[ty][tx] = 0.0;
            }
          }

          
          for(Index_type a=0; a<inumNodesPerElement; ++a)
            {              
              dUdX[0][0] += u_local[0 + local_dim*a]*idNdX(k,q,a,0);
              dUdX[0][1] += u_local[0 + local_dim*a]*idNdX(k,q,a,1);
              dUdX[0][2] += u_local[0 + local_dim*a]*idNdX(k,q,a,2);
              
              dUdX[1][0] += u_local[1 + local_dim*a]*idNdX(k,q,a,0);
              dUdX[1][1] += u_local[1 + local_dim*a]*idNdX(k,q,a,1);
              dUdX[1][2] += u_local[1 + local_dim*a]*idNdX(k,q,a,2);
              
              dUdX[2][0] += u_local[2 + local_dim*a]*idNdX(k,q,a,0);
              dUdX[2][1] += u_local[2 + local_dim*a]*idNdX(k,q,a,1);
              dUdX[2][2] += u_local[2 + local_dim*a]*idNdX(k,q,a,2);
            }

          for(Index_type a=0; a<inumNodesPerElement; ++a)
            {              
              dUhatdX[0][0] += uhat_local[0 + local_dim*a]*idNdX(k,q,a,0);
              dUhatdX[0][1] += uhat_local[0 + local_dim*a]*idNdX(k,q,a,1);
              dUhatdX[0][2] += uhat_local[0 + local_dim*a]*idNdX(k,q,a,2);
              
              dUhatdX[1][0] += uhat_local[1 + local_dim*a]*idNdX(k,q,a,0);
              dUhatdX[1][1] += uhat_local[1 + local_dim*a]*idNdX(k,q,a,1);
              dUhatdX[1][2] += uhat_local[1 + local_dim*a]*idNdX(k,q,a,2);
              
              dUhatdX[2][0] += uhat_local[2 + local_dim*a]*idNdX(k,q,a,0);
              dUhatdX[2][1] += uhat_local[2 + local_dim*a]*idNdX(k,q,a,1);
              dUhatdX[2][2] += uhat_local[2 + local_dim*a]*idNdX(k,q,a,2);
            }
          

          real64 F[local_dim][local_dim];
          real64 Finv[local_dim][local_dim];
          real64 L[local_dim][local_dim];

          //Compute L
          {            
            real64 dvdX[local_dim][local_dim];

            
            for(Index_type ty=0; ty<local_dim; ++ty){
              for(Index_type tx=0; tx<local_dim; ++tx){
                dvdX[ty][tx] = dUhatdX[ty][tx]*(1.0/dt);
              }
            }
            
            //calculate du/dX
            for(Index_type row=0; row<local_dim; ++row){
              for(Index_type col=0; col<local_dim; ++col){
                F[row][col] = dUhatdX[row][col];
              }
            }
            
            for(Index_type row=0; row<local_dim; ++row){
              for(Index_type col=0; col<local_dim; ++col){
                F[row][col] *= 0.5;
              }
            }
            
            for(Index_type row=0; row<local_dim; ++row){
              for(Index_type col=0; col<local_dim; ++col){
                F[row][col] += dUdX[row][col];
              }
            }
            
            for(Index_type tx=0; tx<local_dim; ++tx){
              F[tx][tx] += 1.0;
            }
            
            Finverse(F, Finv);
            

            AijBjk(dvdX,Finv,L);
          }

          //Calculate gradient (end of step)
          for(Index_type ty=0; ty<local_dim; ++ty){
            for(Index_type tx=0; tx<local_dim; ++tx){
              F[ty][tx] = dUhatdX[ty][tx] + dUdX[ty][tx];
            }
          }

          for(int tx=0; tx<local_dim; ++tx){ F[tx][tx] += 1.0;};
          
          real64 detF = det<real64>(F);
          Finverse<real64> (F, Finv);


          real64 Rot[local_dim][local_dim];
          real64 Dadt[local_dim][local_dim];

          //-------------[Hughues Winget]--------------
          HughesWinget(Rot, Dadt, L, dt);
          //-------------------------------------------

          //Write out data to global memory
          detF_ptr(k,q) = detF;
          for(Index_type r = 0; r < local_dim; ++r){
            for(Index_type c = 0; c < local_dim; ++c){
                Dadt_ptr(k,q,r,c) = Dadt[r][c];
                Rot_ptr(k,q,r,c) = Rot[r][c];
                Finv_ptr(k,q,r,c) = Finv[r][c];
            }
          }
                            
        }//end of quadrature

     });

   
}      

///
///Constitutive update. This would normally be a function pointer in a monolithic kernel.
///
template<typename Pol>
RAJA_INLINE void ConstitutiveUpdateKernel(Index_type NoElem, RAJA::TypedListSegment<Index_type> elemList,
                                          geosxData Dadt_ptr, geosxData Rot_ptr, Index_type const * iconstitutiveMap,
                                          geosxData idevStressData, geosxData imeanStress, real64 shearModulus, real64 bulkModulus)
                          
{


#if defined(CUDA)
  RAJA::forall<Pol>
    (elemList,[=] RAJA_DEVICE (Index_type k){
#else
  RAJA::forall<Pol>
    (elemList,[=](Index_type k){
#endif

      for(Index_type q=0; q < inumQuadraturePoints; ++q){
      

        Index_type m = iconstitutiveMap(k,q);
          
        real64 volumeStrain = Dadt_ptr(k,q,0,0) + Dadt_ptr(k,q,1,1) + Dadt_ptr(k,q,2,2);
        imeanStress[m] += volumeStrain * bulkModulus;
        
        real64 temp[local_dim][local_dim];
        for(Index_type i=0; i<3; ++i)
          {
            for(Index_type j=0; j<3; ++j)
              {
                temp[i][j] = Dadt_ptr(k,q,i,j);
            }
            temp[i][i] -= volumeStrain / 3.0;
          }
        
        for(Index_type ty=0; ty<3; ++ty)
          {
            for(Index_type tx=0; tx<3; ++tx)
              {
                temp[ty][tx] *= 2.0 * shearModulus;
              }
          }
        
        real64 localDevStress[local_dim][local_dim];
        idevStressData(k,q,0) += temp[0][0];
        idevStressData(k,q,1) += temp[1][0];
        idevStressData(k,q,2) += temp[1][1];
        idevStressData(k,q,3) += temp[2][0];
        idevStressData(k,q,4) += temp[2][1];
        idevStressData(k,q,5) += temp[2][2];
        
        localDevStress[0][0] = idevStressData(k,q,0);
        localDevStress[1][0] = idevStressData(k,q,1);
        localDevStress[1][1] = idevStressData(k,q,2);
        localDevStress[2][0] = idevStressData(k,q,3);
        localDevStress[2][1] = idevStressData(k,q,4);
        localDevStress[2][2] = idevStressData(k,q,5);
        
        localDevStress[0][1] = localDevStress[1][0];
        localDevStress[0][2] = localDevStress[2][0];
        localDevStress[1][2] = localDevStress[2][1];

        //Make a local copy
        real64 Rot[local_dim][local_dim];
        for(Index_type r=0; r<local_dim; ++r){
          for(Index_type c=0; c<local_dim; ++c){
            Rot[r][c] = Rot_ptr(k,q,r,c); 
          }
        }
        
        //QijAjkQlk
        AijBjk(Rot,localDevStress,temp);
        AijBkj(temp,Rot,localDevStress);
        
        
        idevStressData(k,q,0) = localDevStress[0][0];
        idevStressData(k,q,1) = localDevStress[1][0];
        idevStressData(k,q,2) = localDevStress[1][1];
        idevStressData(k,q,3) = localDevStress[2][0];
        idevStressData(k,q,4) = localDevStress[2][1];
        idevStressData(k,q,5) = localDevStress[2][2];
        
      }//quadrature loop
        
    }); //element loop
  
}


///
///Integration kernel, assumes nodal degrees of freedom and shape function derivatives are
///stored in a array of objects format.
///
template<typename Pol>
RAJA_INLINE void ArrayOfObjects_IntegrationKernel(Index_type NoElem, RAJA::TypedListSegment<Index_type> elemList, real64 dt,
                                                  const Index_type * elemsToNodes,
                                                  geosxData iu, geosxData iuhat, geosxData idNdX,
                                                  Index_type const * iconstitutiveMap, geosxData idevStressData,
                                                  geosxData imeanStress, real64 ishearModulus, real64 ibulkModulus,
                                                  const real64 * idetJ, geosxData iacc, geosxData Dadt_ptr, geosxData Rot_ptr,
                                                  geosxData detF_ptr, geosxData Finv_ptr)
{

  #if defined(CUDA)
  RAJA::forall<Pol>
    (elemList,[=] RAJA_DEVICE (Index_type k){
#else
   RAJA::forall<Pol>
     (elemList,[=](Index_type k){
#endif

       real64 f_local[local_dim*inumNodesPerElement]; for(int i=0; i<local_dim*inumNodesPerElement; ++i) f_local[i] = 0;
  
       for(Index_type q=0; q < inumQuadraturePoints; ++q){
  
         Index_type m = iconstitutiveMap(k,q);
         real64 TotalStress[local_dim][local_dim];
  
         TotalStress[0][0] = idevStressData(k,q,0);
         TotalStress[1][0] = idevStressData(k,q,1);
         TotalStress[1][1] = idevStressData(k,q,2);
         TotalStress[2][0] = idevStressData(k,q,3);
         TotalStress[2][1] = idevStressData(k,q,4);
         TotalStress[2][2] = idevStressData(k,q,5);
         
         TotalStress[0][1] = TotalStress[1][0];
         TotalStress[0][2] = TotalStress[2][0];
         TotalStress[1][2] = TotalStress[2][1];
         
          for(Index_type i=0; i<local_dim; ++i)
            {
              TotalStress[i][i] += imeanStress[m];
            }

          //---------[Integrate - Function]---------------------
          real64 const integrationFactor = idetJ(k,q)*detF_ptr(k,q);
          real64 P[local_dim][local_dim];
          real64 Finv[local_dim][local_dim];

          for(Index_type r=0; r<local_dim; ++r){
            for(Index_type c=0; c<local_dim; ++c){
               Finv[r][c] = Finv_ptr(k,q,r,c);
            }
          }
          
          AijBkj(TotalStress,Finv,P);
          for(Index_type ty=0; ty<local_dim; ++ty){
            for(Index_type tx=0; tx<local_dim; ++tx){
              P[ty][tx] *= integrationFactor;
            }
          } 
          
          //--------------------------------------------------
          for(int a=0; a<inumNodesPerElement; ++a){            
            f_local[0 + local_dim*a] -= P[0][0]*idNdX(k,q,a,0) + P[0][1]*idNdX(k,q,a,1) + P[0][2]*idNdX(k,q,a,2);
            f_local[1 + local_dim*a] -= P[1][0]*idNdX(k,q,a,0) + P[1][1]*idNdX(k,q,a,1) + P[1][2]*idNdX(k,q,a,2);
            f_local[2 + local_dim*a] -= P[2][0]*idNdX(k,q,a,0) + P[2][1]*idNdX(k,q,a,1) + P[2][2]*idNdX(k,q,a,2);
           }

         }//end of quadrature

      
       //-------------[Add Local To Global]----------------
       for(Index_type a=0; a<inumNodesPerElement; ++a)
         {          
           Index_type id = elemsToNodes[a + inumNodesPerElement*k];
           RAJA::atomic::atomicAdd<atomicPol>(&iacc[0 + local_dim*id],f_local[0 + local_dim*a]);
           RAJA::atomic::atomicAdd<atomicPol>(&iacc[1 + local_dim*id],f_local[1 + local_dim*a]);
           RAJA::atomic::atomicAdd<atomicPol>(&iacc[2 + local_dim*id],f_local[2 + local_dim*a]);
         }
                   
     });
      
}


///
///Integration kernel, assumes nodal degrees of freedom and shape function derivatives are
///stored in a array of objects format. Domain is assumed to be structured.
///
template<typename Pol>
RAJA_INLINE void ArrayOfObjects_IntegrationKernel_Structured(Index_type NoElem, RAJA::TypedListSegment<Index_type> elemList, real64 dt,
                                                             Index_type nx, Index_type ny, Index_type nz,
                                                             geosxData iu, geosxData iuhat, geosxData idNdX,
                                                             Index_type const * iconstitutiveMap, geosxData idevStressData,
                                                             geosxData imeanStress, real64 ishearModulus, real64 ibulkModulus,
                                                             const real64 * idetJ, geosxData iacc, geosxData Dadt_ptr, geosxData Rot_ptr,
                                                             geosxData detF_ptr, geosxData Finv_ptr)
{
  
  #if defined(CUDA)
  RAJA::forall<Pol>
    (elemList,[=] RAJA_DEVICE (Index_type k){
#else
   RAJA::forall<Pol>
     (elemList,[=](Index_type k){
#endif

       real64 f_local[local_dim*inumNodesPerElement]; for(int i=0; i<local_dim*inumNodesPerElement; ++i) f_local[i] = 0;

       //Gather element to node list
       Index_type nodeList[8];
       structuredElemToNodes(nodeList,k,nx,ny,nz);

       for(Index_type q=0; q < inumQuadraturePoints; ++q){
  
         Index_type m = iconstitutiveMap(k,q);
         real64 TotalStress[local_dim][local_dim];
  
         TotalStress[0][0] = idevStressData(k,q,0);
         TotalStress[1][0] = idevStressData(k,q,1);
         TotalStress[1][1] = idevStressData(k,q,2);
         TotalStress[2][0] = idevStressData(k,q,3);
         TotalStress[2][1] = idevStressData(k,q,4);
         TotalStress[2][2] = idevStressData(k,q,5);
         
         TotalStress[0][1] = TotalStress[1][0];
         TotalStress[0][2] = TotalStress[2][0];
         TotalStress[1][2] = TotalStress[2][1];
         
          for(Index_type i=0; i<local_dim; ++i)
            {
              TotalStress[i][i] += imeanStress[m];
            }

          //---------[Integrate - Function]---------------------
          real64 const integrationFactor = idetJ(k,q)*detF_ptr(k,q);
          real64 P[local_dim][local_dim];
          real64 Finv[local_dim][local_dim];

          for(Index_type r=0; r<local_dim; ++r){
            for(Index_type c=0; c<local_dim; ++c){
               Finv[r][c] = Finv_ptr(k,q,r,c);
            }
          }
          
          AijBkj(TotalStress,Finv,P);
          for(Index_type ty=0; ty<local_dim; ++ty){
            for(Index_type tx=0; tx<local_dim; ++tx){
              P[ty][tx] *= integrationFactor;
            }
          } 
          
          //--------------------------------------------------
          for(int a=0; a<inumNodesPerElement; ++a){            
            f_local[0 + local_dim*a] -= P[0][0]*idNdX(k,q,a,0) + P[0][1]*idNdX(k,q,a,1) + P[0][2]*idNdX(k,q,a,2);
            f_local[1 + local_dim*a] -= P[1][0]*idNdX(k,q,a,0) + P[1][1]*idNdX(k,q,a,1) + P[1][2]*idNdX(k,q,a,2);
            f_local[2 + local_dim*a] -= P[2][0]*idNdX(k,q,a,0) + P[2][1]*idNdX(k,q,a,1) + P[2][2]*idNdX(k,q,a,2);
           }

         }//end of quadrature

      
       //-------------[Add Local To Global]----------------
       for(Index_type a=0; a<inumNodesPerElement; ++a)
         {          
           Index_type id = nodeList[a];
           RAJA::atomic::atomicAdd<atomicPol>(&iacc[0 + local_dim*id],f_local[0 + local_dim*a]);
           RAJA::atomic::atomicAdd<atomicPol>(&iacc[1 + local_dim*id],f_local[1 + local_dim*a]);
           RAJA::atomic::atomicAdd<atomicPol>(&iacc[2 + local_dim*id],f_local[2 + local_dim*a]);
         }
                   
     });
      
}


///
///Solid mechanics kinematic kernel with nodal degrees of freedom stored and shape function
///derivivatives stored in an array of objects format. This kernel only performs the kinematic step.
///            
template<typename Pol>
RAJA_INLINE void ObjectOfArrays_KinematicKernel(Index_type NoElem, RAJA::TypedListSegment<Index_type> elemList, real64 dt,
                                                const Index_type * elemsToNodes,
                                                geosxData iu_x, geosxData iu_y, geosxData iu_z,
                                                geosxData iuhat_x, geosxData iuhat_y, geosxData iuhat_z,
                                                geosxData idNdX_x, geosxData idNdX_y, geosxData idNdX_z,
                                                Index_type const * iconstitutiveMap, geosxData idevStressData,
                                                geosxData imeanStress, real64 ishearModulus, real64 ibulkModulus,
                                                const real64 * idetJ,
                                                geosxData iacc_x, geosxData iacc_y, geosxData iacc_z,
                                                geosxData Dadt_ptr, geosxData Rot_ptr,
                                                geosxData detF_ptr, geosxData Finv_ptr)
{


#if defined(CUDA)
  RAJA::forall<Pol>
    (elemList,[=] RAJA_DEVICE (Index_type k){
#else
   RAJA::forall<Pol>
     (elemList,[=](Index_type k){
#endif
       
       real64 uhat_local_x[inumNodesPerElement];
       real64 uhat_local_y[inumNodesPerElement];
       real64 uhat_local_z[inumNodesPerElement];
       
       real64 u_local_x[inumNodesPerElement];
       real64 u_local_y[inumNodesPerElement];
       real64 u_local_z[inumNodesPerElement];
       
       //Copy Global to Local
       for(Index_type a=0; a<inumNodesPerElement; ++a)
         {
           Index_type id = elemsToNodes[a + inumNodesPerElement*k];

           u_local_x[a] = iu_x[id];
           u_local_y[a] = iu_y[id];
           u_local_z[a] = iu_z[id];
           
           uhat_local_x[a] = iuhat_x[id];
           uhat_local_y[a] = iuhat_y[id];
           uhat_local_z[a] = iuhat_z[id];
         }
       
       
       //Compute Quadrature
       for(Index_type q=0; q<inumQuadraturePoints; ++q)
         {           
           
           real64 dUhatdX[local_dim][local_dim];
           real64 dUdX[local_dim][local_dim];
           
           for(Index_type ty=0; ty<local_dim; ++ty){
             for(Index_type tx=0; tx<local_dim; ++tx){
               dUhatdX[ty][tx] = 0.0;
               dUdX[ty][tx] = 0.0;
             }
           }
           
           
           for(Index_type a=0; a<inumNodesPerElement; ++a)
             {
               dUdX[0][0] += u_local_x[a]*idNdX_x(k,q,a);
               dUdX[0][1] += u_local_x[a]*idNdX_y(k,q,a);
               dUdX[0][2] += u_local_x[a]*idNdX_z(k,q,a);
               
               dUdX[1][0] += u_local_y[a]*idNdX_x(k,q,a);
               dUdX[1][1] += u_local_y[a]*idNdX_y(k,q,a);
               dUdX[1][2] += u_local_y[a]*idNdX_z(k,q,a);
               
               dUdX[2][0] += u_local_z[a]*idNdX_x(k,q,a);
               dUdX[2][1] += u_local_z[a]*idNdX_y(k,q,a);
               dUdX[2][2] += u_local_z[a]*idNdX_z(k,q,a);
             }
           
           for(Index_type a=0; a<inumNodesPerElement; ++a)
             {
               dUhatdX[0][0] += uhat_local_x[a]*idNdX_x(k,q,a);
               dUhatdX[0][1] += uhat_local_x[a]*idNdX_y(k,q,a);
               dUhatdX[0][2] += uhat_local_x[a]*idNdX_z(k,q,a);
               
               dUhatdX[1][0] += uhat_local_y[a]*idNdX_x(k,q,a);
               dUhatdX[1][1] += uhat_local_y[a]*idNdX_y(k,q,a);
               dUhatdX[1][2] += uhat_local_y[a]*idNdX_z(k,q,a);
               
               dUhatdX[2][0] += uhat_local_z[a]*idNdX_x(k,q,a);
               dUhatdX[2][1] += uhat_local_z[a]*idNdX_y(k,q,a);
               dUhatdX[2][2] += uhat_local_z[a]*idNdX_z(k,q,a);
             }
           
           
           real64 F[local_dim][local_dim];
           real64 Finv[local_dim][local_dim];
           real64 L[local_dim][local_dim];
           
           {            
             real64 dvdX[local_dim][local_dim];
             
             
             for(Index_type ty=0; ty<local_dim; ++ty){
               for(Index_type tx=0; tx<local_dim; ++tx){
                 dvdX[ty][tx] = dUhatdX[ty][tx]*(1.0/dt);
               }
             }
             
             //calculate du/dX
             for(Index_type row=0; row<local_dim; ++row){
               for(Index_type col=0; col<local_dim; ++col){
                F[row][col] = dUhatdX[row][col];
               }
             }
             
             for(Index_type row=0; row<local_dim; ++row){
               for(Index_type col=0; col<local_dim; ++col){
                 F[row][col] *= 0.5;
               }
             }
             
             for(Index_type row=0; row<local_dim; ++row){
               for(Index_type col=0; col<local_dim; ++col){
                 F[row][col] += dUdX[row][col];
               }
             }
             
             for(Index_type tx=0; tx<local_dim; ++tx){
               F[tx][tx] += 1.0;
             }
             
            
             Finverse(F, Finv);
             
             
             AijBjk(dvdX,Finv,L);
           }
           
           //Calculate gradient (end of step)
           for(Index_type ty=0; ty<local_dim; ++ty){
             for(Index_type tx=0; tx<local_dim; ++tx){
               F[ty][tx] = dUhatdX[ty][tx] + dUdX[ty][tx];
             }
           }

           for(int tx=0; tx<local_dim; ++tx){ F[tx][tx] += 1.0;};
           
           real64 detF = det<real64>(F);
           Finverse<real64> (F, Finv);
           
           real64 Rot[local_dim][local_dim];
           real64 Dadt[local_dim][local_dim];
           
           //-------------[Hughues Winget]--------------
           HughesWinget(Rot, Dadt, L, dt);
           //-------------------------------------------

           //Write out data to global memory
           detF_ptr(k,q) = detF;
           for(Index_type r = 0; r < local_dim; ++r){
             for(Index_type c = 0; c < local_dim; ++c){
               Dadt_ptr(k,q,r,c) = Dadt[r][c];
               Rot_ptr(k,q,r,c) = Rot[r][c];
               Finv_ptr(k,q,r,c) = Finv[r][c];
             }
           }           
           
         }//end of quadrature
       
     });
   
   
}      

///
///Solid mechanics kinematic kernel with nodal degrees of freedom stored and shape function
///derivivatives stored in an array of objects format. This kernel only performs the kinematic step.
///Domain is assumed to be structured.      
///            
template<typename Pol>
RAJA_INLINE void ObjectOfArrays_KinematicKernel_Structured(Index_type NoElem, RAJA::TypedListSegment<Index_type> elemList, real64 dt,
                                                           Index_type nx, Index_type ny, Index_type nz, 
                                                           geosxData iu_x, geosxData iu_y, geosxData iu_z,
                                                           geosxData iuhat_x, geosxData iuhat_y, geosxData iuhat_z,
                                                           geosxData idNdX_x, geosxData idNdX_y, geosxData idNdX_z,
                                                           Index_type const * iconstitutiveMap, geosxData idevStressData,
                                                           geosxData imeanStress, real64 ishearModulus, real64 ibulkModulus,
                                                           const real64 * idetJ,
                                                           geosxData iacc_x, geosxData iacc_y, geosxData iacc_z,
                                                           geosxData Dadt_ptr, geosxData Rot_ptr,
                                                           geosxData detF_ptr, geosxData Finv_ptr)
{
      

#if defined(CUDA)
  RAJA::forall<Pol>
    (elemList,[=] RAJA_DEVICE (Index_type k){
#else
   RAJA::forall<Pol>
     (elemList,[=](Index_type k){
#endif
       
       real64 uhat_local_x[inumNodesPerElement];
       real64 uhat_local_y[inumNodesPerElement];
       real64 uhat_local_z[inumNodesPerElement];
       
       real64 u_local_x[inumNodesPerElement];
       real64 u_local_y[inumNodesPerElement];
       real64 u_local_z[inumNodesPerElement];

      //Gather element to node list
      Index_type nodeList[8];
      structuredElemToNodes(nodeList,k,nx,ny,nz);
       
       //Copy Global to Local
       for(Index_type a=0; a<inumNodesPerElement; ++a)
         {

           Index_type id = nodeList[a];
           u_local_x[a] = iu_x[id];
           u_local_y[a] = iu_y[id];
           u_local_z[a] = iu_z[id];
           
           uhat_local_x[a] = iuhat_x[id];
           uhat_local_y[a] = iuhat_y[id];
           uhat_local_z[a] = iuhat_z[id];
         }
       
       
       //Compute Quadrature
       for(Index_type q=0; q<inumQuadraturePoints; ++q)
         {           
           
           real64 dUhatdX[local_dim][local_dim];
           real64 dUdX[local_dim][local_dim];
           
           for(Index_type ty=0; ty<local_dim; ++ty){
             for(Index_type tx=0; tx<local_dim; ++tx){
               dUhatdX[ty][tx] = 0.0;
               dUdX[ty][tx] = 0.0;
             }
           }
           
           
           for(Index_type a=0; a<inumNodesPerElement; ++a)
             {
               dUdX[0][0] += u_local_x[a]*idNdX_x(k,q,a);
               dUdX[0][1] += u_local_x[a]*idNdX_y(k,q,a);
               dUdX[0][2] += u_local_x[a]*idNdX_z(k,q,a);
               
               dUdX[1][0] += u_local_y[a]*idNdX_x(k,q,a);
               dUdX[1][1] += u_local_y[a]*idNdX_y(k,q,a);
               dUdX[1][2] += u_local_y[a]*idNdX_z(k,q,a);
               
               dUdX[2][0] += u_local_z[a]*idNdX_x(k,q,a);
               dUdX[2][1] += u_local_z[a]*idNdX_y(k,q,a);
               dUdX[2][2] += u_local_z[a]*idNdX_z(k,q,a);
             }
           
           for(Index_type a=0; a<inumNodesPerElement; ++a)
             {
               dUhatdX[0][0] += uhat_local_x[a]*idNdX_x(k,q,a);
               dUhatdX[0][1] += uhat_local_x[a]*idNdX_y(k,q,a);
               dUhatdX[0][2] += uhat_local_x[a]*idNdX_z(k,q,a);
               
               dUhatdX[1][0] += uhat_local_y[a]*idNdX_x(k,q,a);
               dUhatdX[1][1] += uhat_local_y[a]*idNdX_y(k,q,a);
               dUhatdX[1][2] += uhat_local_y[a]*idNdX_z(k,q,a);
               
               dUhatdX[2][0] += uhat_local_z[a]*idNdX_x(k,q,a);
               dUhatdX[2][1] += uhat_local_z[a]*idNdX_y(k,q,a);
               dUhatdX[2][2] += uhat_local_z[a]*idNdX_z(k,q,a);
             }
           
           
           real64 F[local_dim][local_dim];
           real64 Finv[local_dim][local_dim];
           real64 L[local_dim][local_dim];
           
           {            
             real64 dvdX[local_dim][local_dim];
             
             
             for(Index_type ty=0; ty<local_dim; ++ty){
               for(Index_type tx=0; tx<local_dim; ++tx){
                 dvdX[ty][tx] = dUhatdX[ty][tx]*(1.0/dt);
               }
             }
             
             //calculate du/dX
             for(Index_type row=0; row<local_dim; ++row){
               for(Index_type col=0; col<local_dim; ++col){
                F[row][col] = dUhatdX[row][col];
               }
             }
             
             for(Index_type row=0; row<local_dim; ++row){
               for(Index_type col=0; col<local_dim; ++col){
                 F[row][col] *= 0.5;
               }
             }
             
             for(Index_type row=0; row<local_dim; ++row){
               for(Index_type col=0; col<local_dim; ++col){
                 F[row][col] += dUdX[row][col];
               }
             }
             
             for(Index_type tx=0; tx<local_dim; ++tx){
               F[tx][tx] += 1.0;
             }
             
            
             Finverse(F, Finv);
             
             
             AijBjk(dvdX,Finv,L);
           }
           
           //Calculate gradient (end of step)
           for(Index_type ty=0; ty<local_dim; ++ty){
             for(Index_type tx=0; tx<local_dim; ++tx){
               F[ty][tx] = dUhatdX[ty][tx] + dUdX[ty][tx];
             }
           }

           for(int tx=0; tx<local_dim; ++tx){ F[tx][tx] += 1.0;};
           
           real64 detF = det<real64>(F);
           Finverse<real64> (F, Finv);
           
           real64 Rot[local_dim][local_dim];
           real64 Dadt[local_dim][local_dim];
           
           //-------------[Hughues Winget]--------------
           HughesWinget(Rot, Dadt, L, dt);
           //-------------------------------------------

           //Write out data to global memory
           detF_ptr(k,q) = detF;
           for(Index_type r = 0; r < local_dim; ++r){
             for(Index_type c = 0; c < local_dim; ++c){
               Dadt_ptr(k,q,r,c) = Dadt[r][c];
               Rot_ptr(k,q,r,c) = Rot[r][c];
               Finv_ptr(k,q,r,c) = Finv[r][c];
             }
           }           
           
         }//end of quadrature
       
     });
   
   
}      

///
///Solid mechanics integration kernel with nodal degrees of freedom stored and shape function
///derivivatives stored in an objects of arrays format. This kernel only performs the kinematic step.
///
template<typename Pol>
RAJA_INLINE void ObjectOfArrays_IntegrationKernel(Index_type NoElem, RAJA::TypedListSegment<Index_type> elemList, real64 dt,
                                                  const Index_type * elemsToNodes,                          
                                                  geosxData iu_x, geosxData iu_y, geosxData iu_z,
                                                  geosxData iuhat_x, geosxData iuhat_y, geosxData iuhat_z,
                                                  geosxData idNdX_x, geosxData idNdX_y, geosxData idNdX_z,
                                                  Index_type const * iconstitutiveMap, geosxData idevStressData,
                                                  geosxData imeanStress, real64 ishearModulus, real64 ibulkModulus,
                                                  const real64 * idetJ, geosxData iacc_x, geosxData iacc_y, geosxData iacc_z,
                                                  geosxData Dadt_ptr, geosxData Rot_ptr,
                                                  geosxData detF_ptr, geosxData Finv_ptr)
{

  #if defined(CUDA)
  RAJA::forall<Pol>
    (elemList,[=] RAJA_DEVICE (Index_type k){
#else
   RAJA::forall<Pol>
     (elemList,[=](Index_type k){
#endif

       real64 f_local_x[inumNodesPerElement]; for(int i=0; i<inumNodesPerElement; ++i) f_local_x[i] = 0;
       real64 f_local_y[inumNodesPerElement]; for(int i=0; i<inumNodesPerElement; ++i) f_local_y[i] = 0;
       real64 f_local_z[inumNodesPerElement]; for(int i=0; i<inumNodesPerElement; ++i) f_local_z[i] = 0;

       
       for(Index_type q=0; q < inumQuadraturePoints; ++q){

         Index_type m = iconstitutiveMap(k,q);
         real64 TotalStress[local_dim][local_dim];
         
         TotalStress[0][0] = idevStressData(k,q,0);
         TotalStress[1][0] = idevStressData(k,q,1);
         TotalStress[1][1] = idevStressData(k,q,2);
         TotalStress[2][0] = idevStressData(k,q,3);
         TotalStress[2][1] = idevStressData(k,q,4);
         TotalStress[2][2] = idevStressData(k,q,5);
         
         TotalStress[0][1] = TotalStress[1][0];
         TotalStress[0][2] = TotalStress[2][0];
         TotalStress[1][2] = TotalStress[2][1];
         
          for(Index_type i=0; i<local_dim; ++i)
            {
              TotalStress[i][i] += imeanStress[m];
            }

          //---------[Integrate - Function]---------------------
          real64 const integrationFactor = idetJ(k,q)*detF_ptr(k,q);
          real64 P[local_dim][local_dim];
          real64 Finv[local_dim][local_dim];

          for(Index_type r=0; r<local_dim; ++r){
            for(Index_type c=0; c<local_dim; ++c){
              Finv[r][c] = Finv_ptr(k,q,r,c);
            }
          }
          
          AijBkj(TotalStress,Finv,P);
          for(Index_type ty=0; ty<local_dim; ++ty){
            for(Index_type tx=0; tx<local_dim; ++tx){
              P[ty][tx] *= integrationFactor;
            }
          } 
          
          //--------------------------------------------------
          for(int a=0; a<inumNodesPerElement; ++a){
            f_local_x[a] -= P[0][0]*idNdX_x(k,q,a) + P[0][1]*idNdX_y(k,q,a) + P[0][2]*idNdX_z(k,q,a);
            f_local_y[a] -= P[1][0]*idNdX_x(k,q,a) + P[1][1]*idNdX_y(k,q,a) + P[1][2]*idNdX_z(k,q,a);
            f_local_z[a] -= P[2][0]*idNdX_x(k,q,a) + P[2][1]*idNdX_y(k,q,a) + P[2][2]*idNdX_z(k,q,a);
           }

        }//end of quadrature

      //-------------[Add Local To Global]----------------
      for(Index_type a=0; a<inumNodesPerElement; ++a)
        {
          Index_type id = elemsToNodes[a + inumNodesPerElement*k];
          RAJA::atomic::atomicAdd<atomicPol>(&iacc_x[id],f_local_x[a]);
          RAJA::atomic::atomicAdd<atomicPol>(&iacc_y[id],f_local_y[a]);
          RAJA::atomic::atomicAdd<atomicPol>(&iacc_z[id],f_local_z[a]);
        }                

                   
     });
      
}

///
///Solid mechanics integration kernel with nodal degrees of freedom stored and shape function
///derivivatives stored in an objects of arrays format. This kernel only performs the kinematic step.
///Domain is assumed to be structured      
///
template<typename Pol>
RAJA_INLINE void ObjectOfArrays_IntegrationKernel_Structured(Index_type NoElem, RAJA::TypedListSegment<Index_type> elemList, real64 dt,
                                                           Index_type nx, Index_type ny, Index_type nz,
                                                           geosxData iu_x, geosxData iu_y, geosxData iu_z,
                                                           geosxData iuhat_x, geosxData iuhat_y, geosxData iuhat_z,
                                                           geosxData idNdX_x, geosxData idNdX_y, geosxData idNdX_z,
                                                           Index_type const * iconstitutiveMap, geosxData idevStressData,
                                                           geosxData imeanStress, real64 ishearModulus, real64 ibulkModulus,
                                                           const real64 * idetJ, geosxData iacc_x, geosxData iacc_y, geosxData iacc_z,
                                                           geosxData Dadt_ptr, geosxData Rot_ptr,
                                                           geosxData detF_ptr, geosxData Finv_ptr)
{

  #if defined(CUDA)
  RAJA::forall<Pol>
    (elemList,[=] RAJA_DEVICE (Index_type k){
#else
   RAJA::forall<Pol>
     (elemList,[=](Index_type k){
#endif

       real64 f_local_x[inumNodesPerElement]; for(int i=0; i<inumNodesPerElement; ++i) f_local_x[i] = 0;
       real64 f_local_y[inumNodesPerElement]; for(int i=0; i<inumNodesPerElement; ++i) f_local_y[i] = 0;
       real64 f_local_z[inumNodesPerElement]; for(int i=0; i<inumNodesPerElement; ++i) f_local_z[i] = 0;


       //Gather element to node list
       Index_type nodeList[8];
       structuredElemToNodes(nodeList,k,nx,ny,nz);
       
       for(Index_type q=0; q < inumQuadraturePoints; ++q){

         Index_type m = iconstitutiveMap(k,q);
         real64 TotalStress[local_dim][local_dim];
         
         TotalStress[0][0] = idevStressData(k,q,0);
         TotalStress[1][0] = idevStressData(k,q,1);
         TotalStress[1][1] = idevStressData(k,q,2);
         TotalStress[2][0] = idevStressData(k,q,3);
         TotalStress[2][1] = idevStressData(k,q,4);
         TotalStress[2][2] = idevStressData(k,q,5);
         
         TotalStress[0][1] = TotalStress[1][0];
         TotalStress[0][2] = TotalStress[2][0];
         TotalStress[1][2] = TotalStress[2][1];
         
          for(Index_type i=0; i<local_dim; ++i)
            {
              TotalStress[i][i] += imeanStress[m];
            }

          //---------[Integrate - Function]---------------------
          real64 const integrationFactor = idetJ(k,q)*detF_ptr(k,q);
          real64 P[local_dim][local_dim];
          real64 Finv[local_dim][local_dim];

          for(Index_type r=0; r<local_dim; ++r){
            for(Index_type c=0; c<local_dim; ++c){
              Finv[r][c] = Finv_ptr(k,q,r,c);
            }
          }
          
          AijBkj(TotalStress,Finv,P);
          for(Index_type ty=0; ty<local_dim; ++ty){
            for(Index_type tx=0; tx<local_dim; ++tx){
              P[ty][tx] *= integrationFactor;
            }
          } 
          
          //--------------------------------------------------
          for(int a=0; a<inumNodesPerElement; ++a){
            f_local_x[a] -= P[0][0]*idNdX_x(k,q,a) + P[0][1]*idNdX_y(k,q,a) + P[0][2]*idNdX_z(k,q,a);
            f_local_y[a] -= P[1][0]*idNdX_x(k,q,a) + P[1][1]*idNdX_y(k,q,a) + P[1][2]*idNdX_z(k,q,a);
            f_local_z[a] -= P[2][0]*idNdX_x(k,q,a) + P[2][1]*idNdX_y(k,q,a) + P[2][2]*idNdX_z(k,q,a);
           }

        }//end of quadrature

      //-------------[Add Local To Global]----------------
      for(Index_type a=0; a<inumNodesPerElement; ++a)
        {
          Index_type id = nodeList[a];
          RAJA::atomic::atomicAdd<atomicPol>(&iacc_x[id],f_local_x[a]);
          RAJA::atomic::atomicAdd<atomicPol>(&iacc_y[id],f_local_y[a]);
          RAJA::atomic::atomicAdd<atomicPol>(&iacc_z[id],f_local_z[a]);
        }                

                   
     });
      
}
      

      
  
}//namespace
      


#endif
