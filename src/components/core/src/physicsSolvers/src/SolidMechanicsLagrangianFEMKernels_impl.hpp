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
#include "../miniApps/SolidMechanicsLagrangianFEM-MiniApp/MatrixMath_impl.hpp"
#include "../miniApps/SolidMechanicsLagrangianFEM-MiniApp/Layout.hpp"
#include "../miniApps/SolidMechanicsLagrangianFEM-MiniApp/ShapeFun_impl.hpp"
#include "../miniApps/SolidMechanicsLagrangianFEM-MiniApp/ConstitutiveUpdate_impl.hpp"
#include "../../core/src/rajaInterface/GEOS_RAJA_Interface.hpp"

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
void ObjectOfArraysKernel(Index_type noElem, RAJA::TypedListSegment<Index_type> elemList, real64 dt,
                          const Index_type * elemsToNodes, geosxData iu_x,
                          geosxData iu_y, geosxData iu_z,
                          geosxData iuhat_x,
                          geosxData iuhat_y, geosxData iuhat_z,
                          geosxData idNdX_x,
                          geosxData idNdX_y,geosxData idNdX_z,
                          Index_type const * iconstitutiveMap, geosxData idevStressData,
                          geosxData imeanStress, real64 ishearModulus, real64 ibulkModulus,
                          const real64 * idetJ, geosxData iacc_x,
                          geosxData iacc_y, geosxData iacc_z, constUpdate updateState_ptr, 
                          Index_type nx=2, Index_type ny=2, Index_type nz=2)
{


  forall_in_set<Pol>(elemList, GEOSX_LAMBDA (Index_type k) {
      
      real64 uhat_local_x[inumNodesPerElement];
      real64 uhat_local_y[inumNodesPerElement];
      real64 uhat_local_z[inumNodesPerElement];

      real64 u_local_x[inumNodesPerElement];
      real64 u_local_y[inumNodesPerElement];
      real64 u_local_z[inumNodesPerElement];
      
      real64 f_local_x[inumNodesPerElement] = {0};
      real64 f_local_y[inumNodesPerElement] = {0};
      real64 f_local_z[inumNodesPerElement] = {0};


#if defined(STRUCTURED_GRID)
       Index_type nodeList[inumNodesPerElement];       
       structuredElemToNodes(nodeList,k,nx,ny,nz);
#else
       const Index_type *nodeList = (&elemsToNodes[inumNodesPerElement*k]);
#endif       
      
       //Copy Global To Local
       GlobalToLocal(nodeList, k, 
                     u_local_x, u_local_y, u_local_z,
                     uhat_local_x, uhat_local_y, uhat_local_z,
                     iu_x, iu_y, iu_z, iuhat_x, iuhat_y, iuhat_z);
              
      //Compute Quadrature
      for(Index_type q=0; q<inumQuadraturePoints; ++q)
        {


          real64 dUhatdX[local_dim][local_dim] = {{0.0}};
          real64 dUdX[local_dim][local_dim] = {{0.0}};

          //Calculate gradient
          CalculateGradient(dUdX, u_local_x, u_local_y, u_local_z,
                            idNdX_x, idNdX_y, idNdX_z, k, q, noElem);
          
          CalculateGradient(dUhatdX, uhat_local_x, uhat_local_y, uhat_local_z,
                            idNdX_x, idNdX_y, idNdX_z, k, q, noElem);                            

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
          //UpdateStatePoint(Dadt,Rot,m, q,k, idevStressData, imeanStress,ishearModulus, ibulkModulus, noElem);
          updateState_ptr(Dadt,Rot,m, q,k, idevStressData, imeanStress, ishearModulus, ibulkModulus, noElem);
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

          Integrate(f_local_x, f_local_y, f_local_z,
                    idetJ(k,q), detF, Finv, TotalStress, idNdX_x, idNdX_y, idNdX_z, k, q, noElem);
        }//end of quadrature

      //Atomic policy
      AddLocalToGlobal<atomicPol>(nodeList,f_local_x, f_local_y, f_local_z, iacc_x, iacc_y, iacc_z);

    });

}

///
///Solid mechanics update kernel with nodal degrees of freedom and
///shape function derivatives stored in object of arrays format. 
///Computations are carried out in a monothilic kernel.
///Computes shape function derivatives on the fly.
///
template<typename Pol>
RAJA_INLINE
void ObjectOfArraysKernel_Shape(Index_type noElem, RAJA::TypedListSegment<Index_type> elemList, real64 dt,
                                const Index_type * elemsToNodes, geosxData iu_x,
                                geosxData iu_y, geosxData iu_z,
                                geosxData iuhat_x,
                                geosxData iuhat_y, geosxData iuhat_z, const real64 * X, P_Wrapper P,
                                Index_type const * iconstitutiveMap, geosxData idevStressData,
                                geosxData imeanStress, real64 ishearModulus, real64 ibulkModulus,
                                const real64 * idetJ, geosxData iacc_x,
                                geosxData iacc_y, geosxData iacc_z, constUpdate updateState_ptr, 
                                Index_type nx=2, Index_type ny=2, Index_type nz=2)
{

  forall_in_set<Pol>(elemList, GEOSX_LAMBDA (Index_type k) {

      real64 uhat_local_x[inumNodesPerElement];
      real64 uhat_local_y[inumNodesPerElement];
      real64 uhat_local_z[inumNodesPerElement];

      real64 u_local_x[inumNodesPerElement];
      real64 u_local_y[inumNodesPerElement];
      real64 u_local_z[inumNodesPerElement];
      
      real64 f_local_x[inumNodesPerElement] = {0};
      real64 f_local_y[inumNodesPerElement] = {0};
      real64 f_local_z[inumNodesPerElement] = {0};

      //Copy element to node list
#if defined(STRUCTURED_GRID)
       Index_type nodeList[inumNodesPerElement];       
       structuredElemToNodes(nodeList,k,nx,ny,nz);
#else
       const Index_type *nodeList = (&elemsToNodes[inumNodesPerElement*k]);
#endif       
       //Copy Global To Local
       GlobalToLocal(nodeList, k, 
                     u_local_x, u_local_y, u_local_z,
                     uhat_local_x, uhat_local_y, uhat_local_z,
                     iu_x, iu_y, iu_z, iuhat_x, iuhat_y, iuhat_z);

       //Compute shape function derivatives
       real64 dNdX_x[inumQuadraturePoints][inumNodesPerElement];
       real64 dNdX_y[inumQuadraturePoints][inumNodesPerElement];
       real64 dNdX_z[inumQuadraturePoints][inumNodesPerElement];

       //Evaluate shape function derivatives at quadrature points
#if defined(PRE_COMPUTE_P)
       make_dNdX(nodeList, X, dNdX_x, dNdX_y, dNdX_z, P, 
                 inumQuadraturePoints, inumNodesPerElement);
#else
       make_dNdX(nodeList, X, dNdX_x, dNdX_y, dNdX_z,
                 inumQuadraturePoints, inumNodesPerElement);
#endif
       
      //Compute Quadrature
      for(Index_type q=0; q<inumQuadraturePoints; ++q)
        {


          real64 dUhatdX[local_dim][local_dim] = {{0.0}};
          real64 dUdX[local_dim][local_dim] = {{0.0}};

          //Calculate gradient
          CalculateGradient(dUdX, u_local_x, u_local_y, u_local_z,
                            dNdX_x, dNdX_y, dNdX_z, k, q);
          
          CalculateGradient(dUhatdX, uhat_local_x, uhat_local_y, uhat_local_z,
                            dNdX_x, dNdX_y, dNdX_z, k, q);                            

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
          //UpdateStatePoint(Dadt,Rot,m, q,k, idevStressData, imeanStress,ishearModulus, ibulkModulus, noElem);
          updateState_ptr(Dadt,Rot,m, q,k, idevStressData, imeanStress, ishearModulus, ibulkModulus, noElem);
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

          Integrate(f_local_x, f_local_y, f_local_z,
                    idetJ(k,q), detF, Finv, TotalStress, dNdX_x, dNdX_y, dNdX_z, q);
        }//end of quadrature

      //Atomic policy
      AddLocalToGlobal<atomicPol>(nodeList,f_local_x, f_local_y, f_local_z, iacc_x, iacc_y, iacc_z);

    });

} 

///
///Solid mechanics update kernel with nodal degrees of freedom and
///shape function derivatives stored in object of arrays format. 
///Computations are carried out in a monothilic kernel.
// Computes shape function derivatives on the fly.
///
template<typename Pol>
RAJA_INLINE void ArrayOfObjectsKernel_Shape(Index_type noElem, RAJA::TypedListSegment<Index_type> elemList, real64 dt,
                                            const Index_type * elemsToNodes, geosxData iu,
                                            geosxData iuhat, const real64 * X, P_Wrapper P,
                                            Index_type const * iconstitutiveMap, geosxData idevStressData,
                                            geosxData imeanStress, real64 ishearModulus, real64 ibulkModulus,
                                            const real64 * idetJ, geosxData iacc, constUpdate updateState_ptr,
                                            Index_type nx=2, Index_type ny=2, Index_type nz=2)
{


  forall_in_set<Pol>(elemList, GEOSX_LAMBDA (Index_type k) {
      
       real64 uhat_local[local_dim*inumNodesPerElement];
       real64 u_local[local_dim*inumNodesPerElement];
       real64 f_local[local_dim*inumNodesPerElement] = {0};

       //Copy element to node list       
#if defined(STRUCTURED_GRID)
       Index_type nodeList[inumNodesPerElement];       
       structuredElemToNodes(nodeList,k,nx,ny,nz);
#else
       const Index_type *nodeList = (&elemsToNodes[inumNodesPerElement*k]);
#endif       
       
       //Copy Global to Local
       GlobalToLocal(nodeList, k, 
                     u_local,  uhat_local, iu, iuhat);

       //Compute shape function derivatives
       real64 dNdX[inumQuadraturePoints][inumNodesPerElement][local_dim];

       //Compute shape function derivatives at quadrature points
#if defined(PRE_COMPUTE_P)       
       make_dNdX(nodeList, X, dNdX, P, inumQuadraturePoints, inumNodesPerElement);
#else       
       make_dNdX(nodeList, X, dNdX, inumQuadraturePoints, inumNodesPerElement);
#endif       
                                          
      //Compute Quadrature
      for(Index_type q=0; q<inumQuadraturePoints; ++q)
        {
          real64 dUdX[local_dim][local_dim] = {{0.0}};          
          real64 dUhatdX[local_dim][local_dim] = {{0.0}};

          //Calculate Gradient
          CalculateGradient(dUdX, u_local, dNdX, q);
          CalculateGradient(dUhatdX, uhat_local, dNdX, q);
          
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
          //UpdateStatePoint(Dadt,Rot,m, q,k, idevStressData,
          //imeanStress,ishearModulus, ibulkModulus,noElem);
          
          updateState_ptr(Dadt,Rot,m, q,k, idevStressData,
                          imeanStress, ishearModulus, ibulkModulus, noElem);
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
          
          Integrate(f_local, idetJ(k,q), detF, Finv, TotalStress, dNdX, q, noElem); 
        }//end of quadrature

      AddLocalToGlobal<atomicPol>(nodeList, f_local, iacc);      
      
     });

}
      
///
///Solid mechanics update kernel with nodal degrees of freedom and
///shape function derivatives stored as an object of structs format.
///All computations are done in a monolithic kernel.
///      
template<typename Pol>
RAJA_INLINE void ArrayOfObjectsKernel(Index_type noElem, RAJA::TypedListSegment<Index_type> elemList, real64 dt,
                                      const Index_type * elemsToNodes, geosxData iu,
                                      geosxData iuhat, geosxData idNdX,
                                      Index_type const * iconstitutiveMap, geosxData idevStressData,
                                      geosxData imeanStress, real64 ishearModulus, real64 ibulkModulus,
                                      const real64 * idetJ, geosxData iacc, constUpdate updateState_ptr,
                                      Index_type nx=2, Index_type ny=2, Index_type nz=2)
{
  
  forall_in_set<Pol>(elemList, GEOSX_LAMBDA (Index_type k) {
      
       real64 uhat_local[local_dim*inumNodesPerElement];
       real64 u_local[local_dim*inumNodesPerElement];
       real64 f_local[local_dim*inumNodesPerElement] = {0};

#if defined(STRUCTURED_GRID)
       Index_type nodeList[inumNodesPerElement];       
       structuredElemToNodes(nodeList,k,nx,ny,nz);
#else
       const Index_type *nodeList = (&elemsToNodes[inumNodesPerElement*k]);
#endif       

       //Copy Global to Local
       GlobalToLocal(nodeList, k, 
                     u_local,  uhat_local, iu, iuhat);              
                                          
      //Compute Quadrature
      for(Index_type q=0; q<inumQuadraturePoints; ++q)

        {

          real64 dUdX[local_dim][local_dim] = {{0.0}};
          
          real64 dUhatdX[local_dim][local_dim] = {{0.0}};

          //Calculate Gradient
          CalculateGradient(dUdX, u_local, idNdX, k, q, noElem);
          CalculateGradient(dUhatdX, uhat_local, idNdX, k, q, noElem);
          
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
          //UpdateStatePoint(Dadt,Rot,m, q,k, idevStressData,
          //imeanStress,ishearModulus, ibulkModulus,noElem);
          
          updateState_ptr(Dadt,Rot,m, q,k, idevStressData,
                          imeanStress, ishearModulus, ibulkModulus, noElem);
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
                    
          Integrate(f_local, idetJ(k,q), detF, Finv, TotalStress, idNdX, k, q, noElem); 
        }//end of quadrature

      AddLocalToGlobal<atomicPol>(nodeList, f_local, iacc);      
      
     });
   
}


///
///Solid mechanics update kernel with nodal degrees of freedom and
///shape function derivatives stored as an object of structs format.
///All computations are done in a monolithic kernel. Only the kinematic step
///is taken here. 
///            
template<typename Pol>
RAJA_INLINE void ArrayOfObjects_KinematicKernel(Index_type noElem, RAJA::TypedListSegment<Index_type> elemList, real64 dt,
                                                const Index_type * elemsToNodes,
                                                geosxData iu, geosxData iuhat, geosxData idNdX,                          
                                                Index_type const * iconstitutiveMap, geosxData idevStressData,
                                                geosxData imeanStress, real64 ishearModulus, real64 ibulkModulus,
                                                const real64 * idetJ, geosxData iacc, geosxData Dadt_ptr, geosxData Rot_ptr,
                                                geosxData detF_ptr, geosxData Finv_ptr,
                                                Index_type nx=2, Index_type ny=2, Index_type nz=3)
{

  forall_in_set<Pol>(elemList, GEOSX_LAMBDA (Index_type k) {
       
       real64 uhat_local[local_dim*inumNodesPerElement];
       real64 u_local[local_dim*inumNodesPerElement];

#if defined(STRUCTURED_GRID)
       Index_type nodeList[inumNodesPerElement];       
       structuredElemToNodes(nodeList,k,nx,ny,nz);
#else
       const Index_type *nodeList = (&elemsToNodes[inumNodesPerElement*k]);
#endif       
       
       //Copy Global to Local
       GlobalToLocal(nodeList, k, 
                     u_local,  uhat_local, iu, iuhat);        
      
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

          //Calculate Gradient
          CalculateGradient(dUdX, u_local, idNdX, k, q, noElem);

          CalculateGradient(dUhatdX, uhat_local, idNdX, k, q, noElem);
          
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
RAJA_INLINE void ConstitutiveUpdateKernel(Index_type noElem, RAJA::TypedListSegment<Index_type> elemList,
                                          geosxData Dadt_ptr, geosxData Rot_ptr, Index_type const * iconstitutiveMap,
                                          geosxData idevStressData, geosxData imeanStress, real64 shearModulus, real64 bulkModulus)
                          
{


  forall_in_set<Pol>(elemList, GEOSX_LAMBDA (Index_type k) {
      
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
RAJA_INLINE void ArrayOfObjects_IntegrationKernel(Index_type noElem, RAJA::TypedListSegment<Index_type> elemList, real64 dt,
                                                  const Index_type * elemsToNodes,
                                                  geosxData iu, geosxData iuhat, geosxData idNdX,
                                                  Index_type const * iconstitutiveMap, geosxData idevStressData,
                                                  geosxData imeanStress, real64 ishearModulus, real64 ibulkModulus,
                                                  const real64 * idetJ, geosxData iacc, geosxData Dadt_ptr, geosxData Rot_ptr,
                                                  geosxData detF_ptr, geosxData Finv_ptr,
                                                  Index_type nx=2, Index_type ny=2, Index_type nz=2)
{
  
  forall_in_set<Pol>(elemList, GEOSX_LAMBDA (Index_type k) {
       
#if defined(STRUCTURED_GRID)
       Index_type nodeList[inumNodesPerElement];       
       structuredElemToNodes(nodeList,k,nx,ny,nz);
#else
       const Index_type *nodeList = (&elemsToNodes[inumNodesPerElement*k]);
#endif       
       
       real64 f_local[local_dim*inumNodesPerElement] = {0};
       
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
          Integrate(f_local, idetJ(k,q), detF_ptr(k,q), Finv_ptr, TotalStress, idNdX, k, q, noElem);
          //------------------------------------------------------
                              
         }//end of quadrature

       AddLocalToGlobal<atomicPol>(nodeList, f_local, iacc);       
     });
      
}


///
///Solid mechanics kinematic kernel with nodal degrees of freedom stored and shape function
///derivivatives stored in an array of objects format. This kernel only performs the kinematic step.
///            
template<typename Pol>
RAJA_INLINE void ObjectOfArrays_KinematicKernel(Index_type noElem, RAJA::TypedListSegment<Index_type> elemList, real64 dt,
                                                const Index_type * elemsToNodes,
                                                geosxData iu_x, geosxData iu_y, geosxData iu_z,
                                                geosxData iuhat_x, geosxData iuhat_y, geosxData iuhat_z,
                                                geosxData idNdX_x, geosxData idNdX_y, geosxData idNdX_z,
                                                Index_type const * iconstitutiveMap, geosxData idevStressData,
                                                geosxData imeanStress, real64 ishearModulus, real64 ibulkModulus,
                                                const real64 * idetJ,
                                                geosxData iacc_x, geosxData iacc_y, geosxData iacc_z,
                                                geosxData Dadt_ptr, geosxData Rot_ptr,
                                                geosxData detF_ptr, geosxData Finv_ptr,
                                                Index_type nx = 2, Index_type ny = 2, Index_type nz=2)
{


  forall_in_set<Pol>(elemList, GEOSX_LAMBDA (Index_type k) {
       
       real64 uhat_local_x[inumNodesPerElement];
       real64 uhat_local_y[inumNodesPerElement];
       real64 uhat_local_z[inumNodesPerElement];
       
       real64 u_local_x[inumNodesPerElement];
       real64 u_local_y[inumNodesPerElement];
       real64 u_local_z[inumNodesPerElement];

#if defined(STRUCTURED_GRID)
       Index_type nodeList[inumNodesPerElement];       
       structuredElemToNodes(nodeList,k,nx,ny,nz);
#else
       const Index_type *nodeList = (&elemsToNodes[inumNodesPerElement*k]);
#endif       

       //Copy array
       GlobalToLocal(nodeList, k, 
                     u_local_x, u_local_y, u_local_z,
                     uhat_local_x, uhat_local_y, uhat_local_z,
                     iu_x, iu_y, iu_z, iuhat_x, iuhat_y, iuhat_z);
       
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
                      
          //Calculate gradient
          CalculateGradient(dUdX, u_local_x, u_local_y, u_local_z,
                            idNdX_x, idNdX_y, idNdX_z, k, q, noElem);
          
          CalculateGradient(dUhatdX, uhat_local_x, uhat_local_y, uhat_local_z,
                            idNdX_x, idNdX_y, idNdX_z, k, q, noElem);           
           
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

           //Write out intermediate data
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
RAJA_INLINE void ObjectOfArrays_IntegrationKernel(Index_type noElem, RAJA::TypedListSegment<Index_type> elemList, real64 dt,
                                                  const Index_type * elemsToNodes,                          
                                                  geosxData iu_x, geosxData iu_y, geosxData iu_z,
                                                  geosxData iuhat_x, geosxData iuhat_y, geosxData iuhat_z,
                                                  geosxData idNdX_x, geosxData idNdX_y, geosxData idNdX_z,
                                                  Index_type const * iconstitutiveMap, geosxData idevStressData,
                                                  geosxData imeanStress, real64 ishearModulus, real64 ibulkModulus,
                                                  const real64 * idetJ, geosxData iacc_x, geosxData iacc_y, geosxData iacc_z,
                                                  geosxData Dadt_ptr, geosxData Rot_ptr,geosxData detF_ptr, geosxData Finv_ptr,
                                                  Index_type nx=2, Index_type ny=2, Index_type nz=2)
                                                  
{

  forall_in_set<Pol>(elemList, GEOSX_LAMBDA (Index_type k) {

#if defined(STRUCTURED_GRID)
       Index_type nodeList[inumNodesPerElement];       
       structuredElemToNodes(nodeList,k,nx,ny,nz);
#else
       const Index_type *nodeList = (&elemsToNodes[inumNodesPerElement*k]);
#endif       
       
       real64 f_local_x[inumNodesPerElement] = {0};
       real64 f_local_y[inumNodesPerElement] = {0};
       real64 f_local_z[inumNodesPerElement] = {0};
       
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
          Integrate(f_local_x, f_local_y, f_local_z,
                    idetJ(k,q), detF_ptr(k,q), Finv_ptr, TotalStress, idNdX_x, idNdX_y, idNdX_z, k, q, noElem);
        }//end of quadrature

       AddLocalToGlobal<atomicPol>(nodeList, f_local_x, f_local_y, f_local_z, iacc_x, iacc_y, iacc_z);
                   
     });
      
}


//
// Time-stepping routines


//3. v^{n+1/2} = v^{n} + a^{n} dt/2

//One array present
template<typename T, typename U>
void OnePoint(T const dydx, U & y,
              real64 const dx,
              localIndex const length){
  
  forall_in_range(0, length, GEOSX_LAMBDA (localIndex a) {
    y[a][0] += dx*dydx[a][0];
    y[a][1] += dx*dydx[a][1];
    y[a][2] += dx*dydx[a][1];
  });
  
}


//Three arrays present
template<typename T, typename U>
void OnePoint( T const dydx_0,
               T const dydx_1,
               T const dydx_2,
               U &  y,
               real64 const dx,
               localIndex const length )
{  
  forall_in_range(0, length, GEOSX_LAMBDA (localIndex a) {
    y[a][0] += dx*dydx_0[a];
    y[a][1] += dx*dydx_1[a];
    y[a][2] += dx*dydx_2[a];
  });
  
}
  
//                     dydx, dy,   y, dx, length
//4. x^{n+1} = x^{n} + v^{n+{1}/{2}} dt (x is displacement)
template<typename T, typename U>
void OnePoint(U dydx, T dy, T y,
              real64 const dx, localIndex const length){

  forall_in_range(0, length, GEOSX_LAMBDA (localIndex a) {
      
      dy[a][0] = dydx[a][0] * dx;
      dy[a][1] = dydx[a][1] * dx;
      dy[a][2] = dydx[a][2] * dx;

      y[a][0] += dy[a][0];
      y[a][1] += dy[a][1];
      y[a][2] += dy[a][2];      
    });  
}

template<typename T, typename U>
void OnePoint(U dydx,
              T dy_1, T dy_2, T dy_3,
              T y_1, T y_2, T y_3,
              real64 const dx, localIndex const length){

  forall_in_range(0, length, GEOSX_LAMBDA (localIndex a) {
      
    dy_1[a] = dydx[a][0] * dx;
    dy_2[a] = dydx[a][1] * dx;
    dy_3[a] = dydx[a][2] * dx;
    
    y_1[a] += dy_1[a];
    y_2[a] += dy_2[a];
    y_3[a] += dy_3[a];
  });
  
}



      
}//namespace
      


#endif
