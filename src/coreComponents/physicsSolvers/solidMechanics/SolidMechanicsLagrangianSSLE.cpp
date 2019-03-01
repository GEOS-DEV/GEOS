/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
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

/**
 * @file SolidMechanicsLagrangianSSLE.hpp
 */

#include "SolidMechanicsLagrangianSSLE.hpp"

#include "codingUtilities/Utilities.hpp"
#include "finiteElement/Kinematics.h"



namespace geosx
{

using namespace constitutive;

SolidMechanicsLagrangianSSLE::SolidMechanicsLagrangianSSLE( string const & name,
                                                            ManagedGroup * const parent ):
  SolidMechanicsLagrangianFEM( name, parent )
{
  this->m_strainTheory = 0;
}

SolidMechanicsLagrangianSSLE::~SolidMechanicsLagrangianSSLE()
{}

template< localIndex NUM_NODES_PER_ELEM, localIndex NUM_QUADRATURE_POINTS >
real64 SolidMechanicsLagrangianSSLE::ExplicitElementKernelWrapper::
Launch( localIndex const er,
        localIndex const esr,
        set<localIndex> const & elementList,
        arrayView2d<localIndex> const & elemsToNodes,
        arrayView3d< R1Tensor > const & dNdX,
        arrayView2d<real64> const & detJ,
        arrayView1d<R1Tensor> const & u,
        arrayView1d<R1Tensor> const & vel,
        arrayView1d<R1Tensor> & acc,
        ElementRegionManager::ConstitutiveRelationAccessor<ConstitutiveBase> constitutiveRelations,
        ElementRegionManager::MaterialViewAccessor< arrayView2d<real64> > const & meanStress,
        ElementRegionManager::MaterialViewAccessor< arrayView2d<R2SymTensor> > const & devStress,
        real64 const dt )
{

  ConstitutiveBase::UpdateFunctionPointer update = constitutiveRelations[er][esr][0]->GetStateUpdateFunctionPointer();
  void * data = nullptr;
  constitutiveRelations[er][esr][0]->SetParamStatePointers( data );
  forall_in_set<elemPolicy>( elementList.values(),
                             elementList.size(),
                             GEOSX_LAMBDA ( localIndex const k)
  {
    R1Tensor v_local[ NUM_NODES_PER_ELEM ];
    R1Tensor u_local[ NUM_NODES_PER_ELEM ];
    R1Tensor f_local[ NUM_NODES_PER_ELEM ];

    real64 c[6][6];
    constitutiveRelations[er][esr][0]->GetStiffness( c );


    CopyGlobalToLocal<NUM_NODES_PER_ELEM, R1Tensor>( elemsToNodes[k],
                                                     u, vel,
                                                     u_local, v_local );

    //Compute Quadrature
    for( localIndex q = 0 ; q<NUM_QUADRATURE_POINTS ; ++q)
    {

      real64 p_stress[6] = {0};
//      real64 p_Cdamp[6] = {0};
      for( localIndex a=0 ; a<NUM_NODES_PER_ELEM ; ++a )
      {
        R1Tensor const & v_a = v_local[a];

        const R1Tensor& dNdXb = dNdX[k][q][a];
        const realT v0_x_dNdXb0 = v_a[0]*dNdXb[0];
        const realT v1_x_dNdXb1 = v_a[1]*dNdXb[1];
        const realT v2_x_dNdXb2 = v_a[2]*dNdXb[2];

        p_stress[0] += ( v0_x_dNdXb0*c[0][0] + v1_x_dNdXb1*c[0][1] + v2_x_dNdXb2*c[0][2] ) * dt;
        p_stress[1] += ( v0_x_dNdXb0*c[1][0] + v1_x_dNdXb1*c[1][1] + v2_x_dNdXb2*c[1][2] ) * dt;
        p_stress[2] += ( v0_x_dNdXb0*c[2][0] + v1_x_dNdXb1*c[2][1] + v2_x_dNdXb2*c[2][2] ) * dt;
        p_stress[3] += ( v_a[2]*dNdXb[1] + v_a[1]*dNdXb[2] )*c[3][3] * dt;
        p_stress[4] += ( v_a[2]*dNdXb[0] + v_a[0]*dNdXb[2] )*c[4][4] * dt;
        p_stress[5] += ( v_a[1]*dNdXb[0] + v_a[0]*dNdXb[1] )*c[5][5] * dt;

//        p_Cdamp[0] += v0_x_dNdXb0*c[0][0] + v1_x_dNdXb1*c[0][1] + v2_x_dNdXb2*c[0][2];
//        p_Cdamp[1] += v0_x_dNdXb0*c[1][0] + v1_x_dNdXb1*c[1][1] + v2_x_dNdXb2*c[1][2];
//        p_Cdamp[2] += v0_x_dNdXb0*c[2][0] + v1_x_dNdXb1*c[2][1] + v2_x_dNdXb2*c[2][2];
//        p_Cdamp[3] += ( v_a[2]*dNdXb[1] + v_a[1]*dNdXb[2] )*c[3][3];
//        p_Cdamp[4] += ( v_a[2]*dNdXb[0] + v_a[0]*dNdXb[2] )*c[4][4];
//        p_Cdamp[5] += ( v_a[1]*dNdXb[0] + v_a[0]*dNdXb[1] )*c[5][5];
      }
      real64 const dMeanStress = ( p_stress[0] + p_stress[1] + p_stress[2] )/3.0;
      meanStress[er][esr][0][k][q] += dMeanStress;

      p_stress[0] -= dMeanStress;
      p_stress[1] -= dMeanStress;
      p_stress[2] -= dMeanStress;

      real64 * const restrict p_devStress = devStress[er][esr][0][k][q].Data();
      p_devStress[0] += p_stress[0];
      p_devStress[2] += p_stress[1];
      p_devStress[5] += p_stress[2];
      p_devStress[4] += p_stress[3];
      p_devStress[3] += p_stress[4];
      p_devStress[1] += p_stress[5];

      for( localIndex a=0 ; a<NUM_NODES_PER_ELEM ; ++a )
      {
        const R1Tensor& dNdXa = dNdX[k][q][a];

        f_local[a][0] -= ( p_devStress[1]*dNdXa[1]
                      + p_devStress[3]*dNdXa[2]
                      + dNdXa[0]*(p_devStress[0] + meanStress[er][esr][0][k][q]) ) * detJ[k][q];
        f_local[a][1] -= ( p_devStress[1]*dNdXa[0]
                      + p_devStress[4]*dNdXa[2]
                      + dNdXa[1]*(p_devStress[2] + meanStress[er][esr][0][k][q]) ) * detJ[k][q];
        f_local[a][2] -= ( p_devStress[3]*dNdXa[0]
                      + p_devStress[4]*dNdXa[1]
                      + dNdXa[2]*(p_devStress[5] + meanStress[er][esr][0][k][q]) ) * detJ[k][q];
      }
    }//quadrature loop


    AddLocalToGlobal<NUM_NODES_PER_ELEM>( elemsToNodes[k], f_local, acc );
  });

  return dt;
}
template
real64 SolidMechanicsLagrangianSSLE::
ExplicitElementKernelWrapper::Launch<8l,8l>( localIndex const er,
                                  localIndex const esr,
                                  set<localIndex> const & elementList,
                                  arrayView2d<localIndex> const & elemsToNodes,
                                  arrayView3d< R1Tensor > const & dNdX,
                                  arrayView2d<real64> const & detJ,
                                  arrayView1d<R1Tensor> const & u,
                                  arrayView1d<R1Tensor> const & vel,
                                  arrayView1d<R1Tensor> & acc,
                                  ElementRegionManager::ConstitutiveRelationAccessor<ConstitutiveBase> constitutiveRelations,
                                  ElementRegionManager::MaterialViewAccessor< arrayView2d<real64> > const & meanStress,
                                  ElementRegionManager::MaterialViewAccessor< arrayView2d<R2SymTensor> > const & devStress,
                                  real64 const dt );

REGISTER_CATALOG_ENTRY( SolverBase, SolidMechanicsLagrangianSSLE, string const &, dataRepository::ManagedGroup * const )
} /* namespace geosx */

