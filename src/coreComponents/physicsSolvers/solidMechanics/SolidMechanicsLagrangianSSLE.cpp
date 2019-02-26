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

#include "SolidMechanicsLagrangianSSLE.hpp"

#include "codingUtilities/Utilities.hpp"
#include "finiteElement/Kinematics.h"



namespace geosx
{

using namespace constitutive;

SolidMechanics_LagrangianSSLE::SolidMechanics_LagrangianSSLE( string const & name,
                                                              ManagedGroup * const parent ):
  SolidMechanicsLagrangianFEM( name, parent )
{
}

SolidMechanics_LagrangianSSLE::~SolidMechanics_LagrangianSSLE()
{}

template< localIndex NUM_NODES_PER_ELEM, localIndex NUM_QUADRATURE_POINTS >
real64 SolidMechanics_LagrangianSSLE::ExplicitElementKernelLaunch( localIndex const er,
                                                            localIndex const esr,
                                                            set<localIndex> const & elementList,
                                                            arrayView2d<localIndex> const & elemsToNodes,
                                                            arrayView3d< R1Tensor > const & dNdX,
                                                            arrayView2d<real64> const & detJ,
                                                            arrayView1d<R1Tensor> const & u,
                                                            arrayView1d<R1Tensor> const & uhat,
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
    R1Tensor uhat_local[ NUM_NODES_PER_ELEM ];
    R1Tensor v_local[ NUM_NODES_PER_ELEM ];
    R1Tensor u_local[ NUM_NODES_PER_ELEM ];
    R1Tensor f_local[ NUM_NODES_PER_ELEM ];

    real64 c[6][6];
    constitutiveRelations[er][esr][0]->GetStiffness( c );

    CopyGlobalToLocal<NUM_NODES_PER_ELEM, R1Tensor>( elemsToNodes[k],
                                           u, uhat, vel,
                                           u_local, uhat_local, v_local );

    //Compute Quadrature
    for( localIndex q = 0 ; q<NUM_QUADRATURE_POINTS ; ++q)
    {

      for( localIndex a=0 ; a<NUM_NODES_PER_ELEM ; ++a )
      {
        R1Tensor const & uhat_a = uhat_local[a];
        R1Tensor const & v_a = v_local[a];

        const R1Tensor& dNdXb = dNdX[k][q][a];
        const realT u0_x_dNdXb0 = uhat_a[0]*dNdXb[0];
        const realT u1_x_dNdXb1 = uhat_a[1]*dNdXb[1];

        const realT v0_x_dNdXb0 = v_a[0]*dNdXb[0];
        const realT v1_x_dNdXb1 = v_a[1]*dNdXb[1];

        const realT u2_x_dNdXb2 = uhat_a[2]*dNdXb[2];
        const realT v2_x_dNdXb2 = v_a[2]*dNdXb[2];

        real64 p_stress[6], p_Cdamp[6];
        p_stress[0] += u0_x_dNdXb0*c[0][0] + u1_x_dNdXb1*c[0][1] + u2_x_dNdXb2*c[0][2];
        p_stress[2] += u0_x_dNdXb0*c[1][0] + u1_x_dNdXb1*c[1][1] + u2_x_dNdXb2*c[1][2];
        p_stress[5] += u0_x_dNdXb0*c[2][0] + u1_x_dNdXb1*c[2][1] + u2_x_dNdXb2*c[2][2];
        p_stress[4] += ( uhat_a[2]*dNdXb[1] + uhat_a[1]*dNdXb[2] )*c[3][3];
        p_stress[3] += ( uhat_a[2]*dNdXb[0] + uhat_a[0]*dNdXb[2] )*c[4][4];
        p_stress[1] += ( uhat_a[1]*dNdXb[0] + uhat_a[0]*dNdXb[1] )*c[5][5];

        p_Cdamp[0] += v0_x_dNdXb0*c[0][0] + v1_x_dNdXb1*c[0][1] + v2_x_dNdXb2*c[0][2];
        p_Cdamp[2] += v0_x_dNdXb0*c[1][0] + v1_x_dNdXb1*c[1][1] + v2_x_dNdXb2*c[1][2];
        p_Cdamp[5] += v0_x_dNdXb0*c[2][0] + v1_x_dNdXb1*c[2][1] + v2_x_dNdXb2*c[2][2];
        p_Cdamp[4] += ( v_a[2]*dNdXb[1] + v_a[1]*dNdXb[2] )*c[3][3];
        p_Cdamp[3] += ( v_a[2]*dNdXb[0] + v_a[0]*dNdXb[2] )*c[4][4];
        p_Cdamp[1] += ( v_a[1]*dNdXb[0] + v_a[0]*dNdXb[1] )*c[5][5];
      }

    }//quadrature loop


    AddLocalToGlobal<NUM_NODES_PER_ELEM>( elemsToNodes[k], f_local, acc );
  });

  return dt;
}
template
real64 SolidMechanics_LagrangianSSLE::
ExplicitElementKernelLaunch<8l,8l>( localIndex const er,
                                  localIndex const esr,
                                  set<localIndex> const & elementList,
                                  arrayView2d<localIndex> const & elemsToNodes,
                                  arrayView3d< R1Tensor > const & dNdX,
                                  arrayView2d<real64> const & detJ,
                                  arrayView1d<R1Tensor> const & u,
                                  arrayView1d<R1Tensor> const & uhat,
                                  arrayView1d<R1Tensor> const & vel,
                                  arrayView1d<R1Tensor> & acc,
                                  ElementRegionManager::ConstitutiveRelationAccessor<ConstitutiveBase> constitutiveRelations,
                                  ElementRegionManager::MaterialViewAccessor< arrayView2d<real64> > const & meanStress,
                                  ElementRegionManager::MaterialViewAccessor< arrayView2d<R2SymTensor> > const & devStress,
                                  real64 const dt );

} /* namespace geosx */

