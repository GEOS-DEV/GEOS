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

#ifndef CORECOMPONENTS_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSLAGRANGIANSSLE_HPP_
#define CORECOMPONENTS_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSLAGRANGIANSSLE_HPP_

#include "SolidMechanicsLagrangianFEM.hpp"

namespace geosx
{

class SolidMechanicsLagrangianSSLE : public SolidMechanicsLagrangianFEM
{
public:
  SolidMechanicsLagrangianSSLE( string const & name,
                                ManagedGroup * const parent );
  virtual ~SolidMechanicsLagrangianSSLE() override;

  static string CatalogName() { return "SolidMechanicsLagrangianSSLE"; }

  EXPLICIT_ELEMENT_KERNEL_LAUNCH(SolidMechanicsLagrangianSSLE::ExplicitElementKernelWrapper,override)

  IMPLICIT_ELEMENT_KERNEL_LAUNCH(SolidMechanicsLagrangianSSLE::ImplicitElementKernelWrapper,override)

  struct ExplicitElementKernelWrapper
  {
    template< localIndex NUM_NODES_PER_ELEM, localIndex NUM_QUADRATURE_POINTS >
    static real64
    Launch( localIndex const er,
            localIndex const esr,
            set<localIndex> const & elementList,
            arrayView2d<localIndex> const & elemsToNodes,
            arrayView3d< R1Tensor > const & dNdX,
            arrayView2d<real64> const & detJ,
            arrayView1d<R1Tensor> const & u,
            arrayView1d<R1Tensor> const & vel,
            arrayView1d<R1Tensor> & acc,
            ElementRegionManager::ConstitutiveRelationAccessor<constitutive::ConstitutiveBase> constitutiveRelations,
            ElementRegionManager::MaterialViewAccessor< arrayView2d<real64> > const & meanStress,
            ElementRegionManager::MaterialViewAccessor< arrayView2d<R2SymTensor> > const & devStress,
            real64 const dt );
  };

  struct ImplicitElementKernelWrapper
  {
    template< localIndex NUM_NODES_PER_ELEM, localIndex NUM_QUADRATURE_POINTS >
    static real64
    Launch( localIndex const numElems,
            real64 const dt,
            arrayView3d<R1Tensor const> const & dNdX,
            arrayView2d<real64 const > const& detJ,
            FiniteElementBase const * const fe,
            arrayView1d<constitutive::ConstitutiveBase *> const & constitutiveRelations,
            arrayView1d< integer const > const & elemGhostRank,
            arrayView2d< localIndex const > const & elemsToNodes,
            arrayView1d< globalIndex const > const & globalDofNumber,
            arrayView1d< R1Tensor const > const & disp,
            arrayView1d< R1Tensor const > const & uhat,
            arrayView1d< R1Tensor const > const & vtilde,
            arrayView1d< R1Tensor const > const & uhattilde,
            arrayView1d< real64 const > const & density,
            arrayView1d< real64 const > const & fluidPressure,
            arrayView1d< real64 const > const & deltaFluidPressure,
            arrayView1d< real64 const > const & biotCoefficient,
            timeIntegrationOption const tiOption,
            real64 const stiffnessDamping,
            real64 const massDamping,
            real64 const newmarkBeta,
            real64 const newmarkGamma,
            Epetra_FECrsMatrix * const dRdU,
            Epetra_FEVector * const residual );
  };
};



template< localIndex NUM_NODES_PER_ELEM, localIndex NUM_QUADRATURE_POINTS >
real64
SolidMechanicsLagrangianSSLE::ExplicitElementKernelWrapper::
Launch( localIndex const er,
        localIndex const esr,
        set<localIndex> const & elementList,
        arrayView2d<localIndex> const & elemsToNodes,
        arrayView3d< R1Tensor > const & dNdX,
        arrayView2d<real64> const & detJ,
        arrayView1d<R1Tensor> const & u,
        arrayView1d<R1Tensor> const & vel,
        arrayView1d<R1Tensor> & acc,
        ElementRegionManager::ConstitutiveRelationAccessor<constitutive::ConstitutiveBase> constitutiveRelations,
        ElementRegionManager::MaterialViewAccessor< arrayView2d<real64> > const & meanStress,
        ElementRegionManager::MaterialViewAccessor< arrayView2d<R2SymTensor> > const & devStress,
        real64 const dt )
{

  constitutive::ConstitutiveBase::UpdateFunctionPointer update = constitutiveRelations[er][esr][0]->GetStateUpdateFunctionPointer();
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

template< localIndex NUM_NODES_PER_ELEM, localIndex NUM_QUADRATURE_POINTS >
real64
SolidMechanicsLagrangianSSLE::
ImplicitElementKernelWrapper::Launch( localIndex const numElems,
                             real64 const dt,
                             arrayView3d<R1Tensor const> const & dNdX,
                             arrayView2d<real64 const > const& detJ,
                             FiniteElementBase const * const fe,
                             arrayView1d<constitutive::ConstitutiveBase *> const & constitutiveRelations,
                             arrayView1d< integer const > const & elemGhostRank,
                             arrayView2d< localIndex const > const & elemsToNodes,
                             arrayView1d< globalIndex const > const & globalDofNumber,
                             arrayView1d< R1Tensor const > const & disp,
                             arrayView1d< R1Tensor const > const & uhat,
                             arrayView1d< R1Tensor const > const & vtilde,
                             arrayView1d< R1Tensor const > const & uhattilde,
                             arrayView1d< real64 const > const & density,
                             arrayView1d< real64 const > const & fluidPressure,
                             arrayView1d< real64 const > const & deltaFluidPressure,
                             arrayView1d< real64 const > const & biotCoefficient,
                             timeIntegrationOption const tiOption,
                             real64 const stiffnessDamping,
                             real64 const massDamping,
                             real64 const newmarkBeta,
                             real64 const newmarkGamma,
                             Epetra_FECrsMatrix * const globaldRdU,
                             Epetra_FEVector * const globalResidual )
{
  constexpr int dim = 3;
  Epetra_LongLongSerialDenseVector  elementLocalDofIndex   (dim*static_cast<int>(NUM_NODES_PER_ELEM));
  Epetra_SerialDenseVector     R     (dim*static_cast<int>(NUM_NODES_PER_ELEM));
  Epetra_SerialDenseMatrix     dRdU  (dim*static_cast<int>(NUM_NODES_PER_ELEM),
                                                dim*static_cast<int>(NUM_NODES_PER_ELEM));
  Epetra_SerialDenseVector     element_dof_np1 (dim*static_cast<int>(NUM_NODES_PER_ELEM));

  Epetra_SerialDenseVector R_InertiaMassDamping(R);
  Epetra_SerialDenseMatrix dRdU_InertiaMassDamping(dRdU);
  Epetra_SerialDenseVector R_StiffnessDamping(R);
  Epetra_SerialDenseMatrix dRdU_StiffnessDamping(dRdU);

  real64 maxForce = 0;

  for( localIndex k=0 ; k<numElems ; ++k )
  {

    R1Tensor u_local[NUM_NODES_PER_ELEM];
    R1Tensor uhat_local[NUM_NODES_PER_ELEM];
    R1Tensor vtilde_local[NUM_NODES_PER_ELEM];
    R1Tensor uhattilde_local[NUM_NODES_PER_ELEM];

    dRdU.Scale(0);
    R.Scale(0);

    dRdU_InertiaMassDamping.Scale(0);
    R_InertiaMassDamping.Scale(0);
    dRdU_StiffnessDamping.Scale(0);
    R_StiffnessDamping.Scale(0);

    real64 c[6][6];
    constitutiveRelations[0]->GetStiffness( c );

    if(elemGhostRank[k] < 0)
    {
      for( localIndex a=0 ; a<NUM_NODES_PER_ELEM ; ++a)
      {

        localIndex localNodeIndex = elemsToNodes[k][a];

        for( int i=0 ; i<dim ; ++i )
        {
          elementLocalDofIndex[static_cast<int>(a)*dim+i] = dim*globalDofNumber[localNodeIndex]+i;

          // TODO must add last solution estimate for this to be valid
          element_dof_np1(static_cast<int>(a)*dim+i) = disp[localNodeIndex][i];
        }
      }

      if( tiOption == timeIntegrationOption::ImplicitDynamic )
      {
        GEOS_ERROR("Option not supported");
        CopyGlobalToLocal< NUM_NODES_PER_ELEM, R1Tensor>( elemsToNodes[k],
                                     disp, uhat, vtilde, uhattilde,
                                     u_local, uhat_local, vtilde_local, uhattilde_local );
      }
      else
      {
        CopyGlobalToLocal<NUM_NODES_PER_ELEM,R1Tensor>( elemsToNodes[k], disp, uhat, u_local, uhat_local );
      }

      R2SymTensor referenceStress;
      if( !fluidPressure.empty() )
      {
        referenceStress.PlusIdentity( - biotCoefficient[0] * (fluidPressure[k] + deltaFluidPressure[k]));
      }


      R1Tensor dNdXa;
      R1Tensor dNdXb;


      for( integer q=0 ; q<NUM_QUADRATURE_POINTS ; ++q )
      {
        const realT detJq = detJ[k][q];
        std::vector<double> const & N = fe->values(q);

        for( integer a=0 ; a<NUM_NODES_PER_ELEM ; ++a )
        {
    //      realT const * const dNdXa = dNdX(q,a).Data();
          dNdXa = dNdX[k][q][a];

          for( integer b=0 ; b<NUM_NODES_PER_ELEM ; ++b )
          {
    //        realT const * const dNdXb = dNdX(q,b).Data();
            dNdXb = dNdX[k][q][b];

            dRdU(a*dim+0,b*dim+0) -= ( c[0][0]*dNdXa[0]*dNdXb[0] + c[5][5]*dNdXa[1]*dNdXb[1] + c[4][4]*dNdXa[2]*dNdXb[2] ) * detJq;
            dRdU(a*dim+0,b*dim+1) -= ( c[5][5]*dNdXa[1]*dNdXb[0] + c[0][1]*dNdXa[0]*dNdXb[1] ) * detJq;
            dRdU(a*dim+0,b*dim+2) -= ( c[4][4]*dNdXa[2]*dNdXb[0] + c[0][2]*dNdXa[0]*dNdXb[2] ) * detJq;

            dRdU(a*dim+1,b*dim+0) -= ( c[0][1]*dNdXa[1]*dNdXb[0] + c[5][5]*dNdXa[0]*dNdXb[1] ) * detJq;
            dRdU(a*dim+1,b*dim+1) -= ( c[5][5]*dNdXa[0]*dNdXb[0] + c[1][1]*dNdXa[1]*dNdXb[1] + c[3][3]*dNdXa[2]*dNdXb[2] ) * detJq;
            dRdU(a*dim+1,b*dim+2) -= ( c[3][3]*dNdXa[2]*dNdXb[1] + c[1][2]*dNdXa[1]*dNdXb[2] ) * detJq;

            dRdU(a*dim+2,b*dim+0) -= ( c[0][2]*dNdXa[2]*dNdXb[0] + c[4][4]*dNdXa[0]*dNdXb[2] ) * detJq;
            dRdU(a*dim+2,b*dim+1) -= ( c[1][2]*dNdXa[2]*dNdXb[1] + c[3][3]*dNdXa[1]*dNdXb[2] ) * detJq;
            dRdU(a*dim+2,b*dim+2) -= ( c[4][4]*dNdXa[0]*dNdXb[0] + c[3][3]*dNdXa[1]*dNdXb[1] + c[2][2]*dNdXa[2]*dNdXb[2] ) * detJq;


            if( tiOption == timeIntegrationOption::ImplicitDynamic )
            {

              real64 integrationFactor = density[k] * N[a] * N[b] * detJq;
              real64 temp1 = ( massDamping * newmarkGamma/( newmarkBeta * dt ) + 1.0 / ( newmarkBeta * dt * dt ) )* integrationFactor;

              for( int i=0 ; i<dim ; ++i )
              {
                realT const acc = 1.0 / ( newmarkBeta * dt * dt ) * ( uhat[b][i] - uhattilde[b][i] );
                realT const velb = vtilde[b][i] + newmarkGamma/( newmarkBeta * dt ) *( uhat[b][i] - uhattilde[b][i] );

                dRdU_InertiaMassDamping(a*dim+i,b*dim+i) -= temp1;
                R_InertiaMassDamping(a*dim+i) -= ( massDamping * velb + acc ) * integrationFactor;
              }
            }
          }
        }
      }



        R1Tensor temp;
        for( integer q=0 ; q<NUM_QUADRATURE_POINTS ; ++q )
        {
          const realT detJq = detJ[k][q];
          R2SymTensor stress0 = referenceStress;
          stress0 *= detJq;
          for( integer a=0 ; a<NUM_NODES_PER_ELEM ; ++a )
          {
            dNdXa = dNdX[k][q][a];

            temp.AijBj(stress0,dNdXa);
            realT maxf = temp.MaxVal();
            if( maxf > maxForce )
            {
              maxForce = maxf;
            }

            R(a*dim+0) -= temp[0];
            R(a*dim+1) -= temp[1];
            R(a*dim+2) -= temp[2];
          }
        }


    // TODO It is simpler to do this...try it.
    //  dRdU.Multiply(dof_np1,R);
      for( integer a=0 ; a<NUM_NODES_PER_ELEM ; ++a )
      {
        realT nodeForce = 0;
        for( integer b=0 ; b<NUM_NODES_PER_ELEM ; ++b )
        {
          for( int i=0 ; i<dim ; ++i )
          {
            for( int j=0 ; j<dim ; ++j )
            {
              R(a*dim+i) += dRdU(a*dim+i,b*dim+j) * u_local[b][j];
            }
          }

          if( tiOption == timeIntegrationOption::ImplicitDynamic )
          {
            for( int i=0 ; i<dim ; ++i )
            {
              for( int j=0 ; j<dim ; ++j )
              {
                R_StiffnessDamping(a*dim+i) += stiffnessDamping * dRdU(a*dim+i,b*dim+j) * ( vtilde[b][j] + newmarkGamma/(newmarkBeta * dt)*(uhat[b][j]-uhattilde[b][j]) );
              }
            }
          }

        }

        nodeForce = std::max( std::max( R(a*dim+0), R(a*dim+1) ),  R(a*dim+2) );
        if( fabs(nodeForce) > maxForce )
        {
          maxForce = fabs(nodeForce);
        }
      }


      if( tiOption == timeIntegrationOption::ImplicitDynamic )
      {
        dRdU_StiffnessDamping = dRdU;
        dRdU_StiffnessDamping.Scale( stiffnessDamping * newmarkGamma / ( newmarkBeta * dt ) );

        dRdU += dRdU_InertiaMassDamping;
        dRdU += dRdU_StiffnessDamping;
        R    += R_InertiaMassDamping;
        R    += R_StiffnessDamping;
      }

      globaldRdU->SumIntoGlobalValues( elementLocalDofIndex,
                                   dRdU);

      globalResidual->SumIntoGlobalValues( elementLocalDofIndex,
                                R);
    }
  }
  return maxForce;
}

} /* namespace geosx */

#endif /* CORECOMPONENTS_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSLAGRANGIANSSLE_HPP_ */
