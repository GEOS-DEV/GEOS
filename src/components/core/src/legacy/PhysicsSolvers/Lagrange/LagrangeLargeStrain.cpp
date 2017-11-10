/*
 * LagrangeLargeStrain.cpp
 *
 *  Created on: Nov 8, 2012
 *      Author: settgast1
 */

#include "PhysicsSolvers/Lagrange/LagrangeHelperFunctions.h"
#include "LagrangeLargeStrain.h"
#include "PhysicsSolvers/SolverFactory.h"
//#include "PhysicsSolvers/PhysicsSolverStrings.h"
#include "DataStructures/VectorFields/ElementRegionT.h"
#include "ElementLibrary/FiniteElement.h"
#include "ElementLibrary/FiniteElementUtilities.h"

#include "Utilities/Kinematics.h"


LagrangeLargeStrain::LagrangeLargeStrain(  const std::string& name,
                                                ProblemManagerT* const pm ):
LagrangeSolverBase(name,pm)
{
  // TODO Auto-generated constructor stub

}


LagrangeLargeStrain::~LagrangeLargeStrain()
{
  // TODO Auto-generated destructor stub
}




void LagrangeLargeStrain::ProcessElementRegion( NodeManager& nodeManager,
                                                     ElementRegionT& elemRegion,
                                                     const realT dt )
{

  R2Tensor A;
  R2Tensor F;
  R2Tensor Finv;
  R2Tensor dUhatdX;
  R2SymTensor totalStress;
  R2SymTensor Dadt;
  R2Tensor Rot;


  static array< R1Tensor > u_local;
  static array< R1Tensor > uhat_local;
  static array<R1Tensor> f_local;
  static array< R1Tensor > x;
  static array< R1Tensor > v;
  static array<R1Tensor> s_dNdx;
  static array<R1Tensor> f_zemc;
  static array<R1Tensor> Q;


  rArray2d& detJ = elemRegion.m_detJ;
  rArray2d& detJ_n = elemRegion.m_detJ_n;
  rArray2d& detJ_np1 = elemRegion.m_detJ_np1;

  Array2dT< R2Tensor >&  dUdX = elemRegion.m_dUdX;
  FiniteElementBase*& finiteElement = elemRegion.m_finiteElement;

  //  MaterialBaseT*& m_materialComputations = elemRegion.m_materialComputations;

  const unsigned int numNodesPerElem = elemRegion.m_numNodesPerElem;

  if( u_local.size() != numNodesPerElem )
  {
    u_local.resize(numNodesPerElem);
    uhat_local.resize(numNodesPerElem);
    f_local.resize(numNodesPerElem);
    x.resize(numNodesPerElem);
    v.resize(numNodesPerElem);
    s_dNdx.resize(numNodesPerElem);
    f_zemc.resize(numNodesPerElem);
    Q.resize(finiteElement->zero_energy_modes());
  }


  const array<R1Tensor>& incrementalDisplacement = nodeManager.GetFieldData<FieldInfo::incrementalDisplacement>();
  const array<R1Tensor>& totalDisplacement = nodeManager.GetFieldData<FieldInfo::displacement>();

//  const array<integer>& attachedToSendingGhostNode = GetFieldData<int>("attachedToSendingGhostNode");


  array<real64>& volume = elemRegion.GetFieldData<FieldInfo::volume>();
  array<real64>& volume_n = elemRegion.GetFieldData<realT>("volume_n");

  volume_n = volume;
  detJ_n = detJ_np1;








  const array<R1Tensor>& referencePosition = nodeManager.GetFieldData<FieldInfo::referencePosition>();
  const array<R1Tensor>& velocity = nodeManager.GetFieldData<FieldInfo::velocity>();
  array<R1Tensor>& force = nodeManager.GetFieldData<FieldInfo::force>();
  array<R1Tensor>& hgforce = nodeManager.GetFieldData<FieldInfo::hgforce> ();

//  hgforce = 0.0;
  array<array<R1Tensor>*> Qstiffness(elemRegion.m_finiteElement->zero_energy_modes(),NULL);
  if( elemRegion.m_finiteElement->zero_energy_modes() >= 1 )
  {
    Qstiffness[0] = &(elemRegion.GetFieldData<R1Tensor>("Qhg1"));
  }
  if( elemRegion.m_finiteElement->zero_energy_modes() >= 2 )
  {
    Qstiffness[1] = &(elemRegion.GetFieldData<R1Tensor>("Qhg2"));
  }
  if( elemRegion.m_finiteElement->zero_energy_modes() >= 3 )
  {
    Qstiffness[2] = &(elemRegion.GetFieldData<R1Tensor>("Qhg3"));
  }
  if( elemRegion.m_finiteElement->zero_energy_modes() >= 4 )
  {
    Qstiffness[3] = &(elemRegion.GetFieldData<R1Tensor>("Qhg4"));
  }

//  const array<integer>& ghostRank = elemRegion.GetFieldData<FieldInfo::ghostRank>();



  elemRegion.m_energy.Zero();




  for( localIndex k=0 ; k<elemRegion.m_numElems ; ++k )
  {
//    if( ghostRank[k] < 0 )
    elemRegion.GetElementCenter(k,nodeManager);
    {

      const localIndex* const elemToNodeMap = elemRegion.m_toNodesRelation[k];


      CopyGlobalToLocal( elemToNodeMap,
                         incrementalDisplacement, totalDisplacement,
                         uhat_local, u_local );


      if( elemRegion.m_finiteElement->zero_energy_modes() )
      {
        CopyGlobalToLocal( elemToNodeMap,
                           referencePosition, velocity,
                           x, v );

        x += u_local;
      }


      volume[k] = 0.0;

      f_local = 0.0;

//      if(0)

      const localIndex paramIndex = elemRegion.m_mat->NumParameterIndex0() > 1 ? k : 0 ;
      const MaterialBaseParameterData& param = *( elemRegion.m_mat->ParameterData(paramIndex) );
      s_dNdx = 0.0;
      realT initVolume = 0.0;

      for( unsigned int q=0 ; q<elemRegion.m_numIntegrationPointsPerElem ; ++q )
      {
        MaterialBaseStateData& state = *(elemRegion.m_mat->StateData(k,q));
        R1Tensor* const dNdX = elemRegion.m_dNdX(k)[q];

        // Velocity Gradient

        // calculate dUhat/dX at beginning of step
        CalculateGradient( dUhatdX ,uhat_local, dNdX );

        // calculate velocity gradient (mid-step)
        R2Tensor L;
	if(dt > 0)
        {
          // calculate dv/dX
          R2Tensor dvdX = dUhatdX;
          dvdX *= 1.0 / dt;

          // calculate du/dX
          F = dUhatdX;
          F *= 0.5;
          F += dUdX(k,q);
          F.PlusIdentity(1.0);

          // calculate dX/du
          Finv.Inverse(F);

          // chain rule: calculate dv/du = dv/dX * dX/du
          L.AijBjk(dvdX, Finv);
        }

        // calculate gradient (end of step)
        dUdX(k,q) += dUhatdX;
        F = dUdX(k,q);
        F.PlusIdentity(1.0);
        realT detF = F.Det();

        // calculate element volume
        detJ_np1(k,q) = detJ(k,q) * detF;
        volume[k] += detJ_np1(k,q);
        initVolume += detJ(k,q);

        Finv.Inverse(F);

        // Calculate Rate of deformation tensor and rotationAxis tensor
        A.AijBjk(dUhatdX,Finv);
        IncrementalKinematics(A,Dadt,Rot);

        const realT rho = param.init_density / fabs(detF);

        // update state before exercising material model
        if(elemRegion.m_mat->NeedsDensity() ){
          state.SetDensity(rho);
        }
        if(elemRegion.m_mat->NeedsSpecificInternalEnergy() ){

          state.TotalStress(totalStress);
          realT ie = state.GetSpecificInternalEnergy();
          realT die = L(0,0) * totalStress(0,0)
                    + (L(0,1) + L(1,0)) * totalStress(0,1)
                    + (L(0,2) + L(2,0)) * totalStress(0,2)
                    + L(1,1)*totalStress(1,1)
                    + (L(1,2) + L(2,1)) * totalStress(1,2)
                    + L(2,2)*totalStress(2,2);
           die *= -dt/(rho+1e-64); // negative because stress is +ve in compression
           ie += die;
          state.SetSpecificInternalEnergy(ie);
        }

        // Material Update

        elemRegion.m_mat->StrainDrivenUpdateMember( k, q,
                                                    Dadt, L, Rot,
                                                    detJ_n[k][q],
                                                    detJ_np1[k][q],
                                                    dt);

        // nodal force calculation
        state.TotalStress(totalStress);

        //const realT rho = param.init_density / fabs(detF);
        const realT soundSpeed = sqrt( state.BulkModulus / rho );
        const realT trD = dt > 0 ? Dadt.Trace() / dt : 0.0;


        realT bulkQ = LagrangeHelperFunctions::BulkQ( rho,
                                                soundSpeed,
                                                this->m_bulkQLinear,
                                                this->m_bulkQQuadratic,
                                                trD,
                                                cbrt( volume[k] ) );

        totalStress.PlusIdentity( bulkQ );
        FiniteElementUtilities::Integrate( totalStress,
                                           elemRegion.m_dNdX(k)[q],
                                           detJ(k,q),
                                           detF,
                                           Finv,
                                           f_local.size(),
                                           f_local.data() );

        for( unsigned int b=0 ; b<numNodesPerElem ; ++b )
        {
          R1Tensor temp;
          temp.AijBi(Finv,dNdX[b]);
          s_dNdx(b) += temp;
//          BB += Dot( s_dNdx(b), s_dNdx(b) ) ;
        }
      }
      s_dNdx /= elemRegion.m_numIntegrationPointsPerElem;
      realT BB = 0.0;
      for( unsigned int b=0 ; b<numNodesPerElem ; ++b )
      {
        BB += Dot( s_dNdx(b), s_dNdx(b) ) ;
      }

      realT thisdt =  LagrangeHelperFunctions::CalculateMaxStableExplicitTimestep( param.init_density / ( volume[k]/initVolume ),
                                                                                  param.Lame + 2*param.init_shearModulus,
                                                                                  BB );

      if( elemRegion.m_ElementDimension == 3 )
      {
        thisdt /= sqrt(2.0);
      }

      if( thisdt < SolverBase::m_stabledt.m_maxdt )
      {
//        timeStep.m_region = this->m_regionName;
//        timeStep.m_index = k;
        SolverBase::m_stabledt.m_maxdt = thisdt;
      }


      if( finiteElement->zero_energy_modes() )
      {

        for( int m=0 ; m<finiteElement->zero_energy_modes() ; ++m )
        {
          Q[m] = (*(Qstiffness[m]))[k];
        }

        finiteElement->zero_energy_mode_control( s_dNdx, volume[k], x, v,
                                                   elemRegion.m_hgDamp,
                                                   elemRegion.m_hgStiff*dt,
                                                   param.init_density,
                                                   param.Lame + 2*param.init_shearModulus ,
                                                   dt, Q, f_zemc );
        for( int m=0 ; m<finiteElement->zero_energy_modes() ; ++m )
        {
          (*(Qstiffness[m]))[k] = Q[m];
        }

        AddLocalToGlobal( elemToNodeMap,
                          f_zemc, f_zemc,
                          force, hgforce);

      }

      AddLocalToGlobal( elemToNodeMap, f_local, force);
    }
  }

}



void LagrangeLargeStrain::ApplyThermalStress( ElementRegionT& elemRegion,
                                              NodeManager& nodeManager,
                                              const localIndex& elementID,
                                              Epetra_SerialDenseVector * rhs)
{

  array<real64>& temperature = elemRegion.GetFieldData<realT>("temperature");
  array<real64>& CTE = elemRegion.GetFieldData<realT>("linearCTE");
  array<real64>& antiThermalStress = elemRegion.GetFieldData<realT>("antiThermalStress");
  array<real64>* refTemperature = elemRegion.GetFieldDataPointer<realT>("refTemperature");


  if (m_useNodalTemperature > 0)
  {
    const array<real64>& nodalTemp = nodeManager.GetFieldData<realT>("Temperature");
    temperature[elementID] = 0.0;
    for (localIndex i = 0; i<elemRegion.m_toNodesRelation.Dimension(1); ++i)
    {
      temperature[elementID] += nodalTemp[elemRegion.m_toNodesRelation[elementID][i]];
    }
    temperature[elementID] /= elemRegion.m_toNodesRelation.Dimension(1);
  }

  {
    const localIndex paramIndex = elemRegion.m_mat->NumParameterIndex0() > 1 ? elementID : 0 ;

    const MaterialBaseParameterData& matParams = *(elemRegion.m_mat->ParameterData(paramIndex));

    const realT K =  matParams.Lame + matParams.init_shearModulus * 2.0 / 3.0;

    R2SymTensor thermalStress;
    thermalStress *= 0.0;

    if (!refTemperature)
    {
      for (int i = 0; i < dim; ++i) thermalStress(i,i) = -(temperature[elementID] - m_refTemperature) * CTE [elementID] * K * 3;  //The factor 3 is because CTE is the linear CTE
    }
    else
    {
      for (int i = 0; i < dim; ++i) thermalStress(i,i) = -(temperature[elementID] - (*refTemperature)[elementID]) * CTE [elementID] * K * 3;  //The factor 3 is because CTE is the linear CTE
    }
    antiThermalStress[elementID] = -thermalStress(0,0);

    CalculateNodalForceFromStress( nodeManager, elemRegion, elementID, thermalStress);

  }

}

void LagrangeLargeStrain::CalculateNodalForceFromStress( NodeManager& nodeManager,
                                                         ElementRegionT& elemRegion,
                                                         const localIndex& elementID,
                                                         R2SymTensor& stress )
{
  R2Tensor F;
  R2Tensor Finv;
  R2Tensor dUhatdX;


  static array< R1Tensor > u_local;
  static array< R1Tensor > uhat_local;
  static array<R1Tensor> f_local;

  rArray2d& detJ = elemRegion.m_detJ;

  Array2dT< R2Tensor >&  dUdX = elemRegion.m_dUdX;

  const unsigned int numNodesPerElem = elemRegion.m_numNodesPerElem;

  if( u_local.size() != numNodesPerElem )
  {
    u_local.resize(numNodesPerElem);
    uhat_local.resize(numNodesPerElem);
    f_local.resize(numNodesPerElem);
  }


  const array<R1Tensor>& incrementalDisplacement = nodeManager.GetFieldData<FieldInfo::incrementalDisplacement>();
  const array<R1Tensor>& totalDisplacement = nodeManager.GetFieldData<FieldInfo::displacement>();


  array<R1Tensor>& force = nodeManager.GetFieldData<FieldInfo::force>();

  const localIndex* const elemToNodeMap = elemRegion.m_toNodesRelation[elementID];


  CopyGlobalToLocal( elemToNodeMap,
                     incrementalDisplacement, totalDisplacement,
                     uhat_local, u_local );

  f_local = 0.0;

  for( unsigned int a=0 ; a<elemRegion.m_numIntegrationPointsPerElem ; ++a )
  {

    R1Tensor* const dNdX = elemRegion.m_dNdX(elementID)[a];

    // Velocity Gradient

    // calculate dUhat/dX at beginning of step
    CalculateGradient( dUhatdX ,uhat_local, dNdX );


    // calculate gradient (end of step)
    dUdX(elementID,a) += dUhatdX;
    F = dUdX(elementID,a);
    F.PlusIdentity(1.0);
    realT detF = F.Det();


    Finv.Inverse(F);
    FiniteElementUtilities::Integrate( stress,
                                       elemRegion.m_dNdX(elementID)[a],
                                       detJ(elementID,a),
                                       detF,
                                       Finv,
                                       f_local.size(),
                                       f_local.data() );
  }


  AddLocalToGlobal( elemToNodeMap, f_local, force);


}

/* Register solver in the solver factory */

SolverRegistrator<LagrangeLargeStrain> reg_LagrangeLargeStrain;
