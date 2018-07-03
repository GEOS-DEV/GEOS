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
 * GEOSX is a free software; you can redistrubute it and/or modify it under
 * the terms of the GNU Lesser General Public Liscense (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file LagrangeDynamicsParallelPlateFlowExplicit.cpp
 * @author settgast1
 * @date Feb 28, 2011
 */

#include "LagrangeDynamicsParallelPlateFlowExplicit.h"
#include "SolverFactory.h"
#include "SurfaceGeneration/FractunatorBase.h"

LagrangeDynamicsParallelPlateFlowExplicit::LagrangeDynamicsParallelPlateFlowExplicit( const std::string& name,
                                                                                      ProblemManagerT* const pm ):
  SolverBase(name,pm),
  m_ldSolve(name,pm),
  m_ppSolve(name,pm)
{}

LagrangeDynamicsParallelPlateFlowExplicit::~LagrangeDynamicsParallelPlateFlowExplicit()
{
  // TODO Auto-generated destructor stub
}

void LagrangeDynamicsParallelPlateFlowExplicit::ReadXML( TICPP::HierarchicalDataNode* const hdn )
{
  SolverBase::ReadXML(hdn);

  m_ppSolve.ReadXML(hdn);
  m_ldSolve.ReadXML(hdn);

  m_kJn = hdn->GetAttributeOrDefault<realT>("normalJointStiffness", 1.0e10);
  m_kJs = hdn->GetAttributeOrDefault<realT>("shearJointStiffness", 1.0e10);
  m_COFJ = fabs(hdn->GetAttributeOrDefault<realT>("COFJoint", 0.5));
  m_fLockedInSIF =  hdn->GetAttributeOrDefault<realT>("lockedInSIFFactor", 1.0);

  m_faceStrengthRandomFactor = hdn->GetAttributeOrDefault<realT>("faceStrengthRandomFactor", 0.0);


}


void LagrangeDynamicsParallelPlateFlowExplicit::RegisterFields( PhysicalDomainT& domain )
{
  m_ppSolve.RegisterFields(domain);
  m_ldSolve.RegisterFields(domain);

  domain.m_feNodeManager.AddKeylessDataField<R1Tensor>("hydroForce", true, true);

  // TODO: remove the following and use interface logic
  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////
  domain.m_feFaceManager.AddKeylessDataField<realT>("initialContactStress", true, true);
  domain.m_feFaceManager.AddKeylessDataField<realT>("delta0N", true, false);
  domain.m_feFaceManager.AddKeylessDataField<R1Tensor>("gapShear0", true, false);
  domain.m_feFaceManager.AddKeylessDataField<R1Tensor>("stressShear0", true, false);
  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////
  domain.m_feFaceManager.AddKeylessDataField<R1Tensor>("contactStress", false, true);
  if (m_faceStrengthRandomFactor > 1e-99)
    domain.m_feFaceManager.AddKeylessDataField<realT>("faceStrengthRandomFactor", true, true);
  if (m_ppSolve.m_leakoffCoef > 0.0)
  {
    domain.m_feFaceManager.AddKeylessDataField<realT>("totalLeakedVolume", true, true);
    domain.m_feFaceManager.AddKeylessDataField<realT>("initialSaturatedTime", true, true);
  }

}

void LagrangeDynamicsParallelPlateFlowExplicit::Initialize(PhysicalDomainT& domain, SpatialPartition& partition)
{
  m_ppSolve.Initialize(domain,partition);
  m_ldSolve.Initialize(domain,partition);

  array<real64>& initialContactStress = domain.m_feFaceManager.GetFieldData<realT>("initialContactStress");
  array<real64>& faceContactStiffness = domain.m_feFaceManager.GetFieldData<realT>("faceContactStiffness");
  array<real64>* faceStrengthRandomFactor = domain.m_feFaceManager.GetFieldDataPointer<realT>("faceStrengthRandomFactor");
  faceContactStiffness = m_kJn;

  for( localIndex kf=0 ; kf<domain.m_feFaceManager.m_numFaces ; ++kf )
  {
    if( domain.m_feFaceManager.m_toElementsRelation[kf].size() > 1 )
    {
      R1Tensor fc;
      R1Tensor fn;

      R2SymTensor stress0;
      R2SymTensor stress1;

      R1Tensor t0, t1;
      R1Tensor temp;

      ElementRegionT& er0 = *(domain.m_feFaceManager.m_toElementsRelation[kf][0].first);
      ElementRegionT& er1 = *(domain.m_feFaceManager.m_toElementsRelation[kf][1].first);
      const localIndex elemIndex0 = domain.m_feFaceManager.m_toElementsRelation[kf][0].second;
      const localIndex elemIndex1 = domain.m_feFaceManager.m_toElementsRelation[kf][1].second;

      er0.m_mat->StateData(elemIndex0,0)->TotalStress(stress0);
      er1.m_mat->StateData(elemIndex1,0)->TotalStress(stress1);

      // normal from away from element 0, into element 1
      domain.m_feFaceManager.FaceCenterAndNormal( domain.m_feNodeManager, kf, fc, fn );

      t0.AijBj(stress0,fn);
      t1.AijBj(stress1,fn);

      realT t0n = Dot(t0,fn);
      realT t1n = Dot(t1,fn);

      // TODO: Add normal stress initialization routine for interfaces and call
      // here
      //////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////
      initialContactStress[kf] = 0.5 * fabs( t0n + t1n );
      //////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////
    }
  }
  if (faceStrengthRandomFactor!=NULL)
  {
    for( localIndex kf=0 ; kf<domain.m_feFaceManager.m_numFaces ; ++kf )
    {
      localIndex iNd = domain.m_feFaceManager.m_toNodesRelation[kf][0];
      realT xFace = (*domain.m_feNodeManager.m_refposition)[iNd][0];
      srand(int(xFace/0.001));
      (*faceStrengthRandomFactor)[kf] = 1.0 + ( (rand()*1.0) / RAND_MAX - 0.5) * m_faceStrengthRandomFactor;
    }
  }

  if (m_ppSolve.m_leakoffCoef > 0.0)
  {
    array<real64>& totalLeakedVolume = domain.m_feFaceManager.GetFieldData<realT>("totalLeakedVolume");
    totalLeakedVolume = 0.0;
    array<real64>& initialSaturatedTime = domain.m_feFaceManager.GetFieldData<realT>("initialSaturatedTime");
    initialSaturatedTime = std::numeric_limits<realT>::max();
  }


}

void LagrangeDynamicsParallelPlateFlowExplicit::InitializeCommunications( PartitionBase& partition )
{

  std::map<PhysicalDomainT::ObjectDataStructureKeys, array<string> > syncedFields;
  syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back(Field<FieldInfo::acceleration>::Name());
  syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back(Field<FieldInfo::velocity>::Name());
  syncedFields[PhysicalDomainT::FiniteElementFaceManager].push_back(Field<FieldInfo::mass>::Name());
  syncedFields[PhysicalDomainT::FiniteElementFaceManager].push_back("tracer");
  partition.SetBufferSizes( syncedFields, CommRegistry::lagrangeParallelPlateFlowSolver);
}


void LagrangeDynamicsParallelPlateFlowExplicit::ApplyForcesFromContact(PhysicalDomainT& domain,
                                                                       StableTimeStep& timeStep,
                                                                       const realT dt )
{
  //-----------------------------
  //external faces
  //-----------------------------
  {
    array<R1Tensor>& contactForce = domain.m_feNodeManager.GetFieldData<FieldInfo::contactForce>();
    contactForce = 0.0;
    array<R1Tensor>& decontactForce = domain.m_discreteElementSurfaceNodes.GetFieldData<FieldInfo::contactForce>();
    decontactForce = 0.0;

    //update nodal positions, velocities, and accelerations before updating face
    // geometry
    //also, reset rotational and translational accelerations as well as forces
    // and moments
    domain.m_discreteElementManager.UpdateNodalStatesZeroForcesAndAccelerations();

    //update face geometry and sort faces if necessary
    const bool resort = domain.m_externalFaces.RecalculateNeighborList(domain.m_feNodeManager,
                                                                       domain.m_discreteElementSurfaceNodes,
                                                                       domain.m_discreteElementManager,
                                                                       dt);

    //if a resort has been triggered, then you also need to update the contact
    // manager
    if (resort)
      domain.m_contactManager.Update(domain.m_externalFaces.m_neighborList);

    {
      array<array<R1Tensor> > xs;
      xs.resize(domain.m_externalFaces.DataLengths());
      domain.m_externalFaces.UpdateGeometricContactProperties(dt, domain, xs);
#ifdef STATES_ON_CONTACTS
      domain.m_externalFaces.UpdateAndApplyContactForcesFromStatesOnCommonPlanes(dt, domain);
#else
      domain.m_externalFaces.UpdateAndApplyContactStresses(timeStep, dt, domain, xs);
#endif
    }

    //for parallel: do NOT need to synchronize nodal force fields across
    // processes for DE nodes
    //before moving to centroid, since each process will do that calculation
    // (redundantly) itself
    // ... there's therefore no need for explicit synchrony
    //as long as nodal states remain synchronous, everything will be fine!
#ifndef DEEFC
    domain.m_discreteElementManager.ApplyDiscreteElementContactForces(45, 0.2);
#else
    domain.m_discreteElementManager.ApplyNodalForces();
#endif
  }
}


double LagrangeDynamicsParallelPlateFlowExplicit::TimeStep( const realT& time,
                                                            const realT& dt,
                                                            const int cycleNumber,
                                                            PhysicalDomainT& domain,
                                                            const array<string>& namesOfSolverRegions,
                                                            SpatialPartition& partition,
                                                            FractunatorBase* const fractunator )
{
  realT dt_return = dt;

  if( fractunator!=NULL )
  {
    if( cycleNumber%fractunator->m_checkInterval==0 )
    {
      fractunator->SeparationDriver( domain.m_feNodeManager,
                                     domain.m_feEdgeManager,
                                     domain.m_feFaceManager,
                                     domain.m_externalFaces,
                                     domain.m_feElementManager,
                                     partition, false, time);
    }
  }

  m_stabledt.m_maxdt = std::numeric_limits<double>::max();
  array<R1Tensor>& hgforce = domain.m_feNodeManager.GetFieldData<FieldInfo::hgforce> ();
  hgforce = 0.0;

  array<integer>& flowFaceType = domain.m_feFaceManager.GetFieldData<int>("flowFaceType");
//  array<lArray1d>& edgeToFlowFaces =
// domain.m_edgeManager.GetVariableOneToManyMap("edgeToFlowFaces");

  array<R1Tensor>* u0 = domain.m_feNodeManager.GetFieldDataPointer<R1Tensor>("displacement0");
  array<R1Tensor>* unet = domain.m_feNodeManager.GetFieldDataPointer<R1Tensor>("netDisplacement");

  m_ppSolve.GenerateParallelPlateGeometricQuantities( domain,time,dt );

  m_ppSolve.CalculateAndApplyMassFlux( dt, domain );
  if (m_ppSolve.m_leakoffCoef > 0.0)
    m_ppSolve.CalculateCarterLeakOff(time, dt, domain);
  m_ppSolve.ApplyFluxBoundaryCondition(time, dt, cycleNumber, partition.m_rank, domain);



  // update nodes
  LagrangeExplicitDynamicsFunctions::LinearPointUpdatePart1( domain.m_feNodeManager, time, dt );

  for( std::map< ElementManagerT::RegKeyType, ElementRegionT >::iterator elementRegionIter = domain.m_feElementManager.m_ElementRegions.begin() ;
       elementRegionIter != domain.m_feElementManager.m_ElementRegions.end() ;
       ++elementRegionIter )
  {
    ElementRegionT& elementRegion = elementRegionIter->second;

    elementRegion.CalculateVelocityGradients(domain.m_feNodeManager);

    elementRegion.MaterialUpdate(dt);

    elementRegion.CalculateNodalForces(domain.m_feNodeManager, m_stabledt, dt);
  }


  BoundaryConditionFunctions::ApplyTractionBoundaryCondition( domain,time);

  // -- CONTACT --
  if (domain.m_externalFaces.m_contactActive)
    ApplyForcesFromContact(domain, this->m_stabledt, dt);

  m_stabledt.m_maxdt *= m_ldSolve.m_courant;

  if (cycleNumber%10000 == 0)
  {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::cout << "t=" << time << "rank: " << rank << "Solid dt: " << m_stabledt.m_maxdt<<", Flow dt: "<<m_ppSolve.m_stabledt.m_maxdt << std::endl;
  }

//  std::cout<<m_stabledt.m_maxdt<<", "<<m_ppSolve.m_stabledt.m_maxdt;

  if(  m_stabledt.m_maxdt > m_ppSolve.m_stabledt.m_maxdt )
    m_stabledt.m_maxdt = m_ppSolve.m_stabledt.m_maxdt;

  // -- CONTACT --
//  if (domain.m_externalFaces.m_contactActive)
//    m_ldSolve.ApplyForcesFromContact(domain, dt);


  m_ppSolve.GenerateParallelPlateGeometricQuantities( domain, time,dt );



  {
    std::map<PhysicalDomainT::ObjectDataStructureKeys, array<string> > syncedFields;

    syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back(Field<FieldInfo::acceleration>::Name());
    syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back(Field<FieldInfo::velocity>::Name());
    syncedFields[PhysicalDomainT::FiniteElementFaceManager].push_back(Field<FieldInfo::mass>::Name());
    syncedFields[PhysicalDomainT::FiniteElementFaceManager].push_back("tracer");
    partition.SynchronizeFields( syncedFields, CommRegistry::lagrangeParallelPlateFlowSolver );
  }

  m_ppSolve.UpdateEOS( time, dt, domain );


  // apply fluid pressure to faces
  const OrderedVariableOneToManyRelation& childFaceIndex = domain.m_feFaceManager.GetVariableOneToManyMap( "childIndices" );
  const array<real64>& faceFluidPressure = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();
//  const array<real64>& initialContactStress =
// domain.m_feFaceManager.GetFieldData<realT>("initialContactStress");
  array<R1Tensor>& nodalForce = domain.m_feNodeManager.GetFieldData<FieldInfo::force>();
  array<R1Tensor>& hydroForce = domain.m_feNodeManager.GetFieldData<R1Tensor>("hydroForce");
  array<R1Tensor>& contactForce = domain.m_feNodeManager.GetFieldData<FieldInfo::contactForce>();
  array<R1Tensor>& contactStress = domain.m_feFaceManager.GetFieldData<R1Tensor>("contactStress");

//  array<R1Tensor>& cohesiveForce =
// domain.m_feNodeManager.GetFieldData<R1Tensor>("cohesiveForce");

//  const array<R1Tensor>& cohesiveTraction =
// domain.m_feFaceManager.GetFieldData<R1Tensor>("cohesiveTraction");

  hydroForce = 0.0;
  contactForce = 0.0;
//  cohesiveForce = 0.0;

  for( localIndex kf=0 ; kf<domain.m_feFaceManager.DataLengths() ; ++kf )
  {
    if( flowFaceType[kf] == 1 && childFaceIndex[kf].size() >= 1 )
    {

      realT pressure = faceFluidPressure[kf];

//      if( pressure < initialContactStress[kf] )
//      {
//        pressure = initialContactStress[kf];
//      }

      localIndex faceIndex[2];
      if (childFaceIndex[kf].size() == 1)
      {
        faceIndex[0] = kf;
        faceIndex[1] = childFaceIndex[kf][0];
      }
      else
      {
        faceIndex[0] = childFaceIndex[kf][0];
        faceIndex[1] = childFaceIndex[kf][1];
      }

      realT stressPen;
      R1Tensor stressShear;
      CalculateContactStress(domain, time, dt, kf, faceIndex, stressPen, stressShear);


      //////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////

      const R1Tensor N[2] = { domain.m_feFaceManager.FaceNormal( domain.m_feNodeManager, faceIndex[0] ),
                              domain.m_feFaceManager.FaceNormal( domain.m_feNodeManager, faceIndex[1] )};

      for( int side=0 ; side<2 ; ++side )
      {
        const localIndex faceID = faceIndex[side];
        const realT area = domain.m_feFaceManager.SurfaceArea( domain.m_feNodeManager, faceID );

        R1Tensor hForce = N[side];
        R1Tensor pForce = N[side];
        R1Tensor cForce, sForce;
        int direction = -1;
        if(side==0)
          direction = 1;

        hForce *= -area * pressure / domain.m_feFaceManager.m_toNodesRelation[faceID].size();
        pForce *= -area * stressPen / domain.m_feFaceManager.m_toNodesRelation[faceID].size();
        sForce = stressShear;
        contactStress[faceID] = pForce / area * domain.m_feFaceManager.m_toNodesRelation[faceID].size();
        contactStress[faceID] += sForce * direction;
        sForce *= direction * area / domain.m_feFaceManager.m_toNodesRelation[faceID].size();

//        cForce = cohesiveTraction[kf];
        cForce *= direction * area / domain.m_feFaceManager.m_toNodesRelation[faceID].size();


        for( lArray1d::const_iterator nodeID=domain.m_feFaceManager.m_toNodesRelation[faceID].begin() ;
             nodeID!=domain.m_feFaceManager.m_toNodesRelation[faceID].end() ; ++nodeID )
        {
          nodalForce[*nodeID] += hForce;
          nodalForce[*nodeID] += pForce;
          nodalForce[*nodeID] += cForce;
          nodalForce[*nodeID] += sForce;

          hydroForce[*nodeID] += hForce;
          contactForce[*nodeID] += pForce;
          contactForce[*nodeID] += sForce;

//          cohesiveForce[*nodeID] += cForce;
        }

      }
    }

  }

//  m_ldSolve.ApplyGapDamping( domain.m_nodeManager,
//                            domain.m_faceManager, dt );



  LagrangeExplicitDynamicsFunctions::LinearPointUpdatePart2( domain.m_feNodeManager, time, dt, m_ldSolve.m_dampingM );

  if (unet != NULL && u0 != NULL)
  {
    for (localIndex i = 0 ; i < domain.m_feNodeManager.DataLengths() ; ++i)
    {
      (*unet)[i] = (*domain.m_feNodeManager.m_displacement)[i];
      (*unet)[i] -= (*u0)[i];
    }
  }

  {
    std::map<PhysicalDomainT::ObjectDataStructureKeys, array<string> > syncedFields;

    syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back(Field<FieldInfo::acceleration>::Name());
    syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back(Field<FieldInfo::velocity>::Name());
    syncedFields[PhysicalDomainT::FiniteElementFaceManager].push_back(Field<FieldInfo::mass>::Name());
    syncedFields[PhysicalDomainT::FiniteElementFaceManager].push_back("tracer");
    {
      array<real64>* initialSaturatedTime = domain.m_feFaceManager.GetFieldDataPointer<realT>("initialSaturatedTime");
      if (initialSaturatedTime != NULL )
        syncedFields[PhysicalDomainT::FiniteElementFaceManager].push_back("initialSaturatedTime");
    }


    partition.SynchronizeFields( syncedFields, CommRegistry::lagrangeParallelPlateFlowSolver );
  }
  return dt_return;

}

void LagrangeDynamicsParallelPlateFlowExplicit::CalculateContactStress(PhysicalDomainT& domain,
                                                                       const realT& time,
                                                                       const realT& dt,
                                                                       localIndex& kf,
                                                                       const localIndex faceIndex[],
                                                                       realT& stressPen,
                                                                       R1Tensor& stressShear)
{
  stressPen = 0;
  stressShear = 0;
  array<real64>& delta0N = domain.m_feFaceManager.GetFieldData<realT>("delta0N");
  array<R1Tensor>& gapShear0 = domain.m_feFaceManager.GetFieldData<R1Tensor>("gapShear0");
  array<R1Tensor>& stressShear0 = domain.m_feFaceManager.GetFieldData<R1Tensor>("stressShear0");
  array<real64>& effectiveStressN = domain.m_feFaceManager.GetFieldData<realT>("effectiveStressN");


  const R1Tensor N[2] = { domain.m_feFaceManager.FaceNormal( domain.m_feNodeManager, faceIndex[0] ),
                          domain.m_feFaceManager.FaceNormal( domain.m_feNodeManager, faceIndex[1] )};

  R1Tensor T[2];
  domain.m_feFaceManager.FaceTangential(domain.m_feNodeManager, faceIndex[0], T[0], T[1]);
  R1Tensor Nbar = N[0];
  Nbar -= N[1];
  Nbar.Normalize();

  const R1Tensor gap = domain.m_feFaceManager.CalculateGapVector( domain.m_feNodeManager, kf );
  const realT gapNormal = Dot(gap, Nbar);
  R1Tensor gapShear = Dot(gap, T[0]) * T[0] + Dot(gap, T[1]) * T[1];
  R1Tensor gapShearInc(gapShear);
  gapShearInc -= gapShear0[kf];


  realT deltaN;
  realT tBuffer = dt * 500;
  if (time == 0)
  {
    deltaN = delta0N[kf];
  }
  else if (time < tBuffer )
  {
    deltaN = delta0N[kf] * (1.0 - time / tBuffer * m_fLockedInSIF);
  }
  else
    deltaN = (1-m_fLockedInSIF) * delta0N[kf];

  if( -gapNormal + deltaN <= 0.0 )
  {
    stressPen = 0.0;
    stressShear = 0.0;
  }
  else
  {
    stressPen = (-gapNormal + deltaN) * m_kJn;

    R1Tensor stressShearPrj = stressShear0[kf];
    stressShearPrj += gapShearInc * m_kJs;

    if ( pow(stressShearPrj.L2_Norm(),2) < pow(stressPen*m_COFJ, 2) )
    {
      stressShear = stressShearPrj;
    }
    else  // We have to calculate plasticity stuff
    {
      realT du[2], t[2], x[2];
      du[0] = Dot(gapShearInc, T[0]);
      du[1] = Dot(gapShearInc, T[1]);
      t[0] = Dot(stressShear0[kf], T[0]);
      t[1] = Dot(stressShear0[kf], T[1]);

      realT a, b, c;
      a = pow(t[0],2) + pow(t[1],2);
      if (a == 0)
        a = pow( t[0]+du[0]*m_kJs,2) + pow(t[1]+ du[1]*m_kJs,2);

      b = -2 * ( t[0]*(t[0] + m_kJs * du[0]) + t[1]*(t[1] + m_kJs * du[1]));
      c = pow(t[0] + m_kJs * du[0], 2) + pow(t[1] + m_kJs * du[1], 2) - pow(stressPen*m_COFJ, 2);

      realT b2_4ac = pow(b,2)- 4 * a * c;

      if (b2_4ac < 0.0 || a==0)
      {
        t[0] = 0;
        t[1] = 0;
      }
      else
      {
        x[0] = (-b - pow(b2_4ac, 0.5))/2/a;
        x[1] = (-b + pow(b2_4ac, 0.5))/2/a;

        if (fabs(x[0]) > fabs(x[1]))
          x[0] = x[1];

        t[0] += m_kJs * du[0] - x[0] * t[0];
        t[1] += m_kJs * du[1] - x[0] * t[1];
      }
      stressShear = t[0] * T[0] + t[1] * T[1];

    }
  }
  effectiveStressN[kf] = stressPen;
  stressShear0[kf] = stressShear;
  gapShear0[kf] = gapShear;
}



void LagrangeDynamicsParallelPlateFlowExplicit::PostProcess (PhysicalDomainT& domain,
                                                             SpatialPartition& partition,
                                                             const array<string>& namesOfSolverRegions)
{
  m_ppSolve.PostProcess(domain, partition, namesOfSolverRegions);
  m_ldSolve.PostProcess(domain, partition, namesOfSolverRegions);
}

/// Register solver in the solver factory
REGISTER_SOLVER( LagrangeDynamicsParallelPlateFlowExplicit )
