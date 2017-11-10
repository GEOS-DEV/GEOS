//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Copyright (c) 2015, Lawrence Livermore National Security, LLC.
//  Produced at the Lawrence Livermore National Laboratory
//
//  GEOS Computational Framework - Core Package, Version 3.0.0
//
//  Written by:
//  Randolph Settgast (settgast1@llnl.gov)
//  Stuart Walsh(walsh24@llnl.gov)
//  Pengcheng Fu (fu4@llnl.gov)
//  Joshua White (white230@llnl.gov)
//  Chandrasekhar Annavarapu Srinivas
//  Eric Herbold
//  Michael Homel
//
//
//  All rights reserved.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
//  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
//  LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
//  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED 
//  AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
//  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//  1. This notice is required to be provided under our contract with the U.S. Department of Energy (DOE). This work was produced at Lawrence Livermore 
//     National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
//  2. Neither the United States Government nor Lawrence Livermore National Security, LLC nor any of their employees, makes any warranty, express or 
//     implied, or assumes any liability or responsibility for the accuracy, completeness, or usefulness of any information, apparatus, product, or 
//     process disclosed, or represents that its use would not infringe privately-owned rights.
//  3. Also, reference herein to any specific commercial products, process, or services by trade name, trademark, manufacturer or otherwise does not 
//     necessarily constitute or imply its endorsement, recommendation, or favoring by the United States Government or Lawrence Livermore National Security, 
//     LLC. The views and opinions of authors expressed herein do not necessarily state or reflect those of the United States Government or Lawrence 
//     Livermore National Security, LLC, and shall not be used for advertising or product endorsement purposes.
//
//  This Software derives from a BSD open source release LLNL-CODE-656616. The BSD  License statment is included in this distribution in src/bsd_notice.txt.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @file LagrangeExplicitDynamicsSolver.cpp
 * @author Randolph Settgast
 * @date created on Sep 13, 2010
 */

#include "LagrangeExplicitDynamicsSolver.h"
#include "SolverFactory.h"
#include "Utilities/Utilities.h"
#include "Utilities/Kinematics.h"
#include "SurfaceGeneration/FractunatorBase.h"

#include "ObjectManagers/TableManager.h"


LagrangeExplicitDynamicsSolver::LagrangeExplicitDynamicsSolver( const std::string& name,
                                                                ProblemManagerT* const pm ):
SolverBase(name,pm)
{
}

LagrangeExplicitDynamicsSolver::~LagrangeExplicitDynamicsSolver()
{
  // TODO Auto-generated destructor stub
}


void LagrangeExplicitDynamicsSolver::ReadXML( TICPP::HierarchicalDataNode* hdn )
{
  SolverBase::ReadXML(hdn);

  m_courant = hdn->GetAttributeValue<realT>("courant");
  m_dampingM = hdn->GetAttributeOrDefault<realT>("dampingM",0.0);
  m_gapdamping = hdn->GetAttributeOrDefault<realT>("gapdamping",0);

  m_tiedNodesFlag = hdn->GetAttributeOrDefault<int> ("tiedNodesFlag", 0);
  m_tiedNodeNormalRuptureStress = hdn->GetAttributeOrDefault<realT> ("tiedNodeNormalRuptureStress", 1.0e99);
  m_tiedNodeShearRuptureStress = hdn->GetAttributeOrDefault<realT> ("tiedNodeShearRuptureStress", 1.0e99);
  m_tiedNodeTolerance = hdn->GetAttributeOrDefault<realT> ("tiedNodeTolerance", 1.0e-8);

}



void LagrangeExplicitDynamicsSolver::RegisterFields( PhysicalDomainT& domain )
{

  // register nodal fields
  domain.m_feNodeManager.AddKeyedDataField<FieldInfo::displacement>();
  domain.m_feNodeManager.AddKeyedDataField<FieldInfo::incrementalDisplacement>();
  domain.m_feNodeManager.AddKeyedDataField<FieldInfo::velocity>();
  domain.m_feNodeManager.AddKeyedDataField<FieldInfo::acceleration>();
  domain.m_feNodeManager.AddKeyedDataField<FieldInfo::force>();
  domain.m_feNodeManager.AddKeyedDataField<FieldInfo::mass>();

  domain.m_feNodeManager.AddKeyedDataField<FieldInfo::hgforce>();
  domain.m_feNodeManager.AddKeyedDataField<FieldInfo::contactForce>();

  domain.m_feNodeManager.AddKeylessDataField<realT>("work", true, false);

  domain.m_feNodeManager.AddKeylessDataField<realT>("nsigma_x", false, true);
  domain.m_feNodeManager.AddKeylessDataField<realT>("nsigma_y", false, true);
  domain.m_feNodeManager.AddKeylessDataField<realT>("nsigma_z", false, true);
  domain.m_feNodeManager.AddKeylessDataField<realT>("nsigma_xy", false, true);
  domain.m_feNodeManager.AddKeylessDataField<realT>("nsigma_yz", false, true);
  domain.m_feNodeManager.AddKeylessDataField<realT>("nsigma_zx", false, true);


  // register discrete element fields - all of the nodal fields + currentPosition + rotation
  domain.m_discreteElementManager.AddKeyedDataField<FieldInfo::displacement>();
  domain.m_discreteElementManager.AddKeyedDataField<FieldInfo::incrementalDisplacement>();
  domain.m_discreteElementManager.AddKeyedDataField<FieldInfo::velocity>();
  domain.m_discreteElementManager.AddKeyedDataField<FieldInfo::acceleration>();
  domain.m_discreteElementManager.AddKeyedDataField<FieldInfo::force>();
  domain.m_discreteElementManager.AddKeyedDataField<FieldInfo::mass>();
  domain.m_discreteElementManager.AddKeylessDataField<realT>("work", true, false);
  //rotation states
  domain.m_discreteElementManager.AddKeyedDataField<FieldInfo::rotationalAxisIncrement>();
  domain.m_discreteElementManager.AddKeyedDataField<FieldInfo::rotationalMagnitudeIncrement>();
  domain.m_discreteElementManager.AddKeyedDataField<FieldInfo::rotationalVelocity>();
  domain.m_discreteElementManager.AddKeyedDataField<FieldInfo::rotationalAcceleration>();
  domain.m_discreteElementManager.AddKeyedDataField<FieldInfo::moment>();



  // register external faces??? ... nothing to do right now

  // trying out ellipsoidal discrete elements
  domain.m_ellipsoidalDiscreteElementManager.AddKeyedDataField<FieldInfo::referencePosition>();
  domain.m_ellipsoidalDiscreteElementManager.AddKeyedDataField<FieldInfo::displacement>();
  domain.m_ellipsoidalDiscreteElementManager.AddKeyedDataField<FieldInfo::incrementalDisplacement>();
  domain.m_ellipsoidalDiscreteElementManager.AddKeyedDataField<FieldInfo::velocity>();
  domain.m_ellipsoidalDiscreteElementManager.AddKeyedDataField<FieldInfo::acceleration>();
  domain.m_ellipsoidalDiscreteElementManager.AddKeyedDataField<FieldInfo::force>();
  domain.m_ellipsoidalDiscreteElementManager.AddKeyedDataField<FieldInfo::mass>();
  domain.m_ellipsoidalDiscreteElementManager.AddKeyedDataField<FieldInfo::currentPosition>();
  domain.m_ellipsoidalDiscreteElementManager.AddKeylessDataField<realT>("work", true, false);

  //rotation states
  domain.m_ellipsoidalDiscreteElementManager.AddKeyedDataField<FieldInfo::rotationAxis>();
  domain.m_ellipsoidalDiscreteElementManager.AddKeyedDataField<FieldInfo::rotationMagnitude>();
  domain.m_ellipsoidalDiscreteElementManager.AddKeyedDataField<FieldInfo::rotationalAxisIncrement>();
  domain.m_ellipsoidalDiscreteElementManager.AddKeyedDataField<FieldInfo::rotationalMagnitudeIncrement>();
  domain.m_ellipsoidalDiscreteElementManager.AddKeyedDataField<FieldInfo::rotationalVelocity>();
  domain.m_ellipsoidalDiscreteElementManager.AddKeyedDataField<FieldInfo::rotationalAcceleration>();
  domain.m_ellipsoidalDiscreteElementManager.AddKeyedDataField<FieldInfo::moment>();
  domain.m_ellipsoidalDiscreteElementManager.AddKeyedDataField<FieldInfo::rotationalInertia>();

}


void LagrangeExplicitDynamicsSolver::Initialize(PhysicalDomainT& domain, SpatialPartition& partition )
{
  for( std::map< std::string, ElementRegionT >::iterator i=domain.m_feElementManager.m_ElementRegions.begin() ;
      i != domain.m_feElementManager.m_ElementRegions.end() ; ++i )
  {
    i->second.CalculateNodalMasses( domain.m_feNodeManager ) ;
  }

  if( m_tiedNodesFlag )
  {
    domain.m_feFaceManager.AddKeylessDataField<int>("numKCBCnodesInFace",true,true);
    BoundaryConditionFunctions::BuildKinematicConstraintBoundaryCondition( domain.m_feNodeManager,
                                                                           domain.m_feFaceManager,
                                                                           m_KinematicConstraintNodes,
                                                                           m_tiedNodeTolerance );
  }
}


void LagrangeExplicitDynamicsSolver::InitializeCommunications( PartitionBase& partition )
{

  m_syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back(Field<FieldInfo::acceleration>::Name());
  m_syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back(Field<FieldInfo::velocity>::Name());
  m_syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back(Field<FieldInfo::displacement>::Name());

  m_syncedFields[PhysicalDomainT::DiscreteElementManager].push_back(Field<FieldInfo::acceleration>::Name());
  m_syncedFields[PhysicalDomainT::DiscreteElementManager].push_back(Field<FieldInfo::velocity>::Name());
  m_syncedFields[PhysicalDomainT::DiscreteElementManager].push_back(Field<FieldInfo::displacement>::Name());
  m_syncedFields[PhysicalDomainT::DiscreteElementManager].push_back(Field<FieldInfo::rotationalAcceleration>::Name());
  m_syncedFields[PhysicalDomainT::DiscreteElementManager].push_back(Field<FieldInfo::rotationalVelocity>::Name());
  m_syncedFields[PhysicalDomainT::DiscreteElementManager].push_back(Field<FieldInfo::rotationMagnitude>::Name());
  m_syncedFields[PhysicalDomainT::DiscreteElementManager].push_back(Field<FieldInfo::rotationAxis>::Name());

  m_syncedFields[PhysicalDomainT::FiniteElementFaceManager].push_back("cohesiveTraction");
  m_syncedFields[PhysicalDomainT::FiniteElementFaceManager].push_back("ruptureState");

  m_syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back("nsigma_x");
  m_syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back("nsigma_y");
  m_syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back("nsigma_z");
  m_syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back("nsigma_xy");
  m_syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back("nsigma_yz");
  m_syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back("nsigma_zx");


  partition.SetBufferSizes(m_syncedFields, CommRegistry::lagrangeSolver01);
}


void LagrangeExplicitDynamicsSolver::ApplyForcesFromContact(PhysicalDomainT& domain,
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

    //update nodal positions, velocities, and accelerations before updating face geometry
    //also, reset rotational and translational accelerations as well as forces and moments
    domain.m_discreteElementManager.UpdateNodalStatesZeroForcesAndAccelerations();

    //update face geometry and sort faces if necessary
    const bool resort = domain.m_externalFaces.RecalculateNeighborList(domain.m_feNodeManager,
                                                                       domain.m_discreteElementSurfaceNodes,
                                                                       domain.m_discreteElementManager,
                                                                       dt);

    //if a resort has been triggered, then you also need to update the contact manager
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

    //for parallel: do NOT need to synchronize nodal force fields across processes for DE nodes
    //before moving to centroid, since each process will do that calculation (redundantly) itself
    // ... there's therefore no need for explicit synchrony
    //as long as nodal states remain synchronous, everything will be fine!
#ifndef DEEFC
    domain.m_discreteElementManager.ApplyDiscreteElementContactForces(45, 0.2);
#else
    domain.m_discreteElementManager.ApplyNodalForces();
#endif
  }


  //-----------------------------
  //ellipsoidal discrete elements
  //-----------------------------
  {
    //update nodal positions, velocities, and accelerations before updating face geometry
    //also, reset rotational and translational accelerations as well as forces and moments
    domain.m_ellipsoidalDiscreteElementManager.UpdateNodalStatesZeroForcesAndAccelerations();

    //update ellipsoidal discrete elements
    bool resort = domain.m_ellipsoidalDiscreteElementManager.RecalculateNeighborList(dt);

    //if a resort has been triggered, then you also need to update the contact manager
    if(resort)
      domain.m_ellipsoidalContactManager.Update(domain.m_ellipsoidalDiscreteElementManager.m_neighborList);

    domain.m_ellipsoidalDiscreteElementManager.UpdateAndApplyContactStresses(timeStep, dt, domain.m_ellipsoidalContactManager);
  }


}

double LagrangeExplicitDynamicsSolver::TimeStep( const realT& time,
                                               const realT& dt,
                                               const int cycleNumber,
                                               PhysicalDomainT& domain,
                                               const array<string>& namesOfSolverRegions,
                                               SpatialPartition& partition,
                                               FractunatorBase* const fractunator )
{
  using namespace BoundaryConditionFunctions;

  realT dt_return = dt;

  m_stabledt.m_maxdt = std::numeric_limits<double>::max();

  array<R1Tensor>& hgforce = domain.m_feNodeManager.GetFieldData<FieldInfo::hgforce> ();
  hgforce = 0.0;

  /////////////////////////////////
  // HACK for testing dry pulling
  if( fractunator!=NULL )
  {
    if( cycleNumber%fractunator->m_checkInterval==0 )
    {
      fractunator->SeparationDriver( domain.m_feNodeManager,
                                     domain.m_feEdgeManager,
                                     domain.m_feFaceManager,
                                     domain.m_externalFaces,
                                     domain.m_feElementManager,
                                     partition, false, time );
    }
  }
  /////////////////////////////////

  // update nodes
  LagrangeExplicitDynamicsFunctions::LinearPointUpdatePart1( domain.m_feNodeManager, time, dt);

  // update discrete elements
  LagrangeExplicitDynamicsFunctions::LinearPointUpdatePart1( domain.m_discreteElementManager, time, dt, false );
  LagrangeExplicitDynamicsFunctions::LinearPointUpdatePart1( domain.m_ellipsoidalDiscreteElementManager, time, dt, false);

  // update discrete element rotations
  LagrangeExplicitDynamicsFunctions::RotationalPointUpdatePart1( domain.m_discreteElementManager, time, dt );
  LagrangeExplicitDynamicsFunctions::RotationalPointUpdatePart1b( domain.m_discreteElementManager);
  LagrangeExplicitDynamicsFunctions::RotationalPointUpdatePart1( domain.m_ellipsoidalDiscreteElementManager, time, dt );

  for( array<string>::const_iterator regionName = namesOfSolverRegions.begin() ;
      regionName != namesOfSolverRegions.end() ; ++regionName )
  {
    //this conditional supports DE, since an element region in DE may not exist in FE; this case should not
    //throw an exception, which occurs in the absence of the following block
    std::map<std::string, ElementRegionT>::iterator iter = domain.m_feElementManager.m_ElementRegions.find(*regionName);
    if(iter == domain.m_feElementManager.m_ElementRegions.end())
      continue;

    ElementRegionT& elementRegion = iter->second; // stlMapLookup( domain.m_elementManager.m_ElementRegions, *regionName );
    elementRegion.CalculateVelocityGradients(domain.m_feNodeManager);
    elementRegion.MaterialUpdate(dt);
  }

  for( array<string>::const_iterator regionName = namesOfSolverRegions.begin() ;
      regionName != namesOfSolverRegions.end() ; ++regionName )
  {
    //this conditional supports DE, since an element region in DE may not exist in FE; this case should not
    //throw an exception, which occurs in the absence of the following block
    std::map<std::string, ElementRegionT>::iterator iter = domain.m_feElementManager.m_ElementRegions.find(*regionName);
    if(iter == domain.m_feElementManager.m_ElementRegions.end())
      continue;

    ElementRegionT& elementRegion = iter->second; // domain.m_elementManager.m_ElementRegions[*regionName];

    elementRegion.CalculateNodalForces(domain.m_feNodeManager, m_stabledt, dt);
  }

  BoundaryConditionFunctions::ApplyTractionBoundaryCondition( domain,time );

  //  ApplyBoundaryCondition(&ApplyTractionBoundaryCondition,
  //                                domain, domain.m_faceManager, "Traction" );

  // -- CONTACT --
  if (domain.m_externalFaces.m_contactActive)
    ApplyForcesFromContact(domain, this->m_stabledt, dt);
  domain.m_ellipsoidalDiscreteElementManager.UpdateCylindricalBoundary(m_stabledt, dt, time);

  if( m_tiedNodesFlag )
  {
    BoundaryConditionFunctions::ApplyKinematicConstraintBoundaryCondition( domain.m_feFaceManager,
                                                                           domain.m_feNodeManager,
                                                                           this->m_KinematicConstraintNodes,
                                                                           m_tiedNodeNormalRuptureStress,
                                                                           m_tiedNodeShearRuptureStress );
  }


  if( domain.m_feNodeManager.HasField<R1Tensor>("cohesiveForce") &&
      domain.m_feFaceManager.HasField<R1Tensor>("cohesiveTraction") &&
      domain.m_feFaceManager.HasField<int>("ruptureState") )
  {
    // apply cohesive forces
    const OrderedVariableOneToManyRelation& childFaceIndex = domain.m_feFaceManager.GetVariableOneToManyMap( "childIndices" );
    array<R1Tensor>& nodalForce = domain.m_feNodeManager.GetFieldData<FieldInfo::force>();
    array<R1Tensor>& cohesiveForce = domain.m_feNodeManager.GetFieldData<R1Tensor>("cohesiveForce");

    const array<R1Tensor>& cohesiveTraction = domain.m_feFaceManager.GetFieldData<R1Tensor>("cohesiveTraction");
    const array<integer>& ruptureState = domain.m_feFaceManager.GetFieldData<int>("ruptureState");

    cohesiveForce = 0.0;

    for( localIndex kf=0 ; kf<domain.m_feFaceManager.DataLengths() ; ++kf )
    {
      if( ruptureState[kf] == 2 && !(childFaceIndex[kf].empty()) )
      {


        const localIndex faceIndex[2] = { kf, childFaceIndex[kf][0] };


        const R1Tensor N[2] = { domain.m_feFaceManager.FaceNormal( domain.m_feNodeManager, faceIndex[0] ),
                                domain.m_feFaceManager.FaceNormal( domain.m_feNodeManager, faceIndex[1] )};

        R1Tensor Nbar = N[0];
        Nbar -= N[1];
        Nbar.Normalize();

        for( int side=0 ; side<2 ; ++side )
        {
          const localIndex faceID = faceIndex[side];
          const realT area = domain.m_feFaceManager.SurfaceArea( domain.m_feNodeManager, faceID );

          R1Tensor cForce;
          int direction = -1;
          if(side==0)
            direction = 1;

          cForce = cohesiveTraction[kf];
          cForce *= direction * area / domain.m_feFaceManager.m_toNodesRelation[faceID].size();


          for( lArray1d::const_iterator nodeID=domain.m_feFaceManager.m_toNodesRelation[faceID].begin() ;
              nodeID!=domain.m_feFaceManager.m_toNodesRelation[faceID].end() ; ++nodeID )
          {
            nodalForce[*nodeID] += cForce;
            cohesiveForce[*nodeID] += cForce;
          }

        }
      }
    }
  }


#ifdef SRC_INTERNAL
  //FIXME: replace the following call with domain.m_externalFaces.GeodynCouplingParallel( domain.m_feNodeManager );
  domain.m_externalFaces.GeodynCoupling( domain.m_feNodeManager );
#endif
  LagrangeExplicitDynamicsFunctions::LinearPointUpdatePart2( domain.m_feNodeManager, time, dt, m_dampingM );
  LagrangeExplicitDynamicsFunctions::LinearPointUpdatePart2( domain.m_discreteElementManager, time, dt );
  LagrangeExplicitDynamicsFunctions::RotationalPointUpdatePart2( domain.m_discreteElementManager, time, dt );

  LagrangeExplicitDynamicsFunctions::LinearPointUpdatePart2( domain.m_ellipsoidalDiscreteElementManager, time, dt );
  LagrangeExplicitDynamicsFunctions::RotationalPointUpdatePart2( domain.m_ellipsoidalDiscreteElementManager, time, dt );


  if(1)
  {

    partition.SynchronizeFields( m_syncedFields, CommRegistry::lagrangeSolver01 );
  }

  m_stabledt.m_maxdt *= this->m_courant;
  return dt_return;
}

void LagrangeExplicitDynamicsSolver::PostProcess (PhysicalDomainT& domain,
                                                  SpatialPartition& partition,
                                                  const array<string>& namesOfSolverRegions)
{
  array<real64>& sigma_x = domain.m_feNodeManager.GetFieldData<realT>("nsigma_x");
  array<real64>& sigma_y = domain.m_feNodeManager.GetFieldData<realT>("nsigma_y");
  array<real64>& sigma_z = domain.m_feNodeManager.GetFieldData<realT>("nsigma_z");
  array<real64>& sigma_xy = domain.m_feNodeManager.GetFieldData<realT>("nsigma_xy");
  array<real64>& sigma_yz = domain.m_feNodeManager.GetFieldData<realT>("nsigma_yz");
  array<real64>& sigma_zx = domain.m_feNodeManager.GetFieldData<realT>("nsigma_zx");

  sigma_x = 0.0;
  sigma_y = 0.0;
  sigma_z = 0.0;
  sigma_xy = 0.0;
  sigma_yz = 0.0;
  sigma_zx = 0.0;

  for( array<string>::const_iterator regionName = namesOfSolverRegions.begin() ;
      regionName != namesOfSolverRegions.end() ; ++regionName )
  {
    std::map<std::string, ElementRegionT>::iterator iter = domain.m_feElementManager.m_ElementRegions.find(*regionName);
    if(iter == domain.m_feElementManager.m_ElementRegions.end())
      continue;

    ElementRegionT& elementRegion = iter->second;


    for( localIndex k=0 ; k<elementRegion.m_numElems ; ++k )
    {
      R2SymTensor s;
      realT pressure = 0.0;
      realT ex, ey, ez, exy, eyz, ezx;
      for( localIndex a=0 ; a<elementRegion.m_numIntegrationPointsPerElem ; ++a )
      {
        const MaterialBaseStateData& state = *(elementRegion.m_mat->StateData(k,a));
        s += state.devStress;
        pressure += state.pressure;
      }
      s /= elementRegion.m_numIntegrationPointsPerElem;
      pressure /= elementRegion.m_numIntegrationPointsPerElem;
      ex =  s(0,0) + pressure;
      ey =  s(1,1) + pressure;
      ez =  s(2,2) + pressure;
      exy =  s(0,1);
      eyz =  s(1,2);
      ezx =  s(0,2);

      for (localIndex j = 0; j < elementRegion.m_toNodesRelation.Dimension(1); ++j)
      {
        localIndex ind = elementRegion.m_toNodesRelation[k][j];
        sigma_x[ind] +=ex;
        sigma_y[ind] +=ey;
        sigma_z[ind] +=ez;
        sigma_xy[ind] +=exy;
        sigma_yz[ind] +=eyz;
        sigma_zx[ind] +=ezx;
      }
    }
  }

  for (localIndex i=0; i<domain.m_feNodeManager.DataLengths(); ++i)
  {
    if (domain.m_feNodeManager.m_toElementsRelation[i].size() > 0)
    {
      sigma_x[i] /= domain.m_feNodeManager.m_toElementsRelation[i].size();
      sigma_y[i] /= domain.m_feNodeManager.m_toElementsRelation[i].size();
      sigma_z[i] /= domain.m_feNodeManager.m_toElementsRelation[i].size();
      sigma_xy[i] /= domain.m_feNodeManager.m_toElementsRelation[i].size();
      sigma_yz[i] /= domain.m_feNodeManager.m_toElementsRelation[i].size();
      sigma_zx[i] /= domain.m_feNodeManager.m_toElementsRelation[i].size();
    }
  }

  {
    std::map<PhysicalDomainT::ObjectDataStructureKeys, array<string>> syncedFields;
    syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back("nsigma_x");
    syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back("nsigma_y");
    syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back("nsigma_z");
    syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back("nsigma_xy");
    syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back("nsigma_yz");
    syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back("nsigma_zx");

    partition.SynchronizeFields( syncedFields, CommRegistry::lagrangeSolver01);
  }



}

void LagrangeExplicitDynamicsSolver::ApplyGapDamping( NodeManager& nodeManager,
                                                      const FaceManagerT& faceManager,
                                                      const realT dt )
{
  if( dt>0.0 )
  {
    const array<R1Tensor>& velocity = nodeManager.GetFieldData<FieldInfo::velocity>();
    array<R1Tensor>& force = nodeManager.GetFieldData<FieldInfo::force>();
    const array<realT>& mass = nodeManager.GetFieldData<FieldInfo::mass>();

    const OrderedVariableOneToManyRelation& childFaces = faceManager.GetVariableOneToManyMap("childIndices");

    R1Tensor faceVelocity[2];
    R1Tensor gapVelocity;
    R1Tensor dampForce;
    for( localIndex k=0 ; k<faceManager.DataLengths() ; ++k )
    {

      if( !(childFaces[k].empty()) )
      {
        const localIndex faceIndex[2] = { childFaces[k][0],
                                          childFaces[k][1] } ;

        for( int i=0 ; i<2 ; ++i )
        {
          faceVelocity[i] = 0;
        }

        realT faceMass = 0.0;
        localIndex numNodesInFace = faceManager.m_toNodesRelation[k].size();

        const lArray1d& nodeIndices0 = faceManager.m_toNodesRelation[faceIndex[0]] ;
        const lArray1d& nodeIndices1 = faceManager.m_toNodesRelation[faceIndex[1]] ;

        for( localIndex a=0 ; a<numNodesInFace ; ++a )
        {
          faceVelocity[0] += velocity[nodeIndices0[a]];
          faceVelocity[1] += velocity[nodeIndices1[a]];

          faceMass += mass[nodeIndices0[a]];
          faceMass += mass[nodeIndices1[a]];

        }
        faceVelocity[0] /= numNodesInFace;
        faceVelocity[1] /= numNodesInFace;
        faceMass /= 2*numNodesInFace;

        gapVelocity  = faceVelocity[1];
        gapVelocity -= faceVelocity[0];

        dampForce = gapVelocity;
        dampForce *= - this->m_gapdamping * faceMass / dt * 0.25;

        for( localIndex a=0 ; a<numNodesInFace ; ++a )
        {
          force[nodeIndices0[a]] -= dampForce;
          force[nodeIndices1[a]] += dampForce;
        }


      }
    }
  }
}



void LagrangeExplicitDynamicsSolver::SetMaxStableTimeStep( const realT& time,
                                                           PhysicalDomainT& domain,
                                                           const array<string>& namesOfSolverRegions,
                                                           SpatialPartition& partition  )
{
  array<R1Tensor>& incrementalDisplacement = domain.m_feNodeManager.GetFieldData<FieldInfo::incrementalDisplacement> ();
  incrementalDisplacement = 0.0;

  for( array<string>::const_iterator regionName = namesOfSolverRegions.begin() ;
      regionName != namesOfSolverRegions.end() ; ++regionName )
  {
    //this conditional supports DE, since an element region in DE may not exist in FE; this case should not
    //throw an exception, which occurs in the absence of the following block
    std::map<std::string, ElementRegionT>::iterator iter = domain.m_feElementManager.m_ElementRegions.find(*regionName);
    if(iter == domain.m_feElementManager.m_ElementRegions.end())
      continue;

    ElementRegionT& elementRegion = iter->second; // stlMapLookup( domain.m_elementManager.m_ElementRegions, *regionName );
    elementRegion.CalculateVelocityGradients(domain.m_feNodeManager);
    elementRegion.CalculateNodalForces(domain.m_feNodeManager, m_stabledt, 0.0);
  }

  m_stabledt.m_maxdt *= this->m_courant;

}



void LagrangeExplicitDynamicsSolver::WriteSiloDerived( SiloFile& siloFile ) const
{

}

void LagrangeExplicitDynamicsSolver::ReadSiloDerived( const SiloFile& siloFile )
{

}


namespace LagrangeExplicitDynamicsFunctions
{
  // *********************************************************************************************************************
  /**
   * @author R. Settgast
   * @param dt time increment of current time step
   *
   * This function is an aggregate operation that multiplies that performs the following
   * operations for each node:
   *
   * push nodal velocity to half step
   * \f[ v_{a}^{n+1/2} = v_{a}^{n} + a_a^{n} * (dt/2) \f]
   *
   * calculate incremental displacement over the step
   * \f[ uhat_{a}^{n+1/2} = v_{a}^{n+1/2} * dt \f]
   *
   * calculate total displacement at end of step
   * \f[ u_{a}^{n+1} = u_{a}^{n} + uhat_{a}^{n+1/2} \f]
   *
   * zero nodal forces
   * \f[ f_{a}^{n+1} = 0 \f]
   *
   *
   */
  void LinearPointUpdatePart1( ObjectDataStructureBaseT& objectManager,
                               const realT& time,
                               const realT& dt,
                               const bool clearForces)
  {
    if( objectManager.DataLengths() > 0 )
    {
      array<R1Tensor>& velocity = objectManager.GetFieldData<FieldInfo::velocity> ();
      const array<R1Tensor>& acceleration = objectManager.GetFieldData<FieldInfo::acceleration> ();
      array<R1Tensor>& incrementalDisplacement = objectManager.GetFieldData<FieldInfo::incrementalDisplacement> ();
      array<R1Tensor>& displacement = objectManager.GetFieldData<FieldInfo::displacement> ();
      array<R1Tensor>& force = objectManager.GetFieldData<FieldInfo::force> ();

      array<real64>& work = objectManager.GetFieldData<realT>("work");


      const realT dtdiv2 = 0.5 * dt;

      for (localIndex a = 0; a < objectManager.DataLengths(); ++a)
      {
        // push nodal velocity forward to the half-step
        velocity[a].plus_cA(dtdiv2, acceleration[a]);
      }

      BoundaryConditionFunctions::ApplyDirichletBoundaryCondition<R1Tensor>(objectManager,
                                                                            Field<FieldInfo::velocity>::Name(),time+0.5*dt);


      BoundaryConditionFunctions::ApplyRigidWallBoundaryCondition( objectManager, dt );

      for (localIndex a = 0; a < objectManager.DataLengths(); ++a)
      {
        // calculate incremental displacements over the step
        incrementalDisplacement[a].cA(dt, velocity[a]);

        work[a] += 0.5 * Dot(force[a],incrementalDisplacement[a]);

        // add incremental displacements to total displacements to bring to end of step
        displacement[a] += incrementalDisplacement[a];

      }


      /*
      const array<lArray1d>& childIndices = objectManager.GetVariableOneToManyMap( "childIndices" );
      for (localIndex a = 0; a < objectManager.DataLengths(); ++a)
      {
        if( !(childIndices[a].empty()) )
        {
          incrementalDisplacement[a] = 0.5*( incrementalDisplacement[childIndices[a][0]] + incrementalDisplacement[childIndices[a][1]] );
          displacement[a] = 0.5*( displacement[childIndices[a][0]] + displacement[childIndices[a][1]] );
        }
      }*/

      // set forces to zero
      if(clearForces)
      {
        force = 0.0;
      }
    }
  }

  // *********************************************************************************************************************
  /**
   * @author R. Settgast
   * @param dt time increment of current time step
   *
   * This function is an aggregate operation that multiplies that performs the following
   * operations for each node:
   *
   * calculate nodal acceleration at end of step
   * \f[ a_{a}^{n+1} = f_{a}^{n+1} / m_a^{n+1} \f]
   *
   * push nodal velocity to end of step
   * \f[ v_{a}^{n+1} = v_{a}^{n+1/2} + a_a^{n+1} * (dt/2) \f]
   *
   */
  void LinearPointUpdatePart2( ObjectDataStructureBaseT& objectManager,
                               const realT& time,
                               const realT& dt,
                               const realT& damping )
  {
    if( objectManager.DataLengths() > 0 )
    {
      array<R1Tensor>& velocity = objectManager.GetFieldData<FieldInfo::velocity> ();
      array<R1Tensor>& acceleration = objectManager.GetFieldData<FieldInfo::acceleration> ();
      array<R1Tensor>& force = objectManager.GetFieldData<FieldInfo::force> ();
      array<realT>& mass = objectManager.GetFieldData<FieldInfo::mass> ();
      array<real64>& work = objectManager.GetFieldData<realT>("work");

      const array<R1Tensor>& incrementalDisplacement = objectManager.GetFieldData<FieldInfo::incrementalDisplacement> ();
      const array<int>* isDetachedFromSolidMesh = objectManager.GetFieldDataPointer<int> ("isDetachedFromSolidMesh");

      const realT dtdiv2 = 0.5 * dt;
      realT dampedMass;
      R1Tensor dampedForce;
      for (localIndex a = 0; a < objectManager.DataLengths(); ++a)
      {

        if ( isDetachedFromSolidMesh == NULL )
        {
          // calculate acceleration
          if( damping > 0.0 )
          {
            dampedMass = mass[a] * ( 1 + 0.5 * dt * damping );
            dampedForce.cA( -damping * mass[a], velocity[a] );
            dampedForce += force[a];
            acceleration[a].Adivc(dampedMass, dampedForce );
          }
          else
          {
            acceleration[a].Adivc(mass[a], force[a]);
          }

          work[a] += 0.5 * Dot(force[a],incrementalDisplacement[a]);
        }
        else
        {
          if ( (*isDetachedFromSolidMesh)[a] == 0 )
          {
            if( damping > 0.0 )
            {
              dampedMass = mass[a] * ( 1 + 0.5 * dt * damping );
              dampedForce.cA( -damping * mass[a], velocity[a] );
              dampedForce += force[a];
              acceleration[a].Adivc(dampedMass, dampedForce );
            }
            else
            {
              acceleration[a].Adivc(mass[a], force[a]);
            }

            work[a] += 0.5 * Dot(force[a],incrementalDisplacement[a]);

          }
          else
          {
            acceleration[a] = 0;
          }
        }
      }

      BoundaryConditionFunctions::ApplyDirichletBoundaryCondition<R1Tensor>(objectManager,
                                                                            Field<FieldInfo::acceleration>::Name(),time);

      for (localIndex a = 0; a < objectManager.DataLengths(); ++a)
      {
        // push velocity forward to end of step
        velocity[a].plus_cA(dtdiv2, acceleration[a]);

        if (mass[a] <= 0) velocity[a] = 0;
      }

      BoundaryConditionFunctions::ApplyDirichletBoundaryCondition<R1Tensor>(objectManager,
                                                                            Field<FieldInfo::velocity>::Name(),time);//+dt);
    }
  }

  /**
   * @brief First part of the rotation update (also updates current positions of de's and nodes)
   * @author Scott Johnson
   *
   * This function is an aggregate operation that multiplies that performs the following
   * operations for each node (after Omelyan, 1998). Note that this interpretation of the
   * quaternion can be parameterized using the rotation, \f[\theta\f], about a unit axis, n:
   * \f[ q[0] = cos(0.5 * \theta) \f]
   * \f[ q[1..3] = sin(0.5 * \theta) * n\f]
   *
   * calculate quaternion time derivative at the step
   * \f[ \dot{q}_{a}^{n} = \frac{1}{2}Q\left(q_{a}^{n}\right)\omega_{a}^{n} \f]
   * \f[ Q\left(q_{a}^{n}\right)=\left[\begin{array}{ccc}-q_{a}^{n}[1] & -q_{a}^{n}[2] & -q_{a}^{n}[3] \\ q_{a}^{n}[0] & -q_{a}^{n}[3] & q_{a}^{n}[2] \\ q_{a}^{n}[3] & q_{a}^{n}[0] & -q_{a}^{n}[1] \\ -q_{a}^{n}[2] & q_{a}^{n}[1] & q_{a}^{n}[0]\end{array}\right]\f]
   *
   * push quaternion to the half step
   * \f[ q_{a}^{n+1/2} = q_{a}^{n} + \dot{q}_{a}^{n} * (dt/2) \f]
   *
   * push body frame nodal rotational velocity to half step
   * \f[ \omega_{a}^{n+1/2} = \omega_{a}^{n} + \alpha_a^{n} * (dt/2) \f]
   *
   * calculate quaternion time derivative at the half step
   * \f[ \dot{q}_{a}^{n+1/2} = \frac{1}{2}Q\left(q_{a}^{n+1/2}\right)\omega_{a}^{n+1/2} \f]
   * \f[ Q\left(q_{a}^{n+1/2}\right)=\left[\begin{array}{ccc}-q_{a}^{n+1/2}[1] & -q_{a}^{n+1/2}[2] & -q_{a}^{n+1/2}[3] \\ q_{a}^{n+1/2}[0] & -q_{a}^{n+1/2}[3] & q_{a}^{n+1/2}[2] \\ q_{a}^{n+1/2}[3] & q_{a}^{n+1/2}[0] & -q_{a}^{n+1/2}[1] \\ -q_{a}^{n+1/2}[2] & q_{a}^{n+1/2}[1] & q_{a}^{n+1/2}[0]\end{array}\right]\f]
   *
   * calculate incremental quaternion over the step
   * \f[ {dq}_{a}^{n+1/2} = \dot{q}_{a}^{n+1/2} * dt \f]
   *
   * calculate quaternion at end of step
   * \f[ q_{a}^{n+1} = q_{a}^{n} + {dq}_{a}^{n+1/2} \f]
   *
   * renormalize quaternion
   * \f[ |q| = 1 \f]
   *
   * @param discreteElementManager manager of the discrete elements to update
   * @param time current simulation time
   * @param dt current timestep length
   */
  void RotationalPointUpdatePart1(DiscreteElementManagerBaseT& discreteElementManager,
                                  const realT& time, const realT& dt)
  {

    if( discreteElementManager.DataLengths() > 0 )
    {

      array<R1Tensor>& rotationalVelocity        = discreteElementManager.GetFieldData<FieldInfo::rotationalVelocity> ();
      array<R1Tensor>& rotationalAcceleration    = discreteElementManager.GetFieldData<FieldInfo::rotationalAcceleration> ();
      array<R1Tensor>& rotationalAxisIncrement   = discreteElementManager.GetFieldData<FieldInfo::rotationalAxisIncrement> ();
      array<realT>& rotationalMagnitudeIncrement = discreteElementManager.GetFieldData<FieldInfo::rotationalMagnitudeIncrement> ();
      array<R1Tensor>& rotationAxis              = discreteElementManager.GetFieldData<FieldInfo::rotationAxis> ();
      array<realT>& rotationMagnitude            = discreteElementManager.GetFieldData<FieldInfo::rotationMagnitude> ();

      const realT dtdiv2 = 0.5 * dt;

      //note: all quantities are in body frame
      for (localIndex a = 0; a < discreteElementManager.DataLengths(); ++a)
      {
        //push the quaternion to the half step
        realT qh[] = {rotationMagnitude[a], rotationAxis[a][0], rotationAxis[a][1], rotationAxis[a][2]};
        realT dqdt[nsdof + 1];
        discreteElementManager.Calculate_dqdt(a, qh, dqdt);
        qh[0] = rotationMagnitude[a] + dtdiv2 * dqdt[0];
        for (unsigned int i = 0; i < nsdof; i++)
          qh[i + 1] = rotationAxis[a][i] + dtdiv2 * dqdt[i + 1];

        // push rotational velocity forward to the half-step
        rotationalVelocity[a].plus_cA(dtdiv2, rotationalAcceleration[a]);

        // calculate dq based on the half step
        discreteElementManager.Calculate_dqdt(a, qh, dqdt);
        rotationalMagnitudeIncrement[a] = dqdt[0] * dt;
        rotationalAxisIncrement[a][0] = dqdt[1] * dt;
        rotationalAxisIncrement[a][1] = dqdt[2] * dt;
        rotationalAxisIncrement[a][2] = dqdt[3] * dt;

        // add incremental rotational displacements to total displacements to bring to end of step
        rotationMagnitude[a] += rotationalMagnitudeIncrement[a];
        rotationAxis[a] += rotationalAxisIncrement[a];

        // renormalize quaternion if necessary
        {
          realT fct = rotationMagnitude[a] * rotationMagnitude[a];
          for (unsigned int i = 0; i < nsdof; i++)
            fct += rotationAxis[a][i] * rotationAxis[a][i];
          if (fct > 0 && !isEqual(fct,1.0) )
          {
            fct = 1. / sqrt(fct);
            rotationMagnitude[a] *= fct;
            for (unsigned int i = 0; i < nsdof; i++)
              rotationAxis[a][i] *= fct;
          }
        }
      }

    }
  }

  void RotationalPointUpdatePart1b(DiscreteElementManagerT& discreteElementManager)
  {
    if( discreteElementManager.DataLengths() > 0 )
    {
      //discrete element positions
      array<R1Tensor>& deCurrentPosition         = discreteElementManager.GetFieldData<FieldInfo::currentPosition> ();
      array<R1Tensor>& deReferencePosition       = discreteElementManager.GetFieldData<FieldInfo::referencePosition> ();
      array<R1Tensor>& deDisplacement            = discreteElementManager.GetFieldData<FieldInfo::displacement> ();

      //nodal positions
      array<R1Tensor>& nodeRelativePosition      = discreteElementManager.m_nodeManager->GetFieldData<FieldInfo::relativePosition> ();
      array<R1Tensor>& nodeCurrentPosition       = discreteElementManager.m_nodeManager->GetFieldData<FieldInfo::currentPosition> ();

      //update nodal and discrete element current positions
      //note: all quantities are in body frame
      for (localIndex a = 0; a < discreteElementManager.DataLengths(); ++a)
      {
        R2Tensor rotation;
        discreteElementManager.RotationTensor(a, rotation);

        //--APPLY THE ROTATION TENSOR--
        deCurrentPosition[a] = deReferencePosition[a];
        deCurrentPosition[a] += deDisplacement[a];
        for(localIndex b=0; b<discreteElementManager.m_discreteElementToExternalNodesMap[a].size(); b++)
        {
          int nn = discreteElementManager.m_discreteElementToExternalNodesMap[a][b];
          //****local to global direction transform (see DiscreteElementManagerBaseT.h)
          //    R2Tensor Rt;
          //    RotationTensorTranspose(a, Rt);
          //    global.AijBj(Rt, local);
          nodeCurrentPosition[nn].AijBi(rotation, nodeRelativePosition[nn]);
          nodeCurrentPosition[nn] += deCurrentPosition[a];
        }
      }
    }
  }

  // *********************************************************************************************************************
  /**
   * @author Scott Johnson
   * @param objectManager manager of objects to operate on
   * @param time current simulation time
   * @param dt time increment of current time step
   *
   * This function is an aggregate operation that multiplies that performs the following
   * operations for each node:
   *
   * calculate body frame nodal rotational acceleration at end of step
   * \f[ \alpha_{a}^{n+1} = \tau_{a}^{n+1} / I_a^{n+1} \f]
   *
   * push body frame nodal rotational velocity to end of step
   * \f[ \omega_{a}^{n+1} = \omega_{a}^{n+1/2} + \alpha_a^{n+1} * (dt/2) \f]
   *
   * zero body frame moments
   * \f[ m_{a}^{n+1} = 0 \f]
   *
   */
  void RotationalPointUpdatePart2(ObjectDataStructureBaseT& objectManager,
                                  const realT& time, const realT& dt)
  {
    if( objectManager.DataLengths() > 0 )
    {
      array<R1Tensor>& rotationalVelocity = objectManager.GetFieldData< FieldInfo::rotationalVelocity> ();
      array<R1Tensor>& rotationalAcceleration = objectManager.GetFieldData< FieldInfo::rotationalAcceleration> ();
      array<R1Tensor>& moment = objectManager.GetFieldData<FieldInfo::moment> ();
      array<R1Tensor>& rotationalInertia = objectManager.GetFieldData<FieldInfo::rotationalInertia> ();

      const realT dtdiv2 = 0.5 * dt;
      for (localIndex a = 0; a < objectManager.DataLengths(); ++a)
      {
        // calculate rotational acceleration
        rotationalAcceleration[a] = moment[a];
        rotationalAcceleration[a] /= rotationalInertia[a];

        // push rotational velocity forward to end of step
        rotationalVelocity[a].plus_cA(dtdiv2, rotationalAcceleration[a]);
      }

      // set moments to zero
      moment = 0.0;
    }
  }
}





/// Register solver in the solver factory
//REGISTER_SOLVER( LagrangeExplicitDynamicsSolver )
