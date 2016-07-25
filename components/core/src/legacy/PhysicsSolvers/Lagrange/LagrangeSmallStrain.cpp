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
 * @file LagrangeSmallStrain.cpp
 * @author Randolph Settgast
 * @date created on Sep 13, 2010
 */

#include "LagrangeSmallStrain.h"
#include "PhysicsSolvers/SolverFactory.h"
#include "PhysicsSolvers/PhysicsSolverStrings.h"

#include "Utilities/Utilities.h"
#include "Utilities/Kinematics.h"
#include "ObjectManagers/ProblemManagerT.h"
#include "ObjectManagers/ElementManagerT.h"
#include "ElementLibrary/FiniteElement.h"
#include "ElementLibrary/LagrangeBasis.h"
#include "ElementLibrary/GaussQuadrature.h"

#include "ObjectManagers/TableManager.h"
#include "BoundaryConditions/BoundaryConditions.h"
#include "BoundaryConditions/ApplyBoundaryConditions.h"




#include "Epetra_Map.h"
#include "Epetra_FECrsGraph.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_LinearProblem.h"

#include "EpetraExt_RowMatrixOut.h"
#include "Teuchos_RCP.hpp"
#include "AztecOO.h"
#include "ml_MultiLevelPreconditioner.h"



LagrangeSmallStrain::LagrangeSmallStrain(  const std::string& name,
                                                ProblemManagerT* const pm ):
LagrangeSolverBase(name,pm),
m_nonContactModulus(0.0),
m_recordIncrementalDisplacement(false)
{}


LagrangeSmallStrain::~LagrangeSmallStrain()
{
  // TODO Auto-generated destructor stub
}



void LagrangeSmallStrain::ReadXML( TICPP::HierarchicalDataNode* const hdn )
{
  LagrangeSolverBase::ReadXML(hdn);
  m_nonContactModulus = hdn->GetAttributeOrDefault("NonContactModulus","2.2 GPa"); // fixme incorrect units
  m_recordIncrementalDisplacement = hdn->GetAttributeOrDefault("RecordIncrementalDisplacement",false); 
  m_doApertureUpdate = hdn->GetAttributeOrDefault("UpdateAperture",false); 
  m_hasReferenceStress = hdn->GetAttributeOrDefault("hasReferenceStress",false);
  // If a model has reference stress defined, this is the stress state in its reference configuration.
  // In the end, the stress in VisIt is the total stress and the displacement is the displacement relative to the reference configuration.
  
  ++m_instances;
}




void LagrangeSmallStrain::RegisterFields( PhysicalDomainT& domain )
{
  LagrangeSolverBase::RegisterFields(domain);
  // register nodal fields
  domain.m_feNodeManager.AddKeyedDataField<FieldInfo::displacement>();
  
  if(m_recordIncrementalDisplacement){
    domain.m_feNodeManager.AddKeyedDataField<FieldInfo::incrementalDisplacement>();
  }
  
  domain.m_feNodeManager.AddKeylessDataField<int>("inContact",true,true);

  domain.m_feNodeManager.AddKeylessDataField<int>(m_trilinosIndexStr,true,true);

  domain.m_contactManager.AddKeylessDataField<int>( "activeInit", true, true);

  if(domain.m_contactManager.m_implicitContactActive == true)
  {
    domain.m_feFaceManager.AddKeylessDataField<R1Tensor>("contactStress", true, true);
    if(domain.m_contactManager.m_sliding_law != 0)
    {
      for (localIndex iLoad=0; iLoad<2; ++iLoad)
      {
        std::string s_base, s_gpNum;
        if(iLoad==0)
          s_base = "uJumpPl_Previous_gp";
        else
          s_base = "uJumpPl_Current_gp";
        for (localIndex iGp = 0; iGp < 4; ++iGp)
        {
          std::stringstream ss;
          ss << iGp;
          s_gpNum = ss.str();
          std::string FieldName = s_base + s_gpNum;
          domain.m_externalFaces.AddKeylessDataField<R1Tensor>(FieldName, true, true);
        }
      }
    }
  }

  if (m_hasReferenceStress)
  {
    for( std::map< ElementManagerT::RegKeyType, ElementRegionT >::iterator i = domain.m_feElementManager.m_ElementRegions.begin() ;
        i != domain.m_feElementManager.m_ElementRegions.end(); ++i )
    {
      ElementRegionT& elementRegion = i->second;
      elementRegion.AddKeylessDataField<R2SymTensor>("referenceStress", true, true);
    }
  }
  for(localIndex i = 0; i < m_thermalRegionNames.size(); ++i)
  {
    std::map<std::string, ElementRegionT>::iterator it = domain.m_feElementManager.m_ElementRegions.find(m_thermalRegionNames[i]);
    if (it == domain.m_feElementManager.m_ElementRegions.end() )
      throw GPException("Cannot find the thermal region specified!");

    ElementRegionT& elemRegion = it -> second;

    elemRegion.AddKeylessDataField<realT>("temperature", true, true);
    elemRegion.AddKeylessDataField<realT>("linearCTE", true, false);
    // This is the linear coefficient of thermal expansion, not the volumetric one.
    elemRegion.AddKeylessDataField<realT>("antiThermalStress", false, false);
    // This is the hydrostatic stress that would cause the same amount of volumetric strain as the thermal expantion/contraction would

    rArray1d& temperature = elemRegion.GetFieldData<realT>("temperature");
    rArray1d& CTE = elemRegion.GetFieldData<realT>("linearCTE");
    temperature = m_refTemperature;
    CTE = m_defaultCTE;

    //We don't register a fied for reference temperature.  If the user needs varying ref temp in the domain, he/she will initialize the field in initial conditions, which will automatically create the field.

  }
}



void LagrangeSmallStrain::Initialize(PhysicalDomainT& domain, SpatialPartition& partition )
{
  LagrangeSolverBase::Initialize(domain, partition);
  using namespace BoundaryConditionFunctions;
  ApplyMultiSetBoundaryCondition<R1Tensor>(this, &LagrangeSmallStrain::NonpenetratingBC_NeighborUpdate,
                                           domain, domain.m_feNodeManager,
                                           NonPenetratingBoundaryCondition::BoundaryConditionName(), 0.0 );
}



void LagrangeSmallStrain::InitializeCommunications( PartitionBase& partition )
{
  m_syncedFields.clear();

  if( m_enableTimeIntegrationOption[QuasiStatic] )
  {
    m_syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back(Field<FieldInfo::displacement>::Name());
//    m_syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back(Field<FieldInfo::incrementalDisplacement>::Name());
//    m_syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back(Field<FieldInfo::velocity>::Name());
    m_syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back(m_trilinosIndexStr);

  }

  if( m_enableTimeIntegrationOption[ImplicitDynamic] || m_enableTimeIntegrationOption[ExplicitDynamic] )
  {
//    m_syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back(m_trilinosIndexStr);
    m_syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back(Field<FieldInfo::acceleration>::Name());
    m_syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back(Field<FieldInfo::velocity>::Name());
    m_syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back(Field<FieldInfo::displacement>::Name());
//    m_syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back(Field<FieldInfo::incrementalDisplacement>::Name());


    m_syncedFields[PhysicalDomainT::FiniteElementFaceManager].push_back("cohesiveTraction");
    m_syncedFields[PhysicalDomainT::FiniteElementFaceManager].push_back("ruptureState");

  }

  if( m_enableTimeIntegrationOption[ExplicitDynamic] )
  {
    m_syncedFields[PhysicalDomainT::DiscreteElementNodeManager].push_back(Field<FieldInfo::acceleration>::Name());
    m_syncedFields[PhysicalDomainT::DiscreteElementNodeManager].push_back(Field<FieldInfo::velocity>::Name());
    m_syncedFields[PhysicalDomainT::DiscreteElementNodeManager].push_back(Field<FieldInfo::displacement>::Name());

  }

  partition.SetBufferSizes(m_syncedFields, CommRegistry::lagrangeSolver02);
}







void LagrangeSmallStrain::SetInitialGuess( const PhysicalDomainT& domain,
                                                realT* const local_solution )
{
/*
  const iArray1d& trilinos_index     = domain.m_feNodeManager.GetFieldData<int>(m_trilinosIndexStr);
  const Array1dT<R1Tensor>& incdisp  = domain.m_feNodeManager.GetFieldData<FieldInfo::incrementalDisplacement>();
  const iArray1d& is_ghost           = domain.m_feNodeManager.GetFieldData<FieldInfo::ghostRank>();

//  std::cout<<"LagrangeSmallStrain::SetInitialGuess"<<std::endl;

  for(unsigned r=0; r<incdisp.size(); ++r)
  {
    if(is_ghost[r] < 0)
    {
      for( int d=0 ; d<dim ; ++d )
      {
        int lid = m_rowMap->LID(dim*trilinos_index[r]+d);
//        local_solution[lid] = incdisp[r][d];
//        std::cout<<local_solution[lid]<<std::endl;
//        local_solution[lid] = 0.0;
      }
    }
  }
  */
}


void LagrangeSmallStrain::PropagateSolution( const realT* const local_solution,
                                                  const realT scalingFactor,
                                                  PhysicalDomainT& domain,
                                                  const localIndex dofOffset )
{
  const iArray1d& trilinos_index = domain.m_feNodeManager.GetFieldData<int>(m_trilinosIndexStr);
  const iArray1d& is_ghost       = domain.m_feNodeManager.GetFieldData<FieldInfo::ghostRank>();
  Array1dT<R1Tensor>& X = domain.m_feNodeManager.GetFieldData<FieldInfo::referencePosition>();
  Array1dT<R1Tensor>& disp = domain.m_feNodeManager.GetFieldData<FieldInfo::displacement>();
//  Array1dT<R1Tensor>& velocity = domain.m_feNodeManager.GetFieldData<FieldInfo::velocity>();
  Array1dT<R1Tensor>& incdisp  = domain.m_feNodeManager.GetFieldData<FieldInfo::incrementalDisplacement>();

  realT maxpos = 0.0;
  realT maxdisp = 0.0;

#if 0
  rArray1d const * const contactOffset = domain.m_feFaceManager.GetFieldDataPointer<realT>("contactOffset" );

  const iArray1d& faceGhostRank = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();
  const OrderedVariableOneToManyRelation& childFaceIndex = domain.m_feFaceManager.GetVariableOneToManyMap( "childIndices" );
  const iArray1d& flowFaceType = domain.m_feFaceManager.GetFieldData<int>("flowFaceType");

  Array1dT<R1Tensor> solnMod(domain.m_feNodeManager.m_numNodes);
  solnMod = 0.0;
  for(localIndex kf = 0 ; kf < domain.m_feFaceManager.m_numFaces ; ++kf)
  {
    if( faceGhostRank[kf] < 0 && flowFaceType[kf]==0 )
    {

      const localIndex faceIndex[2] = { kf, childFaceIndex[kf][0] };
      const R1Tensor N[2] = { domain.m_feFaceManager.FaceNormal( domain.m_feNodeManager, faceIndex[0] ),
                              domain.m_feFaceManager.FaceNormal( domain.m_feNodeManager, faceIndex[1] )};

      R1Tensor Nbar = N[0];
      Nbar -= N[1];
      Nbar.Normalize();
      const localIndex numFaceNodes = domain.m_feFaceManager.m_toNodesRelation[kf].size();

      for( localIndex a=0 ; a<numFaceNodes ; ++a )
      {

        // rhs and mass derivative
        const localIndex aa = a == 0 ? a : numFaceNodes - a;
        const localIndex localNodeIndex1 = domain.m_feFaceManager.m_toNodesRelation[faceIndex[0]][a];
        const localIndex localNodeIndex2 = domain.m_feFaceManager.m_toNodesRelation[faceIndex[1]][aa];

        R1Tensor wvec = disp[localNodeIndex2];
        wvec -= disp[localNodeIndex1];
        realT wnode0 = Dot(wvec,Nbar);


        if( contactOffset!=nullptr )
        {
          if( (*contactOffset)[kf] > m_ppSolve->m_contactOffsetCutoff )
          {
            wnode0 += (*contactOffset)[kf];
          }
          else
          {
            wnode0 += m_ppSolve->m_contactOffsetCutoff;
          }
        }


//        wnode0 += m_ppSolve->zeroApertureOffset();

        bool preContact = false;
        if( wnode0 < 0.0 )
        {
          preContact = true;
        }

        realT wnode1 = wnode0;
        for( int d=0 ; d<dim ; ++d )
        {
          int intDofOffset = dofOffset;
          int lid2 = m_rowMap->LID(intDofOffset+dim*trilinos_index[localNodeIndex2]+d);
          int lid1 = m_rowMap->LID(intDofOffset+dim*trilinos_index[localNodeIndex1]+d);
          realT inc = local_solution[lid2] - local_solution[lid1];
          wnode1 += inc * Nbar[d];
        }
        bool postContact = false;
        if( wnode1 < 0.0 )
        {
          postContact = true;
        }

        if( preContact==true && postContact==false )
        {
          R1Tensor temp(Nbar);
          temp *= wnode1;
          solnMod = temp;
        }
        else if( preContact==false && postContact==true )
        {
          R1Tensor temp(Nbar);
          temp *= wnode1;
          solnMod = temp;
        }
      }
    }
  }

#endif

  for(unsigned r=0; r<disp.size(); ++r)
  {
    if(is_ghost[r] < 0)
    {
      for( int d=0 ; d<dim ; ++d )
      {
        int intDofOffset = dofOffset;
        int lid = m_rowMap->LID(intDofOffset+dim*trilinos_index[r]+d);

        incdisp[r][d] -= scalingFactor*local_solution[lid];
        disp[r][d] -= scalingFactor*local_solution[lid];
        maxpos = std::max( maxpos, fabs(X[r][d]+disp[r][d]) );
        maxdisp = std::max( maxdisp, fabs(disp[r][d]) );
      }
    }
  }
  m_maxDofVal = maxpos;
//  std::cout<<"Maximum DeltaDisplacement, Position = "<<maxinc<<", "<<maxpos<<", "<<maxinc/maxpos<<std::endl;

}






















// Boundary Conditions
///////////////////////


/**
 * @author walsh24
 * @brief Apply traction boundary condition to a given face set
 * 
 */

void LagrangeSmallStrain::TractionBC( PhysicalDomainT& domain,
                                           ObjectDataStructureBaseT& object ,
                                           BoundaryConditionBase* bc,
                                           const lSet& set, realT time )
{
   
  TractionBoundaryCondition* trbc = bc->UpcastActiveBCPointer<TractionBoundaryCondition>(time); //   = dynamic_cast<TractionBoundaryCondition*> ( bc);
  if( trbc )
  {
    iArray1d& trilinos_index = domain.m_feNodeManager.GetFieldData<int>(m_trilinosIndexStr);
//    iArray1d& face_is_ghost  = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();
    iArray1d& node_is_ghost  = domain.m_feNodeManager.GetFieldData<FieldInfo::ghostRank>();
   
    Epetra_IntSerialDenseVector  node_dof(1);
    Epetra_SerialDenseVector     node_rhs(1);
    
//    const R1Tensor& n = trbc->GetDirection(time); // old way
  
    for( lSet::const_iterator fc=set.begin() ; fc!=set.end() ; ++fc )
    {
      if( domain.m_feFaceManager.m_toElementsRelation[*fc].size() > 0)
      {

        const realT area = domain.m_feFaceManager.SurfaceArea( domain.m_feNodeManager, *fc, true );
        
/*  Old way
        realT value = trbc->GetValue(domain.m_feFaceManager,fc,time);
        R1Tensor traction = n;
        if( trbc->IsNormalTraction() )
        {
          traction = domain.m_feFaceManager.FaceNormal( domain.m_feNodeManager, *fc );
          traction *= -1;
        }
        traction *= value * area / domain.m_feFaceManager.m_toNodesRelation[*fc].size();
*/
        R1Tensor traction = trbc->GetTractionOnFace(domain,fc,time) * area / domain.m_feFaceManager.m_toNodesRelation[*fc].size();
   
 /*       R1Tensor traction = trbc->GetTractionOnFace(domain,fc,time); // new way
        traction *= area / domain.m_feFaceManager.m_toNodesRelation[*fc].size();*/

        for( lArray1d::const_iterator nd=domain.m_feFaceManager.m_toNodesRelation[*fc].begin() ;
             nd!=domain.m_feFaceManager.m_toNodesRelation[*fc].end() ; ++nd )
        {
          if (node_is_ghost[*nd] < 0)
          {
            for(int ii =0; ii < dim;++ii)
            {
              node_dof(0) = dim*trilinos_index[*nd]+ii;
              node_rhs(0) = traction(ii);
              m_rhs->SumIntoGlobalValues(node_dof, node_rhs);
            }

          }

        }
      } 
    }
  }
}

/**
 * @author walsh24
 * @brief Apply pressure boundary condition to a given face set
 *
 */

void LagrangeSmallStrain::PressureBC( PhysicalDomainT& domain,
                                               ObjectDataStructureBaseT& ,
                                               BoundaryConditionBase*,
                                               const lSet& set, realT )
{

  iArray1d& trilinos_index = domain.m_feNodeManager.GetFieldData<int>(m_trilinosIndexStr);
  iArray1d& face_is_ghost  = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();
  rArray1d& pressures = domain.m_feFaceManager.GetFieldData<realT>(PS_STR::PressureStr);
  rArray1d* apertures;

  Epetra_IntSerialDenseVector  node_dof(1);
  Epetra_SerialDenseVector     node_rhs(1);

  // If aperture is defined and is <= 0, the pressure BC is not enforced.
  bool isApertureControlled = false;
  if( domain.m_feFaceManager.HasField<realT>(PS_STR::ApertureStr) ){
    isApertureControlled = true;
    apertures = &(domain.m_feFaceManager.GetFieldData<realT>(PS_STR::ApertureStr));
  }

  for( lSet::const_iterator fc=set.begin() ; fc!=set.end() ; ++fc ) {
    if( face_is_ghost[*fc] < 0 ){
      realT pressure = pressures[*fc];

      if(isApertureControlled){
        if ( (*apertures)[*fc] <= 0.0) break;
      }

      const realT area = domain.m_feFaceManager.SurfaceArea( domain.m_feNodeManager, *fc );

      R1Tensor traction = domain.m_feFaceManager.FaceNormal(domain.m_feNodeManager, *fc);
      traction *= (-pressure * area * 0.25); // TODO hardcoded for 4 node face

      for( lArray1d::const_iterator nd=domain.m_feFaceManager.m_toNodesRelation[*fc].begin() ;
          nd!=domain.m_feFaceManager.m_toNodesRelation[*fc].end() ; ++nd )
      {

        for(int ii =0; ii < dim;++ii){
          node_dof(0) = dim*trilinos_index[*nd]+ii;
          node_rhs(0) = traction(ii);
          m_rhs->SumIntoGlobalValues(node_dof, node_rhs);
        }
      }
    }
  }


}

/**
 * @author walsh24
 * @brief  Apply displacement boundary condition to a node set.
 * 
 */

void LagrangeSmallStrain ::
DisplacementBC(PhysicalDomainT& domain, ObjectDataStructureBaseT& object ,
               BoundaryConditionBase* bc, const lSet& set, realT time)
{

  realT LARGE ;

  iArray1d& trilinos_index = domain.m_feNodeManager.GetFieldData<int>(m_trilinosIndexStr);
  iArray1d& node_is_ghost  = domain.m_feNodeManager.GetFieldData<FieldInfo::ghostRank>();
  Array1dT<R1Tensor>& disp = domain.m_feNodeManager.GetFieldData<FieldInfo::displacement>();
  Epetra_IntSerialDenseVector  node_dof(1);
  Epetra_SerialDenseVector     node_rhs(1);

  int component = bc->GetComponent(time);
  if (component == 2 && domain.m_feFaceManager.m_toNodesRelation[0].size() == 2)
    throw GPException("This is a 2D problem but you are applying BC to component 2 (z).");
  for( lSet::const_iterator nd=set.begin() ; nd!=set.end() ; ++nd )
  {
    node_dof(0) = dim*trilinos_index[*nd]+component;

    if( node_is_ghost[*nd] < 0 )
    {
#if USECPP11==1
      LARGE = ClearRow( m_matrix.get(), node_dof(0), 1.0 );
      ClearRow( this->m_system->GetMatrix( EpetraBlock::solidBlock, EpetraBlock::fluidBlock ), node_dof(0),1.0);
#else
      LARGE = ClearRow( m_matrix, node_dof(0), 1.0 );
#endif
      node_rhs(0) = -LARGE*( bc->GetValue(domain.m_feFaceManager,nd,time) - (disp[*nd])[component] );
    }
    else
    {
      LARGE = ClearRow( m_matrix, node_dof(0), 0.0 );
      node_rhs(0) = 0.0;
    }

    if( m_system != NULL )
    {
      if( m_system->GetBlockID(EpetraBlock::solidBlock) != -1 && m_system->GetBlockID(EpetraBlock::fluidBlock) != -1 )
      {
        ClearRow( this->m_system->GetMatrix( EpetraBlock::solidBlock, EpetraBlock::fluidBlock ), node_dof(0),0.0);
      }
    }

    m_rhs->ReplaceGlobalValues(node_dof, node_rhs);
  }
}

/**
 * @author walsh24
 * @brief  Update neighbors for non-penetrating boundary condition. 
 * 
 */

void LagrangeSmallStrain ::
NonpenetratingBC_NeighborUpdate(PhysicalDomainT& domain, ObjectDataStructureBaseT& object ,
                                BoundaryConditionBase* bc, realT time  ){
               	
  NonPenetratingBoundaryCondition* npbc = dynamic_cast<NonPenetratingBoundaryCondition*> (bc);   	
  if(npbc) npbc->UpdateNearestNeighborMaps(domain);         	
               	
}

/**
 * @author walsh24
 * @brief  Fix displacements of penetrating neighbors - caution - assumes nodes are in close proximity at contact. 
 * 
 */

void LagrangeSmallStrain ::
NonpenetratingBC_DetectContact(PhysicalDomainT& domain, ObjectDataStructureBaseT& object, 
                               BoundaryConditionBase* bc, realT time  ){
  NonPenetratingBoundaryCondition* npbc = dynamic_cast<NonPenetratingBoundaryCondition*> (bc);   
  if(npbc){     
  
  	iArray1d& inContact = domain.m_feNodeManager.GetFieldData<int>("inContact");
  	
    Array1dT<R1Tensor>& disp = object.GetFieldData<FieldInfo::displacement> ();
    const Array1dT<R1Tensor>& X = object.GetFieldData<FieldInfo::referencePosition> ();
    
//    iArray1d& face_is_ghost  = domain.m_faceManager.GetFieldData<FieldInfo::ghostRank>();
    
  	lSet&  contactingNodes = npbc->m_contactingNodes;
  	std::map<localIndex,localIndex>&  nearestNodeNeighborMap = npbc->m_nearestNodeNeighborMap;
  	
  	std::map<localIndex,localIndex>::iterator itr = npbc->m_nearestFaceNeighborMap.begin();
  	std::map<localIndex,localIndex>::iterator iend = npbc->m_nearestFaceNeighborMap.end();
  	for(; itr != iend; ++itr){
  	  localIndex fc  = itr->first;
//  	  localIndex fc_nbr = itr->second;
  	 
  	  // loop over nodes in face
  	  lArray1d& nodeList = domain.m_feFaceManager.m_toNodesRelation[fc];
  	  R1Tensor n = domain.m_feFaceManager.FaceNormal(domain.m_feNodeManager, fc);
  	  for(lArray1d::size_type i =0 ; i< nodeList.size(); ++i){
  	 	
  	    localIndex nd = nodeList[i];
  	    localIndex nbr = nearestNodeNeighborMap[nd];
  	   
  	    R1Tensor pa = disp[nd] + X[nd];
  	    R1Tensor l = disp[nbr] + X[nbr] - pa;  // branch vector
  	       	 
  	    if(Dot(n,l) <=  0.0  && (nearestNodeNeighborMap[nbr] == nd) ){
  	 	  contactingNodes.insert(nd);
  	 	  contactingNodes.insert(nbr);
  	 	    
  	 	  l*=0.5;
  	 	  disp[nd] += l;
  	 	  disp[nbr] -= l;
  	 	    
   	      inContact[nd] = nbr;
   	      inContact[nbr] = nd;
   	    }
  	  }
    }
  }                     	
}


void LagrangeSmallStrain ::
NonpenetratingBC_UpdateAperture(PhysicalDomainT& domain, ObjectDataStructureBaseT& object ,
                               BoundaryConditionBase* bc, realT time  ){
  NonPenetratingBoundaryCondition* npbc = dynamic_cast<NonPenetratingBoundaryCondition*> (bc);   
  if(npbc){     
  
  	rArray1d& aperture = domain.m_feFaceManager.GetFieldData<realT>("Aperture");
  	
  	
  	std::map<localIndex,localIndex>::iterator itr = npbc->m_nearestFaceNeighborMap.begin();
  	std::map<localIndex,localIndex>::iterator iend = npbc->m_nearestFaceNeighborMap.end();
  	for(; itr != iend; ++itr){
  	  localIndex fc  = itr->first;
  	  localIndex fc_nbr = itr->second;
  	  
  	  //if(npbc->m_nearestFaceNeighborMap[fc_nbr] == fc){
  	  	
        R1Tensor norm  = domain.m_feFaceManager.FaceNormal(domain.m_feNodeManager, fc);
        R1Tensor normB = domain.m_feFaceManager.FaceNormal(domain.m_feNodeManager, fc_nbr);
        R1Tensor center; domain.m_feFaceManager.FaceCenter(domain.m_feNodeManager, fc,center);
        R1Tensor centerB; domain.m_feFaceManager.FaceCenter(domain.m_feNodeManager, fc_nbr,centerB);
        norm -= normB;
        norm.Normalize();
        
  	    aperture[fc] = std::max((centerB-center)*norm,0.0);
  	  //}
    }
  }                     	
}

/**
 * @author walsh24
 * @brief  Fix displacements of penetrating neighbors - caution - assumes nodes are in close proximity. 
 * 
 */

void LagrangeSmallStrain ::
NonpenetratingBC_Sparsity(PhysicalDomainT& domain, ObjectDataStructureBaseT& object ,
                          BoundaryConditionBase* bc, realT time  ){
  

  NonPenetratingBoundaryCondition* npbc = dynamic_cast<NonPenetratingBoundaryCondition*> (bc);   
  if(npbc){
  	
    iArray1d& trilinos_index = domain.m_feNodeManager.GetFieldData<int>(m_trilinosIndexStr);	
    iArray1d& node_is_ghost  = domain.m_feNodeManager.GetFieldData<FieldInfo::ghostRank>();
    
    std::map<localIndex,localIndex>&  nearestNodeNeighborMap = npbc->m_nearestNodeNeighborMap;
    
    std::vector<int>  node_row_dof(1);
    std::vector<int>  node_row_dofB(dim);
    std::vector<int>  node_col_dofB(2*dim);
    
    // give contacting nodes the same trilinos index, 
    // add other index to dummy trilinos indicies.
    lSet&  contactingNodes = npbc->m_contactingNodes;
    lSet::iterator iend = contactingNodes.end();
    for( lSet::iterator itr = contactingNodes.begin(); itr != iend; ++itr ) {	
      localIndex nd =  *itr;
      localIndex nbr = nearestNodeNeighborMap[nd];
      
      if(trilinos_index[nbr] > trilinos_index[nd]){
          
        int triNbr = trilinos_index[nbr];
        trilinos_index[nbr] = trilinos_index[nd];

        for(int component =0; component < dim; ++component){
          node_row_dof[0] = dim*triNbr+component;
          dummyDof.push_back(node_row_dof[0]);
          m_sparsity->InsertGlobalIndices(node_row_dof.size(),
                                        &node_row_dof.front(),
                                        node_row_dof.size(),
                                        &node_row_dof.front());
      	}
     }
     
    }
    //std::cout << "Number of Contacting Nodes: " <<  contactingNodes.size() <<std::endl;
    
    // add resistance between non contacting nodes.

    std::map<localIndex,localIndex>::iterator iendB = nearestNodeNeighborMap.end();
    for( std::map<localIndex,localIndex>::iterator itr = nearestNodeNeighborMap.begin();
         itr != iendB; ++itr ) {	
      localIndex nd =  itr->first;
      if(node_is_ghost[nd] < 0){
        localIndex nbr = nearestNodeNeighborMap[nd];
        int tri = trilinos_index[nd];
        int triNbr = trilinos_index[nbr];
        if(tri != triNbr){  // ie not in contact
          for(int component =0; component < dim; ++component){
            node_row_dofB[component] = dim*tri+component;
            node_col_dofB[component] = node_row_dofB[component];
            node_col_dofB[dim+component] = dim*triNbr+component;
          }
          m_sparsity->InsertGlobalIndices(node_row_dofB.size(),
                                        &node_row_dofB.front(),
                                        node_col_dofB.size(),
                                        &node_col_dofB.front());

        }
      }
    }
    	         	
  }             	
}

/**
 * @author walsh24
 * @brief  Fix displacements of penetrating neighbors - caution - assumes neighbors in close proximity. 
 * 
 */

void LagrangeSmallStrain ::
NonpenetratingBC_Apply(PhysicalDomainT& domain, ObjectDataStructureBaseT& object ,
               BoundaryConditionBase* bc, realT time  ){


  NonPenetratingBoundaryCondition* npbc = dynamic_cast<NonPenetratingBoundaryCondition*> (bc);   
  if(npbc){
//    const Array1dT<R1Tensor>& disp = object.GetFieldData<FieldInfo::displacement> ();
  	
    iArray1d& trilinos_index = domain.m_feNodeManager.GetFieldData<int>(m_trilinosIndexStr);	
//    iArray1d& node_is_ghost  = domain.m_nodeManager.GetFieldData<FieldInfo::ghostRank>();
    iArray1d& face_is_ghost  = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();
  	
    Epetra_IntSerialDenseVector  node_row_dof(1);
    Epetra_IntSerialDenseVector  node_row_dofB(dim);
    Epetra_IntSerialDenseVector  node_col_dofB(2*dim);
    Epetra_SerialDenseVector     node_rhsB(dim);
    Epetra_SerialDenseMatrix     node_matrix  (1,1);
    Epetra_SerialDenseMatrix     node_matrixB  (dim,2*dim);
  	
    // Set diagonal of dummy dofs to 1
    node_matrix(0,0) = 1.0;
    for( iArray1d::size_type i=0; i < dummyDof.size() ; ++i ) {
      node_row_dof(0) = dummyDof[i];	
      m_matrix->ReplaceGlobalValues(node_row_dof,node_matrix);
    }       	

    // add entries to represent springs between non contacting nodes.
    std::map<localIndex,localIndex>&  nearestFaceNeighborMap = npbc->m_nearestFaceNeighborMap;
    std::map<localIndex,localIndex>&  nearestNodeNeighborMap = npbc->m_nearestNodeNeighborMap;
    std::map<localIndex,localIndex>::iterator iendB = nearestFaceNeighborMap.end();
    for( std::map<localIndex,localIndex>::iterator itr = nearestFaceNeighborMap.begin();
         itr != iendB; ++itr ) {	
      localIndex fc =  itr->first;	
      
      if(face_is_ghost[fc] < 0){
        localIndex fcb =  itr->second;
        R1Tensor norm  = domain.m_feFaceManager.FaceNormal(domain.m_feNodeManager, fc);
        R1Tensor normB = domain.m_feFaceManager.FaceNormal(domain.m_feNodeManager, fcb);
        R1Tensor center; domain.m_feFaceManager.FaceCenter(domain.m_feNodeManager, fc,center);
        R1Tensor centerB; domain.m_feFaceManager.FaceCenter(domain.m_feNodeManager, fcb,centerB);
        norm -= normB;
        norm.Normalize();

        R1Tensor dl = centerB-center;
 
//        realT l = sqrt(Dot(dl,dl));
        realT a = 0.5*(domain.m_feFaceManager.SurfaceArea(domain.m_feNodeManager, fc)
                       + domain.m_feFaceManager.SurfaceArea(domain.m_feNodeManager, fcb));

        realT k =  0.25*a*m_nonContactModulus; // FIXME /l;    // spring stiffness
   
  	    lArray1d& nodeList = domain.m_feFaceManager.m_toNodesRelation[fc];
  	    for( lArray1d::size_type nn =0 ; nn< nodeList.size(); ++nn){
          localIndex nd = nodeList[nn];
          localIndex nbr = nearestNodeNeighborMap[nd];
          
          int tri = trilinos_index[nd];
          int triNbr = trilinos_index[nbr];
          if(tri != triNbr){  // not in contact

           // R1Tensor du = disp[nd] - disp[nbr];
           // R1Tensor f_rhs = norm; f_rhs *= k*Dot(norm,du);

            for(int i =0; i < dim; ++i){
              node_row_dofB(i) = dim*tri+i;
              node_col_dofB(i) = node_row_dofB(i);
              node_col_dofB(dim+i) = dim*triNbr+i;
              for(int j =0; j < dim; ++j){
                realT val = k*(0.1+0.9*norm[i]*norm[j]);
                node_matrixB(i,j) = val;
                node_matrixB(i,j+dim) = -val;
              }
             // node_rhsB(i) = f_rhs[i];
            }
            m_matrix->SumIntoGlobalValues(node_row_dofB,node_col_dofB,node_matrixB);
           // m_rhs->SumIntoGlobalValues(node_row_dofB,node_rhsB);
          }
        }
      }
    }

    if(npbc->m_updatePressure){

      // loop over second set - update pressure based on nearest neighbor from first
      rArray1d& pressures = domain.m_feFaceManager.GetFieldData<realT>(PS_STR::PressureStr);

      const lSet& setB = domain.m_feFaceManager.m_Sets[npbc->m_setNames[1] ];
      for( lSet::const_iterator fcB=setB.begin(); fcB!=setB.end() ; ++fcB ) {
        localIndex nbr = nearestFaceNeighborMap[*fcB];
        pressures[*fcB] = pressures[nbr];
      }
    }

  }             	
}



void LagrangeSmallStrain::FixNodesBC( NodeManagerT const& nodeManager, const lSet& set )
{
  iArray1d const & trilinos_index = nodeManager.GetFieldData<int>(m_trilinosIndexStr);
  iArray1d const & ghostRank  = nodeManager.GetFieldData<FieldInfo::ghostRank>();

  Epetra_IntSerialDenseVector  node_dof(1);
  Epetra_SerialDenseVector     node_rhs(1);

  for( lSet::const_iterator nd=set.begin() ; nd!=set.end() ; ++nd )
  {
    if( ghostRank[*nd] < 0 )
    {
//      std::cout<<"fixing node "<<*nd<<" dof "<<dim*trilinos_index[*nd]<<std::endl;
      for( int d=0 ; d<dim ; ++d )
      {
        node_dof(0) = dim*trilinos_index[*nd]+d;
        node_rhs(0) = 0.0;
#if USECPP11==1
        ClearRow( m_matrix.get(), node_dof(0), 1.0 );
        ClearRow( this->m_system->GetMatrix( EpetraBlock::solidBlock, EpetraBlock::fluidBlock ), node_dof(0),1.0);
#else
        ClearRow( m_matrix, node_dof(0), 1.0 );
        if( m_system != NULL )
        {
          if( m_system->GetBlockID(EpetraBlock::solidBlock) != -1 && m_system->GetBlockID(EpetraBlock::fluidBlock) != -1 )
          {
            ClearRow( this->m_system->GetMatrix( EpetraBlock::solidBlock, EpetraBlock::fluidBlock ), node_dof(0),0.0);
          }
        }
#endif
        m_rhs->ReplaceGlobalValues(node_dof, node_rhs);
      }
    }
  }
}



/* Explicit Instantiations */


/* Register solver in the solver factory */

//SolverRegistrator<LagrangeSmallStrain<3> > reg_LagrangeSmallStrain;

//SolverRegistrator<LagrangeSmallStrain<2> >
//  reg_ImplicitLaplaceSolver2D;



