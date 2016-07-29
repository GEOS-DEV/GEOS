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
 * @file ParallelPlateFlowSolverExplicit.cpp
 * @author settgast1
 * @date Feb 10, 2011
 */

#include "ParallelPlateFlowSolverExplicit.h"
#include "SolverFactory.h"
#include <algorithm>
#include <numeric>

#include "BoundaryConditions/ApplyBoundaryConditions.h"

#ifdef SRC_EXTERNAL
#include "BoundaryConditions/PerforatedCasedWellboreBoundaryCondition.h"
#include "BoundaryConditions/CavityPressureBoundaryCondition.h"
#endif

#include "PhysicsSolverStrings.h"
using namespace PS_STR;
using namespace PPFS;

namespace {
  realT TINY = 1e-64;
}

ParallelPlateFlowSolverExplicit::ParallelPlateFlowSolverExplicit( const std::string& name,
                                                  ProblemManagerT* const pm ):
    ParallelPlateFlowSolverBase(name,pm),
    m_bBarton((1e-3 - 1e-4)/1e7/1e-4),
    m_aBarton(1e-3 * (1e-3 - 1e-4)/1e7/1e-4),
    m_wZeroStress(1e-3)
{}

ParallelPlateFlowSolverExplicit::~ParallelPlateFlowSolverExplicit()
{
}

void ParallelPlateFlowSolverExplicit::ReadXML( TICPP::HierarchicalDataNode* const hdn  )
{
  ParallelPlateFlowSolverBase::ReadXML( hdn );

  // if a separate courant factor is given for the parallel plate flow solver
  m_courant = hdn->GetAttributeOrDefault<realT>("ppcourant",m_courant);
  m_farFieldPorePressure = hdn->GetAttributeOrDefault<realT>("farFieldPorePressre", 0.0);
  m_pressureDependentLeakoff = hdn->GetAttributeOrDefault<int>("pressureDependentLeakoff", 0);
  m_apertureMovingAverageCoeff = hdn->GetAttributeOrDefault<realT>("apertureMovingAverageCoeff", 0.0);

  std::string tempString = hdn->GetAttributeString("leakoffCoefficient");
  if(!tempString.empty())
    throw GPException("Error! leakoffCoefficient is obsolete. You must use the new name (CartersLeakoffCoefficient) instead, which is the coefficient per fracture surface area, consistent with engineering conventions.");
  m_leakoffCoef = hdn->GetAttributeOrDefault<realT>("CartersLeakoffCoefficient", 0.0);

  m_overLeakCompensation = hdn->GetAttributeOrDefault<int>("overLeakCompensation", 0);

  std::string temp = hdn->GetAttributeString("BartonJointParameters"); // aperture at zero effective stress; reference stress; aperture at ref stress
  if( !temp.empty() )
  {
    R1Tensor tempArray;
    tempArray.StrVal( temp );
    m_wZeroStress = tempArray[0];
    realT stressRef = tempArray[1];
    realT wRef = tempArray[2];

    if (wRef < 0.99 * m_wZeroStress)
    {
      m_bBarton = (m_wZeroStress - wRef) / stressRef / wRef;
      m_aBarton = m_wZeroStress * (m_wZeroStress - wRef) / stressRef / wRef;
    }
    else
    {
      m_bBarton = 0;
      m_aBarton = 0;
      m_min_aperture = m_wZeroStress;
    }
  }
  else
  {
    m_bBarton = 0;
    m_aBarton = 0;
  }


}


void ParallelPlateFlowSolverExplicit::RegisterFields( PhysicalDomainT& domain )
{
  ParallelPlateFlowSolverBase::RegisterFields( domain.m_feFaceManager, domain.m_feEdgeManager );

  domain.m_feEdgeManager.AddKeylessDataField<realT>("massRate",true,true);
  domain.m_feFaceManager.AddKeyedDataField<FieldInfo::volume>();
  domain.m_feFaceManager.AddKeylessDataField<realT>("elementFlowRate", false, true);
  domain.m_feFaceManager.AddKeylessDataField<R1Tensor>("flowVector", false, true);

  // domain.m_feFaceManager.AddKeylessDataField<realT>("apertureOffset", false, true);
  // This field was added for a specific project.  Actually we don't really need it here.  The initial condition manager will add the field if it is called in the xml input.

  domain.m_feFaceManager.AddKeylessDataField<int>("tracer", true, true );

  domain.m_feFaceManager.AddKeylessDataField<realT>( "effectiveStressN",true, true );
  domain.m_feFaceManager.AddKeylessDataField<realT>( "permeability",false, true );

  domain.m_feNodeManager.AddKeylessDataField<realT>( "nodalPressure",false, true );

  if( !(m_flowFaceSetName.empty()) )
   {
     domain.m_feFaceManager.AddKeylessDataField<realT>( "stressNOnFace",true, true );
     domain.m_feFaceManager.AddKeylessDataField<R1Tensor>( "stressTOnFace",true, true );
   }




}


void ParallelPlateFlowSolverExplicit::Initialize( PhysicalDomainT& domain, SpatialPartition& partition )
{
  Array1dT<lSet>& edgeToFlowFaces = domain.m_feEdgeManager.GetUnorderedVariableOneToManyMap("edgeToFlowFaces");
  rArray1d& faceArea = domain.m_feFaceManager.GetFieldData<realT>("faceArea");
  iArray1d& flowFaceType = domain.m_feFaceManager.GetFieldData<int>("flowFaceType");
  iArray1d& flowEdgeType = domain.m_feEdgeManager.GetFieldData<int>("flowEdgeType");


  if( !(m_flowFaceSetName.empty()) )
  {

    const lSet& flowfaceset = domain.m_feFaceManager.GetSet(m_flowFaceSetName);

    for( lSet::const_iterator faceID=flowfaceset.begin() ; faceID!=flowfaceset.end() ; ++faceID )
    {
      flowFaceType[*faceID] = 1;
    }

    for( localIndex kf=0 ; kf<domain.m_feFaceManager.DataLengths() ; ++kf )
    {
      if( flowFaceType[kf] == 1 && domain.m_feFaceManager.m_parentIndex[kf] == LOCALINDEX_MAX )
      {
        for( localIndex ke=0 ; ke<domain.m_feFaceManager.m_toEdgesRelation[kf].size() ; ++ke )
        {
          const localIndex edgeIndex = domain.m_feFaceManager.m_toEdgesRelation[kf][ke];
          flowEdgeType[edgeIndex] = 1;

          const localIndex parentEdgeIndex = domain.m_feEdgeManager.GetParentIndex( edgeIndex );
          edgeToFlowFaces[parentEdgeIndex].insert( kf );
        }
      }
    }
  }

  // In 2D, we don't specify flowFaceSet.  We will prepare edgeToFlowFaces in Fractunator2D



  GenerateParallelPlateGeometricQuantities( domain,0,0);
  const rArray1d& fluidVolume  = domain.m_feFaceManager.GetFieldData<FieldInfo::volume>();
  rArray1d& faceFluidMass = domain.m_feFaceManager.GetFieldData<FieldInfo::mass>();
  rArray1d& faceFluidDensity = domain.m_feFaceManager.GetFieldData<FieldInfo::density>();

  for( localIndex kf=0 ; kf<domain.m_feFaceManager.DataLengths() ; ++kf )
  {
    faceFluidMass[kf] = fluidVolume[kf] * faceFluidDensity[kf];
    faceArea[kf] = domain.m_feFaceManager.SurfaceArea( domain.m_feNodeManager, kf );
  }


}





void ParallelPlateFlowSolverExplicit::InitializeCommunications( PartitionBase& partition )
{
  // ParallelPlateFlowSolverBase::InitializeCommunications( partition );

  std::map<PhysicalDomainT::ObjectDataStructureKeys, sArray1d> syncedFields_tmp;
  syncedFields_tmp[PhysicalDomainT::FiniteElementFaceManager].push_back(Field<FieldInfo::mass>::Name());
  syncedFields_tmp[PhysicalDomainT::FiniteElementFaceManager].push_back("fractureNodalPressure");
  syncedFields_tmp[PhysicalDomainT::FiniteElementFaceManager].push_back("isFractureNode");
  // synchronize element fields
  partition.SetBufferSizes( syncedFields_tmp, CommRegistry::parallelPlateFlowSolver);
}




double ParallelPlateFlowSolverExplicit::TimeStep( const realT& time,
                                        const realT& dt,
                                        const int cycleNumber,
                                        PhysicalDomainT& domain,
                                        const sArray1d& namesOfSolverRegions ,
                                        SpatialPartition& partition,
                                        FractunatorBase* const fractunator )
{
  realT dt_return = dt;

  GenerateParallelPlateGeometricQuantities( domain, time,dt );

  CalculateAndApplyMassFlux( dt, domain );

  ApplyFluxBoundaryCondition(time, dt, cycleNumber, partition.m_rank, domain);

  {
    std::map<PhysicalDomainT::ObjectDataStructureKeys, sArray1d> syncedFields_tmp;
    syncedFields_tmp[PhysicalDomainT::FiniteElementFaceManager].push_back(Field<FieldInfo::mass>::Name());

    // synchronize element fields
    partition.SynchronizeFields( syncedFields_tmp, CommRegistry::parallelPlateFlowSolver);
  }

  UpdateEOS( time, dt, domain );
  return dt_return;
}

void ParallelPlateFlowSolverExplicit::PostProcess(PhysicalDomainT & domain,
                                          SpatialPartition& partition,
                                          const sArray1d& namesOfSolverRegions)
{

  rArray1d& aperture = domain.m_feFaceManager.GetFieldData<realT>( ApertureStr );
  const Array1dT<R1Tensor>& faceCenter = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr );
  rArray1d& edgePermeability = domain.m_feEdgeManager.GetFieldData<realT>(PermeabilityStr);
  rArray1d& faceFluidPressure = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();
  rArray1d& faceFluidDensity = domain.m_feFaceManager.GetFieldData<FieldInfo::density>();
  //rArray1d& faceFluidMass = domain.m_feFaceManager.GetFieldData<FieldInfo::mass>();
  // rArray1d& flowRate = domain.m_feFaceManager.GetFieldData<realT>("elementFlowRate");
  // We do this in every timestep anyway.  Because we don't call the flow bc here, the value calculated here will be missing the contribution from the flow bc anyway.

  Array1dT<R1Tensor>& flowVector = domain.m_feFaceManager.GetFieldData<R1Tensor>("flowVector");
  rArray1d& massRate = domain.m_feEdgeManager.GetFieldData<realT>("massRate");
  const iArray1d& flowFaceType = domain.m_feFaceManager.GetFieldData<int>("flowFaceType");
  rArray1d& p = domain.m_feNodeManager.GetFieldData<realT>("nodalPressure");
//  rArray1d& faceArea = domain.m_feFaceManager.GetFieldData<realT>("faceArea");




  rArray1d& edgeLength = domain.m_feEdgeManager.GetFieldData<realT>("length");
  const Array1dT<R1Tensor>& edgeCenter = domain.m_feEdgeManager.GetFieldData<R1Tensor>("center");


  //  const iArray1d& flowFaceType = domain.m_feFaceManager.GetFieldData<int>("flowFaceType");
  const iArray1d& flowEdgeType = domain.m_feEdgeManager.GetFieldData<int>("flowEdgeType");


  const Array1dT<lSet>& edgeToFlowFaces = domain.m_feEdgeManager.GetUnorderedVariableOneToManyMap("edgeToFlowFaces");

  R1Tensor la, lb;


  for( localIndex kf=0 ; kf<domain.m_feFaceManager.DataLengths() ; ++kf )
  {
    flowVector[kf] *= 0.0;
  }


  // Calculate flow vectors.
  for( localIndex ke=0 ; ke<domain.m_feEdgeManager.DataLengths() ; ++ke )
  {

    if( flowEdgeType[ke]== 1 )
    {
      realT w = edgeLength[ke];

      if(w == 0.0){    // Fix bug where edgelength is not calculated in some cases
        w = domain.m_feEdgeManager.EdgeLength( domain.m_feNodeManager, ke );
        edgeLength[ke] = w;
      }


      // we have two different ways of calculating the flow. One for edges attached to only 2
      // faces, and one for more than two faces.
      const unsigned int numFlowFaces = edgeToFlowFaces[ke].size();

      if( numFlowFaces == 2 )
      {
        lSet::const_iterator edgeToFlowFace = edgeToFlowFaces[ke].begin();
        const localIndex faceIndex0 = *edgeToFlowFace;
        const localIndex faceIndex1 = *(++edgeToFlowFace);


        la = edgeCenter[ke];
        la -= faceCenter[faceIndex0];

        lb = edgeCenter[ke];
        lb -= faceCenter[faceIndex1];

        realT norm_la = la.L2_Norm();
        realT norm_lb = lb.L2_Norm();

        realT permeabilityAperture0 = std::min(aperture[faceIndex0], m_max_aperture);
        realT permeabilityAperture1 = std::min(aperture[faceIndex1], m_max_aperture);

        edgePermeability[ke] = CalculatePermeability( norm_la , norm_lb,
                                                      permeabilityAperture0, permeabilityAperture1,
                                                      w,m_mu,m_SHP_FCT);

        realT PRhoGravity = CalculatePRhoGravity(faceCenter[faceIndex0], faceCenter[faceIndex1],
                                                    faceFluidDensity[faceIndex0], faceFluidDensity[faceIndex1],
                                                    m_gravityVector);
        // determine the mass flux across edge at t_n+1/2
        massRate[ke] = edgePermeability[ke] * ( -( faceFluidDensity[faceIndex0] * faceFluidPressure[faceIndex0])
            + ( faceFluidDensity[faceIndex1] * faceFluidPressure[faceIndex1] )
            - PRhoGravity);

        la.Normalize();
        lb.Normalize();

        flowVector[faceIndex0] -= la * massRate[ke] / std::max(m_rho_o, faceFluidDensity[faceIndex0]) * 0.5 / w; //faceArea[faceIndex0] ;
        flowVector[faceIndex1] += lb * massRate[ke] / std::max(m_rho_o, faceFluidDensity[faceIndex1]) * 0.5 / w; //faceArea[faceIndex1] ;

      }
      else if( numFlowFaces > 2 )
      {
        realT rhoP = 0.0;
        realT sumK = 0.0;
        Array1dT<R1Tensor> length(numFlowFaces);
        rArray1d k(numFlowFaces);
        rArray1d q(numFlowFaces);
        rArray1d kRhoP(numFlowFaces);
        rArray1d PRhoGravity(numFlowFaces);

        lSet::const_iterator faceIndex=edgeToFlowFaces[ke].begin();

        for( localIndex kf=0 ; kf<numFlowFaces ; ++kf, ++faceIndex)
        {

          length[kf] = edgeCenter[ke];
          length[kf] -= faceCenter[*faceIndex];

          realT permeabilityAperture = std::min(aperture[*faceIndex], m_max_aperture);

          k[kf] = CalculatePermeability( length[kf].L2_Norm(),
                                         permeabilityAperture,
                                         w, m_mu, m_SHP_FCT );

          PRhoGravity[kf] = CalculatePRhoGravity( faceCenter[*faceIndex], edgeCenter[ke], faceFluidDensity[*faceIndex], m_gravityVector);

          sumK += k[kf];

          kRhoP[kf] = k[kf] * faceFluidDensity[*faceIndex] * faceFluidPressure[*faceIndex];
          rhoP += kRhoP[kf];



        }


        rhoP /= (sumK+TINY);
        faceIndex=edgeToFlowFaces[ke].begin();

        for( localIndex kf=0 ; kf<numFlowFaces ; ++kf, ++faceIndex )
        {
          q[kf] = k[kf] * ( faceFluidDensity[*faceIndex] * faceFluidPressure[*faceIndex] - rhoP + PRhoGravity[kf]);

          length[kf].Normalize();
          flowVector[*faceIndex] = length[kf] * q[kf] / std::max(m_rho_o, faceFluidDensity[kf]) * 0.5 / w; //faceArea[*faceIndex];

        }
      }
    }
  }


  {
    std::map<PhysicalDomainT::ObjectDataStructureKeys, sArray1d> syncedFields_tmp;
    syncedFields_tmp[PhysicalDomainT::FiniteElementFaceManager].push_back("flowVector");

    partition.SynchronizeFields( syncedFields_tmp, CommRegistry::parallelPlateFlowSolver);
  }


  // Get field values of child faces from their parents.
  for( localIndex kf=0 ; kf<domain.m_feFaceManager.DataLengths() ; ++kf )
  {
    if (flowFaceType[kf] == 1)
    {
      for (localIndex i = 0; i < domain.m_feFaceManager.m_childIndices[kf].size(); ++i)
      {
        faceFluidPressure[domain.m_feFaceManager.m_childIndices[kf][i]] = faceFluidPressure[kf];
        flowVector[domain.m_feFaceManager.m_childIndices[kf][i]] = flowVector[kf];
        aperture[domain.m_feFaceManager.m_childIndices[kf][i]] = aperture[kf];
      }
    }
  }

  //Calculate nodal pressure
  iArray1d nFace(domain.m_feNodeManager.DataLengths());
  nFace = 0;
  p = 0.0;
  for( localIndex kf=0 ; kf<domain.m_feFaceManager.DataLengths() ; ++kf )
  {
    if (flowFaceType[kf] == 1)
    {
      for( lArray1d::iterator j = domain.m_feFaceManager.m_toNodesRelation[kf].begin() ;
          j!=domain.m_feFaceManager.m_toNodesRelation[kf].end() ; ++j )
      {
        p[*j] += faceFluidPressure[kf];
        nFace[*j]++;
      }
    }
  }

  for (localIndex i = 0; i<domain.m_feNodeManager.DataLengths(); ++i)
  {
    if (nFace[i]>0)
    {
      p[i] /= nFace[i];
    }
  }

  // Get nodal pressure of child nodes from their parents.
  for (localIndex i = 0; i<domain.m_feNodeManager.DataLengths(); ++i)
  {
    localIndex ancestor = i;
    while (domain.m_feNodeManager.m_parentIndex[ancestor] < domain.m_feNodeManager.DataLengths())
    {
      ancestor = domain.m_feNodeManager.m_parentIndex[ancestor];
    }
    if (ancestor != i)
    {
      p[i] = p[ancestor];
    }
  }



  {
    std::map<PhysicalDomainT::ObjectDataStructureKeys, sArray1d> syncedFields_tmp;
    syncedFields_tmp[PhysicalDomainT::FiniteElementNodeManager].push_back("nodalPressure");

    partition.SynchronizeFields( syncedFields_tmp, CommRegistry::parallelPlateFlowSolver);
  }



}


void ParallelPlateFlowSolverExplicit::GenerateParallelPlateGeometricQuantities( PhysicalDomainT& domain,
        realT time,realT dt )
{
  const iArray1d& flowFaceType = domain.m_feFaceManager.GetFieldData<int>("flowFaceType");
  //  const iArray1d& flowEdgeType = domain.m_edgeManager.GetFieldData<int>("flowEdgeType");
  rArray1d& fluidVolume  = domain.m_feFaceManager.GetFieldData<FieldInfo::volume>();
  rArray1d& faceArea = domain.m_feFaceManager.GetFieldData<realT>("faceArea");
  const Array1dT<R1Tensor>& faceNormal = domain.m_feFaceManager.GetFieldData<R1Tensor>("faceNormal0");
  rArray1d& effectiveStressN = domain.m_feFaceManager.GetFieldData<realT>("effectiveStressN");
  rArray1d& faceFluidPressure = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();
  rArray1d* stressNOnFace = domain.m_feFaceManager.GetFieldDataPointer<realT>("stressNOnFace");

  rArray1d* apertureOffset = domain.m_feFaceManager.GetFieldDataPointer<realT>("apertureOffset");
  rArray1d& faceFluidDensity = domain.m_feFaceManager.GetFieldData<FieldInfo::density>();

  if (m_updateFaceArea)
  {
    for( localIndex kf=0 ; kf<domain.m_feFaceManager.DataLengths() ; ++kf )
    {
      if( flowFaceType[kf] == 1 )
      {
        faceArea[kf] = domain.m_feFaceManager.SurfaceArea( domain.m_feNodeManager, kf );
      }
    }
  }

  rArray1d& aperture = domain.m_feFaceManager.GetFieldData<realT>( ApertureStr );
  Array1dT<R1Tensor>& faceCenter = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr );

  for( localIndex kf=0 ; kf<domain.m_feFaceManager.DataLengths() ; ++kf )
  {
    if( flowFaceType[kf] == 1 )
    {
      localIndex numChildren = domain.m_feFaceManager.m_childIndices[kf].size();
      realT aperture0 = aperture[kf];

      if (numChildren >= 1)  // This is actually an open fracture
      {
        domain.m_feFaceManager.FaceCenter( domain.m_feNodeManager, kf, faceCenter[kf] );
        R1Tensor gap;
        R1Tensor N;

        if (numChildren <= 1)
        {
          N = faceNormal( kf );
        }
        else
        {
          N = faceNormal( domain.m_feFaceManager.m_childIndices[kf][0] );
        }
        gap = domain.m_feFaceManager.CalculateGapVector( domain.m_feNodeManager, kf );
        aperture[kf] = Dot(gap,N) ;
      }
      else if (!(m_flowFaceSetName.empty()) && numChildren == 0)  // This is an embedded flow path.
      {
        domain.m_feFaceManager.CalculateStressOnFace(domain.m_feElementManager, domain.m_feNodeManager, kf);
        effectiveStressN[kf] = std::max(0.0 , -(*stressNOnFace)[kf] - faceFluidPressure[kf]);

      }

      if (false)
      {
        if( aperture[kf]<m_min_aperture )
          aperture[kf] = m_min_aperture;
        // PFU: I commented this out because we only need to put an upper limit on the permeability aperture, not the storage aperture.
        //      else if( aperture[kf] > m_max_aperture )
        //        aperture[kf] = m_max_aperture;
      }
      else
      {
        if (aperture[kf] < 0.0 || (!(m_flowFaceSetName.empty()) && numChildren == 0) )
        {
          if (m_bBarton != 0.0)
          {
            aperture[kf] = m_wZeroStress - m_aBarton * effectiveStressN[kf] / (1 + m_bBarton * effectiveStressN[kf]);
          }
          else
          {
            aperture[kf] = m_min_aperture;
          }
        }
        else
        {
          if( numChildren > 0){ // don't set aperture of unsplit faces
            aperture[kf] += m_wZeroStress;
          }
        }
      }

      if( apertureOffset != NULL )
      {
        aperture[kf] += (*apertureOffset)[kf];
      }

      if (aperture0 > m_min_aperture)  //Otherwise it is the first time step and initial aperture is zero
      {
        aperture[kf] = aperture0 * m_apertureMovingAverageCoeff + aperture[kf] * (1.0 - m_apertureMovingAverageCoeff);
      }

      fluidVolume[kf] = aperture[kf] * faceArea[kf];

      if( isZero( faceFluidDensity[kf] ) )
      {
        faceFluidDensity[kf] =  this->m_rho_o;
      }


    }
  }




  const iArray1d& flowEdgeType = domain.m_feEdgeManager.GetFieldData<int>("flowEdgeType");
  rArray1d& edgeLength = domain.m_feEdgeManager.GetFieldData<realT>("length");
  Array1dT<R1Tensor>& edgeCenter = domain.m_feEdgeManager.GetFieldData<R1Tensor>("center");

  for( localIndex ke=0 ; ke<domain.m_feEdgeManager.DataLengths() ; ++ke )
  {
    // only do work on edges that are active flow edges
    if( flowEdgeType[ke]== 1 )
    {
      domain.m_feEdgeManager.EdgeCenter( domain.m_feNodeManager, ke, edgeCenter[ke] );
      edgeLength[ke] = domain.m_feEdgeManager.EdgeLength( domain.m_feNodeManager, ke );
    }
  }


}


void ParallelPlateFlowSolverExplicit::CalculateAndApplyMassFlux( const realT dt, PhysicalDomainT& domain )
{

  const rArray1d& aperture = domain.m_feFaceManager.GetFieldData<realT>( ApertureStr );
  const Array1dT<R1Tensor>& faceCenter = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr );

  rArray1d& edgePermeability = domain.m_feEdgeManager.GetFieldData<realT>(PermeabilityStr);
  rArray1d& faceFluidPressure = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();
  rArray1d& faceFluidDensity = domain.m_feFaceManager.GetFieldData<FieldInfo::density>();
  rArray1d& faceFluidMass = domain.m_feFaceManager.GetFieldData<FieldInfo::mass>();
  rArray1d& flowRate = domain.m_feFaceManager.GetFieldData<realT>("elementFlowRate");
//  iArray1d& flowFaceType = domain.m_feFaceManager.GetFieldData<int>("flowFaceType");
  rArray1d& facePerm = domain.m_feFaceManager.GetFieldData<realT>("permeability");
  facePerm = -1.0;

  iArray1d& tracer = domain.m_feFaceManager.GetFieldData<int>( "tracer" );


  rArray1d& massRate = domain.m_feEdgeManager.GetFieldData<realT>("massRate");


  const rArray1d& edgeLength = domain.m_feEdgeManager.GetFieldData<realT>("length");
  const Array1dT<R1Tensor>& edgeCenter = domain.m_feEdgeManager.GetFieldData<R1Tensor>("center");

  //  const iArray1d& flowFaceType = domain.m_feFaceManager.GetFieldData<int>("flowFaceType");
  const iArray1d& flowEdgeType = domain.m_feEdgeManager.GetFieldData<int>("flowEdgeType");


  const Array1dT<lSet>& edgeToFlowFaces = domain.m_feEdgeManager.GetUnorderedVariableOneToManyMap("edgeToFlowFaces");

  R1Tensor la, lb;

  m_stabledt.m_maxdt = std::numeric_limits<double>::max();



  for( localIndex kf=0 ; kf<domain.m_feFaceManager.DataLengths() ; ++kf )
  {
    flowRate[kf] = 0.0;
  }


  // loop over all edges
  for( localIndex ke=0 ; ke<domain.m_feEdgeManager.DataLengths() ; ++ke )
  {

    // only do work on edges that are active flow edges
    if( flowEdgeType[ke]== 1 )
    {
      realT w = edgeLength[ke];


      // we have two different ways of calculating the flow. One for edges attached to only 2
      // faces, and one for more than two faces.
      const unsigned int numFlowFaces = edgeToFlowFaces[ke].size();

      if( numFlowFaces == 2 )
      {
        lSet::const_iterator edgeToFlowFace = edgeToFlowFaces[ke].begin();
        const localIndex faceIndex0 = *edgeToFlowFace;
        const localIndex faceIndex1 = *(++edgeToFlowFace);


        la = edgeCenter[ke];
        la -= faceCenter[faceIndex0];

        lb = edgeCenter[ke];
        lb -= faceCenter[faceIndex1];

        realT norm_la = la.L2_Norm();
        realT norm_lb = lb.L2_Norm();

        realT permeabilityAperture0 = std::min(aperture[faceIndex0], m_max_aperture);
        realT permeabilityAperture1 = std::min(aperture[faceIndex1], m_max_aperture);

        edgePermeability[ke] = CalculatePermeability( norm_la , norm_lb,
                                                      permeabilityAperture0, permeabilityAperture1,
                                                      w,m_mu,m_SHP_FCT);

        realT PRhoGravity = CalculatePRhoGravity(faceCenter[faceIndex0], faceCenter[faceIndex1],
                                                    faceFluidDensity[faceIndex0], faceFluidDensity[faceIndex1],
                                                    m_gravityVector);
         // determine the mass flux across edge at t_n+1/2
        if( facePerm[faceIndex0] < 0.0 )
        {
          facePerm[faceIndex0] = faceFluidDensity[faceIndex0] * CalculatePermeability( norm_la, permeabilityAperture0, w, m_mu, m_SHP_FCT );
        }
        if( facePerm[faceIndex1] < 0.0 )
        {
          facePerm[faceIndex1] = faceFluidDensity[faceIndex1] * CalculatePermeability( norm_lb, permeabilityAperture1, w, m_mu, m_SHP_FCT );
        }

        // determine the mass flux across edge at t_n+1/2
        massRate[ke] = edgePermeability[ke] * ( -( faceFluidDensity[faceIndex0] * faceFluidPressure[faceIndex0])
            + ( faceFluidDensity[faceIndex1] * faceFluidPressure[faceIndex1] )
            - PRhoGravity );

        if (domain.m_feFaceManager.m_toEdgesRelation[faceIndex0].size() == 2 )
        {//2D
          if (domain.m_feFaceManager.m_toEdgesRelation[faceIndex0][0] == ke)
          {
            flowRate[faceIndex0] += 0.5*massRate[ke] / m_rho_o;
          }
          else
          {
            flowRate[faceIndex0] -= 0.5*massRate[ke] / m_rho_o;
          }

          if (domain.m_feFaceManager.m_toEdgesRelation[faceIndex1][0] == ke)
          {
            flowRate[faceIndex1] -= 0.5*massRate[ke] / m_rho_o;
          }
          else
          {
            flowRate[faceIndex1] += 0.5*massRate[ke] / m_rho_o;
          }
        }
        else
        {//3D
          flowRate[faceIndex0] += 0.5 * fabs(massRate[ke]) / std::max(m_rho_o, faceFluidDensity[faceIndex0]);
          flowRate[faceIndex1] += 0.5 * fabs(massRate[ke]) / std::max(m_rho_o, faceFluidDensity[faceIndex1]);
        }


        faceFluidMass[faceIndex0] +=  massRate[ke] * dt;
        faceFluidMass[faceIndex1] -=  massRate[ke] * dt;

        if( tracer[faceIndex0] == 1 )
        {
          tracer[faceIndex1] = 1;
        }
        if( tracer[faceIndex1] == 1 )
        {
          tracer[faceIndex0] = 1;
        }


        realT thisdt = 6 * m_mu *( norm_la + norm_lb ) * ( norm_la + norm_lb )
                     / ( m_bulk_modulus * 0.25 * (permeabilityAperture0+permeabilityAperture1)* (permeabilityAperture0+permeabilityAperture1) ) ;

        if( thisdt < m_stabledt.m_maxdt )
          m_stabledt.m_maxdt = thisdt;



      }
      else if( numFlowFaces > 2 )
      {
        realT rhoP = 0.0;
        realT sumK = 0.0;
        Array1dT<R1Tensor> length(numFlowFaces);
        rArray1d k(numFlowFaces);
        rArray1d q(numFlowFaces);
        rArray1d kRhoP(numFlowFaces);
        rArray1d PRhoGravity(numFlowFaces);

        lSet::const_iterator faceIndex=edgeToFlowFaces[ke].begin();

        for( localIndex kf=0 ; kf<numFlowFaces ; ++kf, ++faceIndex)
        {

          length[kf] = edgeCenter[ke];
          length[kf] -= faceCenter[*faceIndex];

          realT permeabilityAperture = std::min(aperture[*faceIndex], m_max_aperture);

          k[kf] = CalculatePermeability( length[kf].L2_Norm(),
                                         permeabilityAperture,
                                         w, m_mu, m_SHP_FCT );

          PRhoGravity[kf] = CalculatePRhoGravity( faceCenter[*faceIndex], edgeCenter[ke], faceFluidDensity[*faceIndex], m_gravityVector);

          sumK += k[kf];

          kRhoP[kf] = k[kf] * faceFluidDensity[*faceIndex] * faceFluidPressure[*faceIndex] + k[kf] * PRhoGravity[kf];
          rhoP += kRhoP[kf];



          realT thisdt = 6 * m_mu * 4 * Dot(length[kf],length[kf])
          / ( m_bulk_modulus * permeabilityAperture * permeabilityAperture ) ;

          if( thisdt < m_stabledt.m_maxdt )
            m_stabledt.m_maxdt = thisdt;
        }


        rhoP /= sumK ;
        faceIndex=edgeToFlowFaces[ke].begin();

        int tracerActive = 0;
        for( localIndex kf=0 ; kf<numFlowFaces ; ++kf, ++faceIndex )
        {
          q[kf] = k[kf] * ( faceFluidDensity[*faceIndex] * faceFluidPressure[*faceIndex] - rhoP + PRhoGravity[kf]);

          faceFluidMass[*faceIndex] -= q[kf] * dt;

          if (domain.m_feFaceManager.m_toEdgesRelation[*faceIndex].size() == 2 )
          {//2D
            if (domain.m_feFaceManager.m_toEdgesRelation[*faceIndex][0] == ke)
            {
              flowRate[*faceIndex] -= 0.5*q[kf] / m_rho_o;
            }
            else
            {
              flowRate[*faceIndex] += 0.5*q[kf] / m_rho_o;
            }
          }
          else //3D
          {
            flowRate[*faceIndex] += 0.5* fabs(q[kf]) / std::max(m_rho_o, faceFluidDensity[*faceIndex]);
          }




          if( tracer[*faceIndex] == 1 )
          {
            tracerActive = 1;
          }

        }

        for( faceIndex=edgeToFlowFaces[ke].begin() ; faceIndex!=edgeToFlowFaces[ke].end() ; ++faceIndex )
        {
          if( tracerActive == 1 && tracer[*faceIndex] == 0 )
          {
            tracer[*faceIndex] = 1;
          }
        }

      }
    }
  }



  if (m_stabledt.m_maxdt != std::numeric_limits<double>::max())
  {
    m_stabledt.m_maxdt *= this->m_courant;
  }
}

void ParallelPlateFlowSolverExplicit::CalculateCarterLeakOff( const realT time,
                                                      const realT dt,
                                                      PhysicalDomainT& domain )
{
  rArray1d* initialSaturatedTime = domain.m_feFaceManager.GetFieldDataPointer<realT>("initialSaturatedTime");
  rArray1d* totalLeakedVolume = domain.m_feFaceManager.GetFieldDataPointer<realT>("totalLeakedVolume");
  iArray1d& flowFaceType = domain.m_feFaceManager.GetFieldData<int>("flowFaceType");
  rArray1d& faceFluidPressure = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();
  rArray1d& faceFluidMass = domain.m_feFaceManager.GetFieldData<FieldInfo::mass>();
  rArray1d& faceArea = domain.m_feFaceManager.GetFieldData<realT>("faceArea");



  if (totalLeakedVolume != NULL)
  {
    for( localIndex kf=0 ; kf<domain.m_feFaceManager.DataLengths() ; ++kf )
    {
      if (flowFaceType[kf] == 1 &&
          (*initialSaturatedTime)[kf] != std::numeric_limits<realT>::max() &&
          (faceFluidPressure[kf] > 0.0 || m_overLeakCompensation) )
      {
        localIndex face0 = kf;
        localIndex face1 = kf;
        if (domain.m_feFaceManager.m_childIndices[kf].size() == 0 && !(m_flowFaceSetName.empty())) ////Embedded flow plane
        {
          face0 = kf;
          face1 = kf;
        }
        else
        {
          if (domain.m_feFaceManager.m_toEdgesRelation[kf].size() == 2 ) //2D
          {
            face0 = domain.m_feFaceManager.m_childIndices[kf][0];
            face1 = domain.m_feFaceManager.m_childIndices[kf][1];
          }
          else if (domain.m_feFaceManager.m_toEdgesRelation[kf].size() > 2 ) //3D open fracture
          {
            face0 = kf;
            face1 = domain.m_feFaceManager.m_childIndices[kf][0];
          }

        }

        realT leakOffVelocity = m_leakoffCoef / sqrt( time + dt * 0.5 - (*initialSaturatedTime)[kf]);
        if (m_pressureDependentLeakoff == 1)
        {
          leakOffVelocity *= faceFluidPressure[kf] - m_farFieldPorePressure;
        }
        realT leakOffMassInc = leakOffVelocity * dt *  faceArea[kf] * m_rho_o;
        leakOffMassInc *= 2.0; //Counting for both sides.

        if (leakOffMassInc < faceFluidMass[kf] || m_overLeakCompensation)
        {
          faceFluidMass[kf] -= leakOffMassInc;
          (*totalLeakedVolume)[face0] += leakOffVelocity * dt * faceArea[kf];

          if (domain.m_feFaceManager.m_toEdgesRelation[kf].size() > 0) (*totalLeakedVolume)[face1] += leakOffVelocity * dt * faceArea[kf];

        }
        else
        {
          (*totalLeakedVolume)[face0] += faceFluidMass[kf] * 0.5 / m_rho_o;
          if (domain.m_feFaceManager.m_toEdgesRelation[kf].size() > 0) (*totalLeakedVolume)[face1] += faceFluidMass[kf] * 0.5 / m_rho_o;
          faceFluidMass[kf] = 0.0;
        }


      }
    }
  }


}

void ParallelPlateFlowSolverExplicit::CalculateMatrixFlowLeakOff( const realT time,
                                                      const realT dt,
                                                      PhysicalDomainT& domain )
{
  iArray1d& flowFaceType = domain.m_feFaceManager.GetFieldData<int>("flowFaceType");
  rArray1d& faceToMatrixLeakOffRate = domain.m_feFaceManager.GetFieldData<realT>("faceToMatrixLeakOffRate");
  rArray1d& faceFluidMass = domain.m_feFaceManager.GetFieldData<FieldInfo::mass>();
  rArray1d& faceArea = domain.m_feFaceManager.GetFieldData<realT>("faceArea");


  for( localIndex kf=0 ; kf<domain.m_feFaceManager.DataLengths() ; ++kf )
  {
    if (flowFaceType[kf] == 1)
    {
      faceFluidMass[kf] -= faceToMatrixLeakOffRate[kf] * dt * m_rho_o * faceArea[kf];
      faceFluidMass[kf] = std::max(0.0, faceFluidMass[kf]);
    }
  }

}




// We know the combined flow rate and have to distribute the given flow rate among the faces in the set.
void ParallelPlateFlowSolverExplicit::ApplyFluxBoundaryCondition( const realT time,
                                                          const realT dt,
                                                          const int cycleNumber,
                                                          const int rank,
                                                          PhysicalDomainT& domain )
{
  m_dT = dt;

  rArray1d& faceFluidMass = domain.m_feFaceManager.GetFieldData<FieldInfo::mass>();
  const iArray1d& flowFaceType = domain.m_feFaceManager.GetFieldData<int>("flowFaceType");
  //const rArray1d& edgeLength = domain.m_feEdgeManager.GetFieldData<realT>("length");
  const Array1dT<R1Tensor>& edgeCenter = domain.m_feEdgeManager.GetFieldData<R1Tensor>("center");
  const Array1dT<R1Tensor>& faceCenter = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr );
  const rArray1d& aperture = domain.m_feFaceManager.GetFieldData<realT>( ApertureStr );
  const rArray1d& faceFluidPressure = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();
  const rArray1d& faceFluidVolume = domain.m_feFaceManager.GetFieldData<FieldInfo::volume>();
  const rArray1d& faceFluidDensity = domain.m_feFaceManager.GetFieldData<FieldInfo::density>();
  const iArray1d& isGhost = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();
  rArray1d& flowRate = domain.m_feFaceManager.GetFieldData<realT>("elementFlowRate");
  rArray1d& flowRateBC = domain.m_feFaceManager.GetFieldData<realT>("flowRateBC");
  flowRateBC = 0.0;

  for( Array1dT<BoundaryConditionBase*>::const_iterator bcItr=domain.m_feFaceManager.m_bcData.begin() ; bcItr!=domain.m_feFaceManager.m_bcData.end() ; ++ bcItr )
  {
    // check to see if the requested field has a boundary condition applied to it.
    BoundaryConditionBase* bc = *bcItr;
#ifdef SRC_EXTERNAL
    if( streq( bc->GetBoundaryConditionName(), PerforatedCasedWellboreBoundaryCondition::BoundaryConditionName()) )
    {
      PerforatedCasedWellboreBoundaryCondition *cpbc = static_cast<PerforatedCasedWellboreBoundaryCondition *>(bc);
      cpbc->Apply(domain, time, dt);
    }
    else if( streq( bc->GetBoundaryConditionName(), CavityPressureBoundaryCondition::BoundaryConditionName()) )
    {
      CavityPressureBoundaryCondition *cpbc = static_cast<CavityPressureBoundaryCondition *>(bc);
      realT dmass = 0., dmassAll;
      realT dmassc = 0;
      const realT pressure = cpbc->GetPressure();
      for(sArray1d::const_iterator its=cpbc->m_setNamesHydro.begin(); its!=cpbc->m_setNamesHydro.end(); ++its)
      {
        std::map< std::string, lSet >::const_iterator setMap = domain.m_feFaceManager.m_Sets.find( *its );
        if( setMap != domain.m_feFaceManager.m_Sets.end() )
        {
          const lSet& set = setMap->second;
          for( lSet::const_iterator a=set.begin() ; a!=set.end() ; ++a)
          {
            const localIndex aa = *a;
            if (flowFaceType[aa] == 1 && isGhost[aa]<0)
            {
              if (0)
              {
                // Assume pressure drop between cavity and fracture
                realT k = 0;
                {
                  const localIndex ke = domain.m_feFaceManager.m_toEdgesRelation[aa][0];
                  R1Tensor length(edgeCenter[ke]);
                  length -= faceCenter[aa];
                  const realT w = 1.0;
                  const realT permeabilityAperture = std::min(aperture[aa], m_max_aperture);
                  k = CalculatePermeability( length.L2_Norm(),
                                             permeabilityAperture,
                                             w, m_mu, m_SHP_FCT );
                }
                const realT q = k * (pressure - faceFluidPressure[aa]);
                dmassc = q * dt * std::max(faceFluidDensity[aa], m_rho_o);
              }
              else 
              {
                // No pressure drop between cavity and fracture
                dmassc = faceFluidVolume[aa]*m_rho_o*(pressure - faceFluidPressure[aa])/m_bulk_modulus;
                // dmassc = std::max(0.0, dmassc);
              }
              faceFluidMass[aa] += dmassc;
              dmass += dmassc;
            }
          }
        }
      }
      MPI_Allreduce(&dmass, &dmassAll, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      cpbc->UpdateMass(dmassAll, dt, time);
      //std::cout <<"Cavity BC: dmassAll=" << dmassAll << ", t=" << time << ", P_inj=" << pressure <<" mass=" << ret[0] << " volume=" << ret[1] << " Kf=" << ret[2] << " rho0=" << ret[3] << std::endl;
    }
    else
#endif
      if( streq( bc->GetFieldName(time), "combinedFlowRate") )
      {
        realT qTotal = 0.0;
        realT sumK_KP[] = {0.0, 0.0};
        rArray1d k;
        //This actually should be a map, since the same face can be in different facesets.  But I am being lazy...
        int nInlets = 0;
        realT pInj;

        localIndex kf=0;

        for(localIndex i =0; i < bc->m_setNames.size(); ++i)
        {
          std::map< std::string, lSet >::iterator setMap = domain.m_feFaceManager.m_Sets.find( bc->m_setNames[i] );
          if( setMap != domain.m_feFaceManager.m_Sets.end() )
          {
            lSet& set = setMap->second;

            lSet::const_iterator b=set.begin();

            qTotal = bc->GetValue(domain.m_feFaceManager, b, time); // The iterator b in this function is not doing anything.


            R1Tensor length;
            realT kP;

            for( lSet::const_iterator a=set.begin() ; a!=set.end() ; ++a)
            {
              if (flowFaceType[*a] == 1 )
              {
                localIndex ke = domain.m_feFaceManager.m_toEdgesRelation[*a][0];
                length = edgeCenter[ke];
                length -= faceCenter[*a];
                realT w = length.L2_Norm() * 2;
                realT permeabilityAperture = std::min(aperture[*a], m_max_aperture);
                k.push_back( CalculatePermeability( length.L2_Norm(),
                                                    permeabilityAperture,
                                                    w, m_mu, m_SHP_FCT ));

                if (isGhost[*a]<0)
                {
                  nInlets++;
                  sumK_KP[0] += k[kf];

                  kP = k[kf] * faceFluidPressure[*a];
                  sumK_KP[1] += kP;
                }
                ++kf;

              }
            }
          }
        }

        int myNInlets = nInlets;
        realT mySumK_KP[] = {sumK_KP[0], sumK_KP[1]};
        int myRank0 = nInlets * rank;
        int rank0;

        MPI_Allreduce(&myNInlets, &nInlets, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&mySumK_KP, &sumK_KP, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        if (cycleNumber%100 == 0) MPI_Allreduce(&myRank0, &rank0, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD);

        if (nInlets >=1)
        {
          sumK_KP[1] += qTotal;
          pInj = sumK_KP[1] / sumK_KP[0];

          kf = 0;


          for(localIndex i =0; i < bc->m_setNames.size(); ++i)
          {
            std::map< std::string, lSet >::iterator setMap = domain.m_feFaceManager.m_Sets.find( bc->m_setNames[i] );
            if( setMap != domain.m_feFaceManager.m_Sets.end() )
            {
              lSet& set = setMap->second;

              lSet::const_iterator b=set.begin();

              for( lSet::const_iterator a=set.begin() ; a!=set.end() ; ++a)
              {
                if (flowFaceType[*a] == 1 )
                {
                  realT q = k[kf] * (pInj - faceFluidPressure[*a]);
                  faceFluidMass[*a] += q * dt * std::max(faceFluidDensity[*a], m_rho_o);
                  flowRate[*a] += 0.5* fabs(q) ;
                  flowRateBC[*a] += q;
                  ++kf;
                }
              }
            }
          }

          if (cycleNumber%100 == 0 && myRank0 == rank0)
          {
            std::cout <<"Flow rate BC applied to set: " << bc->m_setNames[0] << " etc. , t=" << time << ", P_inj=" << pInj <<" rank: " <<rank << std::endl;
          }
        }
      }
  }


  //  BoundaryConditionFunctions::ApplyCombinedFluxBoundaryCondition<realT>(domain.m_feFaceManager,
  //                                                                        "appliedFlowRateBC",time);
  //
  //
  //  for( localIndex kf=0 ; kf<domain.m_feFaceManager.DataLengths() ; ++kf )
  //  {
  //    if( flowFaceType[kf] == 1 )
  //    {
  //      faceFluidMass[kf] += appliedFlowRateBC[kf] * dt * m_rho_o;
  //    }
  //  }

  // Edge fluxes
  BoundaryConditionFunctions::ApplyBoundaryCondition<realT>(this, &ParallelPlateFlowSolverExplicit::FlowControlledBoundaryCondition,
                                domain, domain.m_feEdgeManager, "FixedFlowRate", time );

}


// constant flow into faces
void ParallelPlateFlowSolverExplicit::FlowControlledBoundaryCondition( PhysicalDomainT& domain,
                                                               ObjectDataStructureBaseT& object ,
                                                               BoundaryConditionBase* bc ,
                                                               const lSet& set,
                                                               realT time )
{


	//iArray1d& edge_is_ghost = domain.m_feEdgeManager.GetFieldData<FieldInfo::ghostRank>();


	rArray1d& faceFluidMass = domain.m_feFaceManager.GetFieldData<FieldInfo::mass>();
	const iArray1d& flowFaceType = domain.m_feFaceManager.GetFieldData<int>("flowFaceType");
	const rArray1d& edgeLengths = domain.m_feEdgeManager.GetFieldData<realT>("length");
	//const Array1dT<R1Tensor>& edgeCenter = domain.m_feEdgeManager.GetFieldData<R1Tensor>("center");
	//const Array1dT<R1Tensor>& faceCenter = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr );
	const rArray1d& apertures = domain.m_feFaceManager.GetFieldData<realT>( ApertureStr );
	const rArray1d& faceFluidDensity = domain.m_feFaceManager.GetFieldData<FieldInfo::density>();
	//const iArray1d& isGhost = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();

	const Array1dT<lSet>& edgeToFlowFaces = domain.m_feEdgeManager.GetUnorderedVariableOneToManyMap("edgeToFlowFaces");

	/*Epetra_IntSerialDenseVector  face_dof(1);
	Epetra_SerialDenseVector     face_rhs(1);
	Epetra_SerialDenseMatrix     face_matrix(1,1);*/

	// loop over edges, apply flux to faces

	for( lSet::const_iterator eg=set.begin(); eg != set.end(); ++eg  ){

		realT qdt = bc->GetValue(domain.m_feEdgeManager,eg,time)*m_dT;

		for( lSet::const_iterator fc=edgeToFlowFaces[*eg].begin() ; fc!=edgeToFlowFaces[*eg].end() ; ++fc )
		{
		    if(flowFaceType[*fc] == 1){
			    const realT area = edgeLengths[*eg]*BoundedAperture(apertures[*fc]);

			    const realT massFlux = qdt*area*std::max(faceFluidDensity[*fc], m_rho_o);
			    faceFluidMass[*fc] += massFlux;
		    }
		}
	}
}


void ParallelPlateFlowSolverExplicit::UpdateEOS( const realT time,
                                         const realT dt ,
                                         PhysicalDomainT& domain )
{
  rArray1d& faceFluidPressure = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();
  rArray1d& faceFluidDensity = domain.m_feFaceManager.GetFieldData<FieldInfo::density>();
  rArray1d& faceFluidMass = domain.m_feFaceManager.GetFieldData<FieldInfo::mass>();
  const rArray1d& fluidVolume  = domain.m_feFaceManager.GetFieldData<FieldInfo::volume>();

  const iArray1d& flowFaceType = domain.m_feFaceManager.GetFieldData<int>("flowFaceType");

  rArray1d* initialSaturatedTime = domain.m_feFaceManager.GetFieldDataPointer<realT>("initialSaturatedTime");

  int rank, size ;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  for( localIndex kf=0 ; kf<domain.m_feFaceManager.DataLengths() ; ++kf )
  {
    if( flowFaceType[kf] == 1 )
    {
      faceFluidDensity[kf] = faceFluidMass[kf] / fluidVolume[kf];
    }
  }

  BoundaryConditionFunctions::ApplyDirichletBoundaryCondition<realT>(domain.m_feFaceManager,
                                                                     Field<FieldInfo::density>::Name(),time);


  //  double mass_input = 0;
  //  double mass_output = 0;

  //  lSet& source = domain.m_faceManager.m_Sets["FF0"];
  //  lSet& sink = domain.m_faceManager.m_Sets["FF1"];

  /*
  for( lSet::const_iterator k=source.begin() ; k!=source.end() ; ++k )
  {
    const double oldMass = faceFluidMass[*k];
    const double newMass = faceFluidDensity[*k] * fluidVolume[*k];
    mass_input += newMass - oldMass;
  }

  for( lSet::const_iterator k=sink.begin() ; k!=sink.end() ; ++k )
  {
    const double oldMass = faceFluidMass[*k];
    const double newMass = faceFluidDensity[*k] * fluidVolume[*k];
    mass_output += newMass - oldMass;
  }
   */

  for( localIndex kf=0 ; kf<domain.m_feFaceManager.DataLengths() ; ++kf )
  {
    if( flowFaceType[kf] == 1 )
    {
      faceFluidMass[kf] = faceFluidDensity[kf] * fluidVolume[kf];

      //      faceFluidPressure[kf] = m_bulk_modulus * ( 1.0 - 1.0 / faceFluidDensity[kf] );
      // faceFluidPressure[kf] = m_bulk_modulus * ( faceFluidDensity[kf] - 1.0 ); // wrong - does not account for units of rho_o
      faceFluidPressure[kf] =P_EOS(faceFluidDensity[kf], m_bulk_modulus, m_rho_o, m_pressureCap);
      if( faceFluidPressure[kf] < 0.0 ) faceFluidPressure[kf] = 0.0;


    }
  }

  BoundaryConditionFunctions::ApplyDirichletBoundaryCondition<realT>(domain.m_feFaceManager,
                                                                     Field<FieldInfo::pressure>::Name(),time);

  for( localIndex kf=0 ; kf<domain.m_feFaceManager.DataLengths() ; ++kf )
  {
    if( faceFluidPressure[kf] < 0.0 ) faceFluidPressure[kf] = 0.0;

    if (flowFaceType[kf] < 1) faceFluidPressure[kf] = 0.0;

    if( flowFaceType[kf] == 1 && faceFluidPressure[kf] > 0.0)
    {
      faceFluidDensity[kf] = Inverse_EOS(faceFluidPressure[kf], m_bulk_modulus, m_rho_o, m_pressureCap);
      // faceFluidMass[kf] = faceFluidDensity[kf] * fluidVolume[kf];
      // We face a dilemma here.  We  skip the mass update above so that the negative mass in this cell attached to a pressure boundary indicate the amount of fluid that has flow into the system.
      // If we need to use the switch boundary condition to initialize a proper amount of fluid in natural fractures, we have to use the density bc instead of the pressure bc.
    }
  }

  //Simple Carter's leakoff model; initializing face
  if (initialSaturatedTime != NULL)
  {
    for( localIndex kf=0 ; kf<domain.m_feFaceManager.DataLengths() ; ++kf )
    {
      if (flowFaceType[kf] == 1 && faceFluidPressure[kf] > 0.0 && (*initialSaturatedTime)[kf] == std::numeric_limits<realT>::max())
      {
        (*initialSaturatedTime)[kf] = time - dt/2;
      }
    }
  }





  //  std::cout<<"  Mass In-Out = Net: "<<mass_input<<" - "<<-mass_output<<" = "<<mass_input+mass_output<<std::endl;
}


//Calculating (through a simple averaging) nodal pressure from face pressure for the coupling with the porous media solver.
void ParallelPlateFlowSolverExplicit::CalculateNodalPressure ( PhysicalDomainT& domain, SpatialPartition& partition)
{
//  const iArray1d& ghostRank = domain.m_feNodeManager.GetFieldData<FieldInfo::ghostRank>();
  rArray1d& faceFluidPressure = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();
  rArray1d& fractureNodalPressure = domain.m_feNodeManager.GetFieldData<realT>("fractureNodalPressure");
  iArray1d& isFractureNode = domain.m_feNodeManager.GetFieldData<int>("isFractureNode");
  iArray1d& flowFaceType = domain.m_feFaceManager.GetFieldData<int>("flowFaceType");

  isFractureNode = 0;
  fractureNodalPressure = 0.0;

  for( localIndex kf=0 ; kf<domain.m_feFaceManager.DataLengths() ; ++kf )
  {
    if( flowFaceType[kf] == 1 )
    {
      for (localIndex i = 0; i < domain.m_feFaceManager.m_toNodesRelation[kf].size(); ++i)
      {
        localIndex nd = domain.m_feNodeManager.GetParentIndex(domain.m_feFaceManager.m_toNodesRelation[kf][i]);
        isFractureNode[nd]++;
        fractureNodalPressure[nd] += faceFluidPressure[kf];
      }
    }
  }

  for( localIndex nd=0 ; nd<domain.m_feNodeManager.DataLengths() ; ++nd )
  {
    if( isFractureNode[nd] >= 1 )
    {
      fractureNodalPressure[nd] /= isFractureNode[nd];
    }
  }

  {
    std::map<PhysicalDomainT::ObjectDataStructureKeys, sArray1d> syncedFields_tmp;

    syncedFields_tmp[PhysicalDomainT::FiniteElementFaceManager].push_back("fractureNodalPressure");
    syncedFields_tmp[PhysicalDomainT::FiniteElementFaceManager].push_back("isFractureNode");
    partition.SynchronizeFields( syncedFields_tmp, CommRegistry::parallelPlateFlowSolver );
  }

}




/// Register solver in the solver factory
REGISTER_SOLVER( ParallelPlateFlowSolverExplicit )
