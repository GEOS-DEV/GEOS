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
 * @file HydroStaticParallelPlateFlowSolver.cpp
 * @author fu4
 * @date March 26, 2014
 */

#include "HydroStaticParallelPlateFlowSolver.h"
#include "SolverFactory.h"
#include <algorithm>
#include <numeric>

#include "BoundaryConditions/ApplyBoundaryConditions.h"
#include "PhysicsSolverStrings.h"
#include "BoundaryConditions/BoundaryConditions.h"
using namespace PS_STR;
using namespace PPFS;


HydroStaticParallelPlateFlowSolver::HydroStaticParallelPlateFlowSolver( const std::string& name,
                                                  ProblemManagerT* const pm ):
    ParallelPlateFlowSolverExplicit(name,pm)
{
}

HydroStaticParallelPlateFlowSolver::~HydroStaticParallelPlateFlowSolver()
{
}


void HydroStaticParallelPlateFlowSolver::ReadXML( TICPP::HierarchicalDataNode* const hdn  )
{
  ParallelPlateFlowSolverExplicit::ReadXML( hdn );
  m_cavityVolume = hdn->GetAttributeOrDefault<realT>("cavityVolume", 0.0);
  m_cavityMass = hdn->GetAttributeOrDefault<realT>("initialMass", m_rho_o * m_cavityVolume);
  m_boreholeSetNames = hdn->GetStringVector("boreholeSetNames");
}


void HydroStaticParallelPlateFlowSolver::PostProcess(PhysicalDomainT & domain,
                                          SpatialPartition& partition,
                                          const sArray1d& namesOfSolverRegions)
{


  rArray1d& aperture = domain.m_feFaceManager.GetFieldData<realT>( ApertureStr );
  rArray1d& faceFluidPressure = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();
  const iArray1d& flowFaceType = domain.m_feFaceManager.GetFieldData<int>("flowFaceType");

  // Get field values of child faces from their parents.
  for( localIndex kf=0 ; kf<domain.m_feFaceManager.DataLengths() ; ++kf )
  {
    if (flowFaceType[kf] == 1)
    {
      for (localIndex i = 0; i < domain.m_feFaceManager.m_childIndices[kf].size(); ++i)
      {
        faceFluidPressure[domain.m_feFaceManager.m_childIndices[kf][i]] = faceFluidPressure[kf];
        aperture[domain.m_feFaceManager.m_childIndices[kf][i]] = aperture[kf];
      }
    }
  }
}





void HydroStaticParallelPlateFlowSolver::CalculateAndApplyMassFlux( const realT dt, PhysicalDomainT& domain )
{

  m_stabledt.m_maxdt = std::numeric_limits<double>::max();
}





void HydroStaticParallelPlateFlowSolver::UpdateEOS( const realT time,
                                         const realT dt ,
                                         PhysicalDomainT& domain )
{
  rArray1d& faceFluidPressure = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();
  rArray1d& faceFluidDensity = domain.m_feFaceManager.GetFieldData<FieldInfo::density>();
  rArray1d& faceFluidMass = domain.m_feFaceManager.GetFieldData<FieldInfo::mass>();
  const rArray1d& fluidVolume  = domain.m_feFaceManager.GetFieldData<FieldInfo::volume>();

  const iArray1d& flowFaceType = domain.m_feFaceManager.GetFieldData<int>("flowFaceType");

  rArray1d* initialSaturatedTime = domain.m_feFaceManager.GetFieldDataPointer<realT>("initialSaturatedTime");
  const iArray1d& isGhost = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();


  int rank, size ;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  realT totalFluidMass(0.0), totalFluidVolume(0.0), partitionFluidMass(0.0), partitionFluidVolume(0.0);

  for( localIndex kf=0 ; kf<domain.m_feFaceManager.DataLengths() ; ++kf )
  {
    if( flowFaceType[kf] == 1 && isGhost[kf] < 0)
    {
      partitionFluidMass += faceFluidMass[kf];
      partitionFluidVolume += fluidVolume[kf];
    }
  }

  MPI_Allreduce(&partitionFluidMass, &totalFluidMass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&partitionFluidVolume, &totalFluidVolume, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  realT meanDensity = (totalFluidMass + m_cavityMass)/ (totalFluidVolume + m_cavityVolume);
  realT meanPressure = P_EOS(meanDensity, m_bulk_modulus, m_rho_o, m_pressureCap);

  for( Array1dT<BoundaryConditionBase*>::const_iterator bcItr=domain.m_feFaceManager.m_bcData.begin() ; bcItr!=domain.m_feFaceManager.m_bcData.end() ; ++ bcItr )
  {
    // check to see if the requested field has a boundary condition applied to it.
    BoundaryConditionBase* bc = *bcItr;
    if( streq( bc->GetFieldName(time), Field<FieldInfo::pressure>::Name()) )
    {
      for(localIndex i =0; i < bc->m_setNames.size(); ++i)
      {
        std::map< std::string, lSet >::iterator setMap = domain.m_feFaceManager.m_Sets.find( bc->m_setNames[i] );
        if( setMap != domain.m_feFaceManager.m_Sets.end() )
        {
          lSet& set = setMap->second;
          meanPressure = bc->GetValue(domain.m_feFaceManager,set.begin(),time);
         }
      }
    }
  }

  meanDensity = Inverse_EOS(meanPressure, m_bulk_modulus, m_rho_o, m_pressureCap);


  for( localIndex kf=0 ; kf<domain.m_feFaceManager.DataLengths() ; ++kf )
  {
    if( flowFaceType[kf] == 1 )
    {
      faceFluidPressure[kf] = meanPressure;
      faceFluidMass[kf] = meanDensity * fluidVolume[kf];
      faceFluidDensity[kf] = meanDensity;
    }
  }

  m_cavityMass = meanDensity * m_cavityVolume;


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

  if (m_boreholeSetNames.size() > 0)
  {
    ApplyBoreholePressure( time, dt, meanPressure, domain );
  }

  //  std::cout<<"  Mass In-Out = Net: "<<mass_input<<" - "<<-mass_output<<" = "<<mass_input+mass_output<<std::endl;
}


void HydroStaticParallelPlateFlowSolver::ApplyBoreholePressure( const realT time, const realT dt, realT pressure, PhysicalDomainT& domain )
{
  FaceManagerT& faceManager = domain.m_feFaceManager;
  NodeManagerT& nodeManager = domain.m_feNodeManager;
  Array1dT<R1Tensor>& force = nodeManager.GetFieldData<FieldInfo::force> ();
  for(sArray1d::size_type i =0; i < m_boreholeSetNames.size(); ++i)
  {
    std::string setName = m_boreholeSetNames[i];
    std::map< std::string, lSet >::const_iterator setMap = faceManager.m_Sets.find( setName );
    if( setMap != faceManager.m_Sets.end() )
    {
      const lSet& set = setMap->second;

      for ( lSet::const_iterator k=set.begin() ; k!=set.end() ; ++k )
      {
        if (domain.m_feFaceManager.m_toElementsRelation[*k].size() == 1)
        {
          const int numNodesOnFace = faceManager.m_toNodesRelation[*k].size();
          const realT area = faceManager.SurfaceArea( nodeManager, *k );

          R1Tensor traction;
          faceManager.FaceNormal(nodeManager, *k, traction);
          traction *= -pressure;
          traction *= area / numNodesOnFace;

          for( lArray1d::const_iterator a=faceManager.m_toNodesRelation[*k].begin() ;
              a!=faceManager.m_toNodesRelation[*k].end() ; ++a )
          {
            force[*a] += traction;
          }
        }
      }
    }
  }
}



/// Register solver in the solver factory
REGISTER_SOLVER( HydroStaticParallelPlateFlowSolver )
