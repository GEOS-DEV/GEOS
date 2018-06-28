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
{}

HydroStaticParallelPlateFlowSolver::~HydroStaticParallelPlateFlowSolver()
{}


void HydroStaticParallelPlateFlowSolver::ReadXML( TICPP::HierarchicalDataNode* const hdn  )
{
  ParallelPlateFlowSolverExplicit::ReadXML( hdn );
  m_cavityVolume = hdn->GetAttributeOrDefault<realT>("cavityVolume", 0.0);
  m_cavityMass = hdn->GetAttributeOrDefault<realT>("initialMass", m_rho_o * m_cavityVolume);
  m_boreholeSetNames = hdn->GetStringVector("boreholeSetNames");
}


void HydroStaticParallelPlateFlowSolver::PostProcess(PhysicalDomainT & domain,
                                                     SpatialPartition& partition,
                                                     const array<string>& namesOfSolverRegions)
{


  array<real64>& aperture = domain.m_feFaceManager.GetFieldData<realT>( ApertureStr );
  array<real64>& faceFluidPressure = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();
  const array<integer>& flowFaceType = domain.m_feFaceManager.GetFieldData<int>("flowFaceType");

  // Get field values of child faces from their parents.
  for( localIndex kf=0 ; kf<domain.m_feFaceManager.DataLengths() ; ++kf )
  {
    if (flowFaceType[kf] == 1)
    {
      for (localIndex i = 0 ; i < domain.m_feFaceManager.m_childIndices[kf].size() ; ++i)
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
                                                    const realT dt,
                                                    PhysicalDomainT& domain )
{
  array<real64>& faceFluidPressure = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();
  array<real64>& faceFluidDensity = domain.m_feFaceManager.GetFieldData<FieldInfo::density>();
  array<real64>& faceFluidMass = domain.m_feFaceManager.GetFieldData<FieldInfo::mass>();
  const array<real64>& fluidVolume  = domain.m_feFaceManager.GetFieldData<FieldInfo::volume>();

  const array<integer>& flowFaceType = domain.m_feFaceManager.GetFieldData<int>("flowFaceType");

  array<real64>* initialSaturatedTime = domain.m_feFaceManager.GetFieldDataPointer<realT>("initialSaturatedTime");
  const array<integer>& isGhost = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();


  int rank, size;
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

  for( array<BoundaryConditionBase*>::const_iterator bcItr=domain.m_feFaceManager.m_bcData.begin() ; bcItr!=domain.m_feFaceManager.m_bcData.end() ; ++bcItr )
  {
    // check to see if the requested field has a boundary condition applied to
    // it.
    BoundaryConditionBase* bc = *bcItr;
    if( streq( bc->GetFieldName(time), Field<FieldInfo::pressure>::Name()) )
    {
      for(localIndex i =0 ; i < bc->m_setNames.size() ; ++i)
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

  //  std::cout<<"  Mass In-Out = Net: "<<mass_input<<" - "<<-mass_output<<" =
  // "<<mass_input+mass_output<<std::endl;
}


void HydroStaticParallelPlateFlowSolver::ApplyBoreholePressure( const realT time, const realT dt, realT pressure, PhysicalDomainT& domain )
{
  FaceManagerT& faceManager = domain.m_feFaceManager;
  NodeManager& nodeManager = domain.m_feNodeManager;
  array<R1Tensor>& force = nodeManager.GetFieldData<FieldInfo::force> ();
  for(array<string>::size_type i =0 ; i < m_boreholeSetNames.size() ; ++i)
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
