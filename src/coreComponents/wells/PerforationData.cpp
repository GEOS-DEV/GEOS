/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/*
 * @file PerforationData.cpp
 *
 */

#include "PerforationData.hpp"

#include "mpiCommunications/MpiWrapper.hpp"
#include "Perforation.hpp"

#include "managers/DomainPartition.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"

namespace geosx
{

using namespace dataRepository;

PerforationData::PerforationData(string const & name, Group * const parent)
  : ObjectManagerBase(name, parent),
    m_numPerforationsGlobal(0)
{
  registerWrapper( viewKeyStruct::numPerforationsGlobalString, &m_numPerforationsGlobal, false );
 
  registerWrapper( viewKeyStruct::reservoirElementRegionString,    &m_toMeshElements.m_toElementRegion,    false );
  registerWrapper( viewKeyStruct::reservoirElementSubregionString, &m_toMeshElements.m_toElementSubRegion, false );
  registerWrapper( viewKeyStruct::reservoirElementIndexString,     &m_toMeshElements.m_toElementIndex,     false );

  registerWrapper( viewKeyStruct::wellElementIndexString, &m_wellElementIndex, false );
  registerWrapper( viewKeyStruct::locationString, &m_location, false );
  registerWrapper( viewKeyStruct::wellPeacemanIndexString, &m_wellPeacemanIndex, false );
}

PerforationData::~PerforationData()
{
}

  
void PerforationData::ConnectToWellElements( InternalWellGenerator const & wellGeometry,
                                             unordered_map<globalIndex,localIndex> const & globalToLocalWellElemMap, 
                                             globalIndex elemOffsetGlobal )
{
  arrayView1d<globalIndex const> const & perfElemIndexGlobal = wellGeometry.GetPerfElemIndex();

  for (localIndex iperfLocal = 0; iperfLocal < size(); ++iperfLocal)
  {
    globalIndex const iwelemGlobal = perfElemIndexGlobal[m_localToGlobalMap[iperfLocal]];
    globalIndex const ielemGlobal  = elemOffsetGlobal + iwelemGlobal;
    GEOS_ASSERT( globalToLocalWellElemMap.count( ielemGlobal ) > 0 );
    m_wellElementIndex[iperfLocal] = globalToLocalWellElemMap.at( ielemGlobal ); 
  }

  //DebugLocalPerforations();
}
 
void PerforationData::InitializePostInitialConditions_PreSubGroups( Group * const GEOSX_UNUSED_ARG( problemManager ) )
{
  for (localIndex iperf = 0; iperf < size(); ++iperf)
  {
    if (m_wellPeacemanIndex[iperf] < 0.0)
    {
      // TODO: compute wellPeacemanIndex internally
      GEOS_ERROR( "Invalid well Peaceman index value: " << m_wellPeacemanIndex[iperf] );
    }
  }
}


void PerforationData::DebugLocalPerforations() const
{
  if (size() == 0)
  {
    return;
  } 

  if ( MpiWrapper::Comm_rank( MPI_COMM_GEOSX ) != 1)
  {
    return;
  } 

  std::cout << std::endl;
  std::cout << "++++++++++++++++++++++++++" << std::endl;
  std::cout << "PerforationData = " << getName() << " of " << getParent()->getName() << std::endl;
  std::cout << "MPI rank = " << MpiWrapper::Comm_rank( MPI_COMM_GEOSX ) << std::endl;
  std::cout << "Number of local perforations = " << size() << std::endl;

  for (localIndex iperf = 0; iperf < size(); ++iperf) 
  {
    std::cout << "m_toMeshElements.m_toElementRegion[" << iperf << "] = "
              << m_toMeshElements.m_toElementRegion[iperf]
              << std::endl;
    std::cout << "m_toMeshElements.m_toElementSubRegion[" << iperf << "] = "
              << m_toMeshElements.m_toElementSubRegion[iperf]
              << std::endl;
    std::cout << "m_toMeshElements.m_toElementIndex[" << iperf << "] = "
              << m_toMeshElements.m_toElementIndex[iperf]
              << std::endl;

    std::cout << "m_wellElementIndexLocal[" << iperf << "] = " << m_wellElementIndex[iperf] 
              << std::endl;
    std::cout << "m_location[" << iperf << "] = " << m_location[iperf]
              << std::endl;
    std::cout << "m_wellPeacemanIndex[" << iperf << "] = " << m_wellPeacemanIndex[iperf]
              << std::endl;
  }
}

} //namespace geosx
