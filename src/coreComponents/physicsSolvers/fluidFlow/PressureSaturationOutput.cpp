/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file PressureSaturationOutput.cpp
 */

#include "PressureSaturationOutput.hpp"

#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseExtrinsicData.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBaseExtrinsicData.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseExtrinsicData.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBase.hpp"
#include "mainInterface/ProblemManager.hpp"

namespace geosx
{

using namespace constitutive;
using namespace dataRepository;

PressureSaturationOutput::PressureSaturationOutput( const string & name,
                                                    Group * const parent ):
  TaskBase( name, parent ),
  m_flowSolverName()
{
  registerWrapper( viewKeyStruct::flowSolverNameString(), &m_flowSolverName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the flow solver" );
}

PressureSaturationOutput::~PressureSaturationOutput()
{}

void PressureSaturationOutput::postProcessInput()
{
  ProblemManager & problemManager = this->getGroupByPath< ProblemManager >( "/Problem" );
  PhysicsSolverManager & physicsSolverManager = problemManager.getPhysicsSolverManager();

  GEOSX_THROW_IF( !physicsSolverManager.hasGroup( m_flowSolverName ),
                  GEOSX_FMT( "Task {}: physics solver named {} not found",
                             getName(), m_flowSolverName ),
                  InputError );

  m_flowSolver = &physicsSolverManager.getGroup< CompositionalMultiphaseBase >( m_flowSolverName );
}

bool PressureSaturationOutput::execute( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                        real64 const GEOSX_UNUSED_PARAM( dt ),
                                        integer const cycleNumber,
                                        integer const GEOSX_UNUSED_PARAM( eventCounter ),
                                        real64 const GEOSX_UNUSED_PARAM( eventProgress ),
                                        DomainPartition & domain )
{
  m_flowSolver->forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                              MeshLevel & mesh,
                                                              arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions( regionNames, [&]( localIndex const,
                                                        ElementSubRegionBase & subRegion )
    {
      integer const myRank = MpiWrapper::commRank();

      std::ofstream pressureFile;
      std::ofstream saturationFile;
      size_t n_zero = 4;
      auto old_str = std::to_string( cycleNumber );
      auto new_str = std::string(n_zero - std::min(n_zero, old_str.length()), '0') + old_str;
      pressureFile.open ( "pressure_"+ std::to_string( myRank ) + "_" + new_str + ".txt", std::ios_base::app );
      saturationFile.open ( "saturation_"+ std::to_string( myRank ) + "_" + new_str + ".txt", std::ios_base::app );

      arrayView1d< real64 const > const pressure =
        subRegion.getExtrinsicData< extrinsicMeshData::flow::pressure >();
      arrayView2d< real64 const, compflow::USD_PHASE > const phaseVolFrac =
        subRegion.getExtrinsicData< extrinsicMeshData::flow::phaseVolumeFraction >();
      arrayView2d< real64 const > const elemCenter =
        subRegion.getReference< array2d< real64 > >( ElementSubRegionBase::viewKeyStruct::elementCenterString() );
      arrayView1d< integer const > const elemGhostRank = subRegion.ghostRank();

      forAll< serialPolicy >( subRegion.size(), [&] GEOSX_HOST_DEVICE ( localIndex const ei )
      {
        if( elemGhostRank[ei] < 0 )
        {
          pressureFile << subRegion.localToGlobalMap()[ei] << " " << elemCenter[ei][0] << " " << elemCenter[ei][1] << " " << elemCenter[ei][2] << " " << pressure[ei] << std::endl;
          saturationFile << subRegion.localToGlobalMap()[ei] << " " << elemCenter[ei][0] << " " << elemCenter[ei][1] << " " << elemCenter[ei][2] << " " << phaseVolFrac[ei][0] << " " <<
            phaseVolFrac[ei][1] << std::endl;
        }
      } );

    } );
  } );

  return false;
}

REGISTER_CATALOG_ENTRY( TaskBase,
                        PressureSaturationOutput,
                        string const &, dataRepository::Group * const )

} /* namespace geosx */
