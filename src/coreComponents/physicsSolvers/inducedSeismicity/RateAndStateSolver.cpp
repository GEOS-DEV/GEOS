/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file RateAndState.cpp
 */

#include "RateAndState.hpp"

#include "dataRepository/InputFlags.hpp"
#include "mainInterface/GeosxState.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "rateAndStateFields.hpp"


namespace geos
{

using namespace dataRepository;
using namespace constitutive;

RateAndState::RateAndState( const string & name,
                                Group * const parent ):
  SolverBase( name, parent ),
  m_stressSolver( nullptr )
{
  
}

void RateAndState::postInputInitialization()
{

  SolverBase::postInputInitialization();
}

RateAndState::~RateAndState()
{
  // TODO Auto-generated destructor stub
}

void RateAndState::registerDataOnMesh( Group & meshBodies )
{
  SolverBase::registerDataOnMesh( meshBodies );

  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & mesh,
                                                    arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< SurfaceElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   SurfaceElementSubRegion & subRegion )
    {
      subRegion.registerField< fields::rateAndState::slipRate >( getName() );

      subRegion.registerField< fields::rateAndState::stateVariable >( getName() );
    } );
  } );
}

real64 RateAndState::solverStep( real64 const & time_n,
                                 real64 const & dt,
                                 const int cycleNumber,
                                 DomainPartition & domain )
{
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                    MeshLevel & mesh,
                                                    arrayView1d< string const > const & regionNames )
  {
  ElementRegionManager & elemManager = mesh.getElemManager();
  elemManager.forElementSubRegions< SurfaceElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   SurfaceElementSubRegion & subRegion )
  {

    

  } ); 
  } ) ;
  

  return dt;
}

REGISTER_CATALOG_ENTRY( SolverBase, RateAndState, string const &, dataRepository::Group * const )
} // namespace geos
