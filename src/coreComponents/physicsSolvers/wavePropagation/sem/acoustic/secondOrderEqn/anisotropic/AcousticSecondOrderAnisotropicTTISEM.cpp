/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A / TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "AcousticSecondOrderAnisotropicTTISEM.hpp"
#include "physicsSolvers/wavePropagation/sem/acoustic/shared/AcousticMatricesSEMKernel.hpp"
#include "finiteElement/FiniteElementDispatch.hpp"
#include "AcousticTTIFields.hpp"


namespace geos
{

using namespace fields;

void AcousticSecondOrderAnisotropicTTISEM::registerDataOnMesh( Group & meshBodies )
{
  // Call the parent to register all common anisotropic fields
  AcousticSecondOrderAnisotropicWaveEquationSEM::registerDataOnMesh( meshBodies );

  // Add only the TTI-specific fields
  forDiscretizationOnMeshTargets( meshBodies, [&]( std::string const &, MeshLevel & mesh, string_array const & )
  {
    NodeManager & nodeManager = mesh.getNodeManager();
    // Add TTI-specific nodal fields
    nodeManager.registerField< acousticttifields::AcousticDofTilt,
                               acousticttifields::AcousticDofAzimuth >( getName());

    ElementRegionManager & elemManager = mesh.getElemManager();
    elemManager.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion & subRegion )
    {
      // Add TTI-specific element fields
      subRegion.registerField< acousticttifields::AcousticDipX >( getName());
      subRegion.registerField< acousticttifields::AcousticDipY >( getName());
    } );
  } );
}


void AcousticSecondOrderAnisotropicTTISEM::postInputInitialization()
{
  AcousticSecondOrderAnisotropicWaveEquationSEM::postInputInitialization();
}


void AcousticSecondOrderAnisotropicTTISEM::initializePostInitialConditionsPreSubGroups()
{
  // Call parent for common initialization
  AcousticSecondOrderAnisotropicWaveEquationSEM::initializePostInitialConditionsPreSubGroups();

  // Zero out TTI-specific DOF arrays
  DomainPartition & domain = getGroupByPath< DomainPartition >( "/Problem/domain" );
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &, MeshLevel & mesh, string_array const & )
  {
    NodeManager & nodeManager = mesh.getNodeManager();
    arrayView1d< real32 > const dofTilt = nodeManager.getField< acousticttifields::AcousticDofTilt >();
    arrayView1d< real32 > const dofAzimuth = nodeManager.getField< acousticttifields::AcousticDofAzimuth >();
    dofTilt.zero();
    dofAzimuth.zero();
  } );
}

} // namespace geos
