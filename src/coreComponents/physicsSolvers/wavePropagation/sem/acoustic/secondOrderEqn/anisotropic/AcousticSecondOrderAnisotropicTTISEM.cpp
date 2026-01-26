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
#include "DampingMatrixComputers.hpp"


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

template< typename DampingComputer >
void AcousticSecondOrderAnisotropicTTISEM::initializeMatricesTemplate( MeshLevel & mesh, string_array const & regionNames )
{
  // First do all the common VTI stuff (mass + VTI DOF arrays + damping)
  AcousticSecondOrderAnisotropicWaveEquationSEM::initializeMatricesTemplate< DampingComputer >( mesh, regionNames );

  // Now add TTI-specific DOF arrays computation
  NodeManager & nodeManager = mesh.getNodeManager();
  ElementRegionManager & elemManager = mesh.getElemManager();

  arrayView2d< wsCoordType const, nodes::REFERENCE_POSITION_USD > const nodeCoords = nodeManager.getField< fields::referencePosition32 >().toViewConst();

  // TTI-specific DOF arrays
  arrayView1d< real32 > const dofTilt = nodeManager.getField< acousticttifields::AcousticDofTilt >();
  arrayView1d< real32 > const dofAzimuth = nodeManager.getField< acousticttifields::AcousticDofAzimuth >();
  dofTilt.zero();
  dofAzimuth.zero();

  // Get source coordinates for DOF arrays computation
  arrayView2d< real64 const > const sourceCoordinates = m_sourceCoordinates.toViewConst();

  elemManager.forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                              CellElementSubRegion & elementSubRegion )
  {
    finiteElement::FiniteElementBase const &
    fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );

    arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes = elementSubRegion.nodeList();
    arrayView1d< real32 const > const tti_dipx = elementSubRegion.getField< acousticttifields::AcousticDipX >();
    arrayView1d< real32 const > const tti_dipy = elementSubRegion.getField< acousticttifields::AcousticDipY >();
    arrayView1d< real32 > const dofOrder   = nodeManager.getField< acousticvtifields::AcousticDofOrder >();

    finiteElement::FiniteElementDispatchHandler< SEM_FE_TYPES >::dispatch3D( fe, [&] ( auto const finiteElement )
    {
      using FE_TYPE = TYPEOFREF( finiteElement );

      // TTI-specific DOF arrays computation
      AcousticMatricesSEM::DofArrays< FE_TYPE > kernelDebug( finiteElement );
      kernelDebug.template computeDofArraysTTI< EXEC_POLICY, ATOMIC_POLICY >( elementSubRegion.size(),
                                                                              nodeCoords,
                                                                              elemsToNodes,
                                                                              tti_dipx,
                                                                              tti_dipy,
                                                                              dofTilt,
                                                                              dofAzimuth,
                                                                              dofOrder,
                                                                              sourceCoordinates,
                                                                              m_radiusIsoAroundSource );
    } );
  } );
}

// Explicit template instantiations
template void AcousticSecondOrderAnisotropicTTISEM::initializeMatricesTemplate< ZhangDampingComputer >( MeshLevel & mesh, string_array const & regionNames );
template void AcousticSecondOrderAnisotropicTTISEM::initializeMatricesTemplate< FletcherDampingComputer >( MeshLevel & mesh, string_array const & regionNames );

} // namespace geos
