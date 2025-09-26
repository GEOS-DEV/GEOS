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

#include "AcousticSecondOrderTTIFletcher.hpp"

// For the kernel factory used in applyStiffnessKernels
#include "AcousticTTIFletcherWaveEquationSEMKernel.hpp"
#include "constitutive/NullModel.hpp"

// For the refactored matrix initialization
#include "DampingMatrixComputers.hpp"

// Field includes
#include "mesh/MeshFields.hpp"
#include "physicsSolvers/wavePropagation/sem/acoustic/shared/AcousticFields.hpp"
#include "AcousticVTIFields.hpp"
#include "AcousticTTIFields.hpp"


namespace geos
{

using namespace fields;

AcousticSecondOrderTTIFletcher::AcousticSecondOrderTTIFletcher( const std::string & name, Group * const parent )
  : AcousticSecondOrderAnisotropicTTISEM( name, parent )
{}

void AcousticSecondOrderTTIFletcher::registerDataOnMesh( Group & meshBodies )
{
  AcousticSecondOrderAnisotropicTTISEM::registerDataOnMesh( meshBodies );

  // Now, register the single field that is specific to the Fletcher approximation.
  forDiscretizationOnMeshTargets( meshBodies, [&]( std::string const &, MeshLevel & mesh, string_array const & )
  {
    mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion & subRegion )
    {
      subRegion.registerField< acousticvtifields::AcousticSigma >( getName());
    } );
  } );
}

void AcousticSecondOrderTTIFletcher::initializeMatrices( MeshLevel & mesh, string_array const & regionNames )
{
  // Use the TTI template which does: common VTI stuff + TTI DOF arrays + Fletcher damping
  initializeMatricesTemplate< FletcherDampingComputer >( mesh, regionNames );
}

void AcousticSecondOrderTTIFletcher::applyStiffnessKernels( const real64 & dt, MeshLevel & mesh, const string_array & regionNames )
{
  auto kernelFactory = acousticTTIFletcherWaveEquationSEMKernels::ExplicitAcousticTTIFletcherSEMFactory( dt );

  finiteElement::regionBasedKernelApplication< EXEC_POLICY,
                                               constitutive::NullModel,
                                               CellElementSubRegion >( mesh,
                                                                       regionNames,
                                                                       getDiscretizationName(),
                                                                       "",
                                                                       kernelFactory );
}

REGISTER_CATALOG_ENTRY( PhysicsSolverBase, AcousticSecondOrderTTIFletcher, string const &, dataRepository::Group * const )

} // namespace geos
