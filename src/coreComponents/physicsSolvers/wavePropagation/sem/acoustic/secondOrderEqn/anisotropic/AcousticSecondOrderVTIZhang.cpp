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

#include "AcousticSecondOrderVTIZhang.hpp"

// For the kernel factory used in applyStiffnessKernels
#include "AcousticVTIZhangWaveEquationSEMKernel.hpp"
#include "constitutive/NullModel.hpp"

// For the refactored matrix initialization
#include "DampingMatrixComputers.hpp"

// Field includes
#include "mesh/MeshFields.hpp"
#include "physicsSolvers/wavePropagation/sem/acoustic/shared/AcousticFields.hpp"
#include "AcousticVTIFields.hpp"

namespace geos
{

using namespace fields;

AcousticSecondOrderVTIZhang::AcousticSecondOrderVTIZhang( const std::string & name, Group * const parent )
  : AcousticSecondOrderAnisotropicVTISEM( name, parent )
{}

void AcousticSecondOrderVTIZhang::registerDataOnMesh( Group & meshBodies )
{
  AcousticSecondOrderAnisotropicVTISEM::registerDataOnMesh( meshBodies );

  // Zhang VTI doesn't need to register any additional fields beyond what the parent registers
  // (no equivalent to AcousticSigma like Fletcher has)
}

void AcousticSecondOrderVTIZhang::initializeMatrices( MeshLevel & mesh, string_array const & regionNames )
{
  // Use the base template which does: common VTI stuff (mass + VTI DOF arrays + Zhang damping)
  initializeMatricesTemplate< ZhangDampingComputer >( mesh, regionNames );
}

void AcousticSecondOrderVTIZhang::applyStiffnessKernels( const real64 & dt, MeshLevel & mesh, const string_array & regionNames )
{
  auto kernelFactory = acousticVTIZhangWaveEquationSEMKernels::ExplicitAcousticVTIZhangSEMFactory( dt );

  finiteElement::regionBasedKernelApplication< EXEC_POLICY,
                                               constitutive::NullModel,
                                               CellElementSubRegion >( mesh,
                                                                       regionNames,
                                                                       getDiscretizationName(),
                                                                       "",
                                                                       kernelFactory );
}

REGISTER_CATALOG_ENTRY( PhysicsSolverBase, AcousticSecondOrderVTIZhang, string const &, dataRepository::Group * const )

} // namespace geos
