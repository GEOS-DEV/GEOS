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
 * @file ConstantPlusParallelPlatesPermeability.cpp
 */

#include "ConstantPlusParallelPlatesPermeability.hpp"

#include "constitutive/permeability/PermeabilityExtrinsicData.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{


ConstantPlusParallelPlatesPermeability::ConstantPlusParallelPlatesPermeability( string const & name, Group * const parent ):
  PermeabilityBase( name, parent )
{
  registerExtrinsicData( extrinsicMeshData::permeability::dPerm_dDispJump{}, &m_dPerm_dDispJump );

  registerWrapper( viewKeyStruct::defaultConductivityString(), &m_defaultConductivity ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Value of the default conductivity C_{f,0} for the fracture." );
}

std::unique_ptr< ConstitutiveBase >
ConstantPlusParallelPlatesPermeability::deliverClone( string const & name,
                                                      Group * const parent ) const
{
  return ConstitutiveBase::deliverClone( name, parent );
}

void ConstantPlusParallelPlatesPermeability::allocateConstitutiveData( dataRepository::Group & parent,
                                                                       localIndex const numConstitutivePointsPerParentIndex )
{
  // NOTE: enforcing 1 quadrature point
  m_dPerm_dDispJump.resize( 0, 1, 3, 3 );

  PermeabilityBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
}


REGISTER_CATALOG_ENTRY( ConstitutiveBase, ConstantPlusParallelPlatesPermeability, string const &, Group * const )

}
} /* namespace geosx */
