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

/**
* @file PoreVolumeCompressibleSolid.cpp
*/

#include "PoreVolumeCompressibleSolid.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace cxx_utilities;

namespace constitutive
{


PoreVolumeCompressibleSolid::PoreVolumeCompressibleSolid( std::string const & name, Group * const parent ):
  ConstitutiveBase( name, parent )
{
  registerWrapper( viewKeys.compressibility.Key(), &m_compressibility, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Solid compressibility");

  registerWrapper( viewKeys.referencePressure.Key(), &m_referencePressure, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Reference pressure for fluid compressibility");

  registerWrapper( viewKeyStruct::poreVolumeMultiplierString, &m_poreVolumeMultiplier, false )->
    setDefaultValue(1.0);
  registerWrapper( viewKeyStruct::dPVMult_dPresString, &m_dPVMult_dPressure, false );
}

PoreVolumeCompressibleSolid::~PoreVolumeCompressibleSolid() = default;

void
PoreVolumeCompressibleSolid::DeliverClone( string const & name,
                                           Group * const parent,
                                           std::unique_ptr<ConstitutiveBase> & clone ) const
{
  std::unique_ptr<PoreVolumeCompressibleSolid> newConstitutiveRelation =
    std::make_unique<PoreVolumeCompressibleSolid>( name, parent );

  newConstitutiveRelation->m_compressibility   = this->m_compressibility;
  newConstitutiveRelation->m_referencePressure  = this->m_referencePressure;

  newConstitutiveRelation->m_poreVolumeRelation  = this->m_poreVolumeRelation;

  clone = std::move( newConstitutiveRelation );
}

void PoreVolumeCompressibleSolid::AllocateConstitutiveData( dataRepository::Group * const parent,
                                                            localIndex const numConstitutivePointsPerParentIndex )
{
  ConstitutiveBase::AllocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  this->resize( parent->size() );

  m_poreVolumeMultiplier.resize( parent->size(), numConstitutivePointsPerParentIndex );
  m_dPVMult_dPressure.resize( parent->size(), numConstitutivePointsPerParentIndex );
  m_poreVolumeMultiplier = 1.0;
}

void PoreVolumeCompressibleSolid::PostProcessInput()
{
  if( m_compressibility < 0.0 )
  {
    string const message = "An invalid value of fluid bulk modulus (" + std::to_string(m_compressibility) + ") is specified";
    GEOSX_ERROR(message);
  }
  m_poreVolumeRelation.SetCoefficients( m_referencePressure, 1.0, m_compressibility );

}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, PoreVolumeCompressibleSolid, std::string const &, Group * const )
}
} /* namespace geosx */
