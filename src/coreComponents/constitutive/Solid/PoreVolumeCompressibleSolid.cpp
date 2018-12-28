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
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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


PoreVolumeCompressibleSolid::PoreVolumeCompressibleSolid( std::string const & name, ManagedGroup * const parent ):
  ConstitutiveBase( name, parent ),
  m_poreVolumeRelation( ExponentApproximationType::Linear )
{
  RegisterViewWrapper( viewKeys.compressibility.Key(), &m_compressibility, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Solid compressibility");

  RegisterViewWrapper( viewKeys.referencePressure.Key(), &m_referencePressure, false )->
    setDefaultValue(0.0)->
    setDescription("Reference pressure for fluid compressibility");

  RegisterViewWrapper( viewKeyStruct::poreVolumeMultiplierString, &m_poreVolumeMultiplier, false );
  RegisterViewWrapper( viewKeyStruct::dPVMult_dPresString, &m_dPVMult_dPressure, false );
}

PoreVolumeCompressibleSolid::~PoreVolumeCompressibleSolid() = default;

std::unique_ptr<ConstitutiveBase>
PoreVolumeCompressibleSolid::DeliverClone( string const & name,
                                           ManagedGroup * const parent ) const
{
  std::unique_ptr<PoreVolumeCompressibleSolid> newConstitutiveRelation =
    std::make_unique<PoreVolumeCompressibleSolid>( name, parent );

  newConstitutiveRelation->m_compressibility   = this->m_compressibility;
  newConstitutiveRelation->m_referencePressure  = this->m_referencePressure;

  newConstitutiveRelation->m_poreVolumeRelation  = this->m_poreVolumeRelation;

  std::unique_ptr<ConstitutiveBase> rval = std::move( newConstitutiveRelation );

  return rval;
}

void PoreVolumeCompressibleSolid::AllocateConstitutiveData( dataRepository::ManagedGroup * const parent,
                                                            localIndex const numConstitutivePointsPerParentIndex )
{
  ConstitutiveBase::AllocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  this->resize( parent->size() );

  m_poreVolumeMultiplier.resize( parent->size(), numConstitutivePointsPerParentIndex );
  m_dPVMult_dPressure.resize( parent->size(), numConstitutivePointsPerParentIndex );
  m_poreVolumeMultiplier = 1.0;
}

void PoreVolumeCompressibleSolid::ProcessInputFile_PostProcess()
{
  if( m_compressibility < 0.0 )
  {
    string const message = "An invalid value of fluid bulk modulus (" + std::to_string(m_compressibility) + ") is specified";
    GEOS_ERROR(message);
  }
}

void PoreVolumeCompressibleSolid::FinalInitializationPreSubGroups( ManagedGroup *const parent )
{
  m_poreVolumeRelation.SetCoefficients( m_referencePressure, 1.0, m_compressibility );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, PoreVolumeCompressibleSolid, std::string const &, ManagedGroup * const )
}
} /* namespace geosx */
