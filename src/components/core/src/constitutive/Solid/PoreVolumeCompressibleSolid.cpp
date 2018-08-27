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
  RegisterViewWrapper( viewKeys.compressibility.Key(), &m_compressibility, 0 );
  RegisterViewWrapper( viewKeys.referencePressure.Key(), &m_referencePressure, 0 );

  RegisterViewWrapper( viewKeyStruct::poreVolumeMultiplierString, &m_poreVolumeMultiplier, 0 );
  RegisterViewWrapper( viewKeyStruct::dPVMult_dPresString, &m_dPVMult_dPressure, 0 );

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

void PoreVolumeCompressibleSolid::FillDocumentationNode()
{

  DocumentationNode * const docNode = this->getDocumentationNode();

  docNode->setName( this->CatalogName());
  docNode->setSchemaType( "Node" );
  docNode->setShortDescription( "Slightly compressible single phase fluid equation of state" );

  docNode->AllocateChildNode( viewKeys.compressibility.Key(),
                              viewKeys.compressibility.Key(),
                              -1,
                              "real64",
                              "real64",
                              "Fluid Bulk Modulus",
                              "Fluid Bulk Modulus",
                              "-1",
                              "",
                              1,
                              1,
                              0 );

  docNode->AllocateChildNode( viewKeys.referencePressure.Key(),
                              viewKeys.referencePressure.Key(),
                              -1,
                              "real64",
                              "real64",
                              "Reference pressure",
                              "Reference pressure",
                              "0",
                              "",
                              1,
                              1,
                              0 );
}

void PoreVolumeCompressibleSolid::ReadXML_PostProcess()
{
  if( m_compressibility <= 0.0 )
  {
    string const message = "An invalid value of fluid bulk modulus (" + std::to_string(m_compressibility) + ") is specified";
    GEOS_ERROR(message);
  }
}

void PoreVolumeCompressibleSolid::PoreVolumeMultiplierCompute(real64 const & pres,
                                                              localIndex const i,
                                                              real64 & poro,
                                                              real64 & dPVMult_dPres)
{
  m_poreVolumeRelation.Compute( pres, poro, dPVMult_dPres );
}

void PoreVolumeCompressibleSolid::FinalInitialization( ManagedGroup *const parent )
{
  m_poreVolumeRelation.SetCoefficients( m_referencePressure, 1.0, m_compressibility );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, PoreVolumeCompressibleSolid, std::string const &, ManagedGroup * const )
}
} /* namespace geosx */
