/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
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
  * @file ProppantSlurryFluid.cpp
  */

#include "ProppantSlurryFluid.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace cxx_utilities;

namespace constitutive
{

ProppantSlurryFluid::ProppantSlurryFluid( std::string const & name, Group * const parent ):
  SlurryFluidBase( name, parent )
{
  registerWrapper( viewKeyStruct::compressibilityString, &m_compressibility, false )->
    setApplyDefaultValue(0.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Fluid compressibility");

  registerWrapper( viewKeyStruct::referenceProppantDensityString, &m_referenceProppantDensity, false )->
    setApplyDefaultValue(1400.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Reference proppant density");

    registerWrapper( viewKeyStruct::referenceFluidDensityString, &m_referenceFluidDensity, false )->
    setApplyDefaultValue(1000.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Reference fluid density");

  registerWrapper( viewKeyStruct::referencePressureString, &m_referencePressure, false )->
    setApplyDefaultValue(1e5)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Reference pressure");

  registerWrapper( viewKeyStruct::referenceProppantVolumeFractionString, &m_referenceProppantVolumeFraction, false )->
    setApplyDefaultValue(0.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Reference proppant volume fraction");  

  registerWrapper( viewKeyStruct::referenceDensityString, &m_referenceDensity, false )->
    setApplyDefaultValue(1000.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Reference fluid density");

  registerWrapper( viewKeyStruct::referenceViscosityString, &m_referenceViscosity, false )->
    setApplyDefaultValue(0.001)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Reference fluid viscosity");

  registerWrapper( viewKeyStruct::maxProppantVolumeFractionString, &m_maxProppantVolumeFraction, false )->
    setApplyDefaultValue(0.6)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Maximum proppant volume fraction");  

}

ProppantSlurryFluid::~ProppantSlurryFluid() = default;

void ProppantSlurryFluid::AllocateConstitutiveData( dataRepository::Group * const parent,
                                                             localIndex const numConstitutivePointsPerParentIndex )
{
  SlurryFluidBase::AllocateConstitutiveData(parent, numConstitutivePointsPerParentIndex);

  m_density = m_referenceDensity;
  m_fluidDensity = m_referenceFluidDensity;  
  m_viscosity = m_referenceViscosity;
}

void
ProppantSlurryFluid::DeliverClone( string const & name,
                                            Group * const parent,
                                            std::unique_ptr<ConstitutiveBase> & clone ) const
{
  if( !clone )
  {
    clone = std::make_unique<ProppantSlurryFluid>( name, parent );
  }
  SlurryFluidBase::DeliverClone( name, parent, clone );
  ProppantSlurryFluid * const newConstitutiveRelation = dynamic_cast<ProppantSlurryFluid *>(clone.get());


  newConstitutiveRelation->m_compressibility      = this->m_compressibility;
  newConstitutiveRelation->m_referenceProppantDensity        = this->m_referenceProppantDensity;
  newConstitutiveRelation->m_referenceFluidDensity        = this->m_referenceFluidDensity;  
  newConstitutiveRelation->m_referencePressure    = this->m_referencePressure;
  newConstitutiveRelation->m_referenceProppantVolumeFraction    = this->m_referenceProppantVolumeFraction;  
  newConstitutiveRelation->m_referenceDensity     = this->m_referenceDensity;
  newConstitutiveRelation->m_referenceViscosity   = this->m_referenceViscosity;
  newConstitutiveRelation->m_maxProppantVolumeFraction   = this->m_maxProppantVolumeFraction;

}

void ProppantSlurryFluid::PostProcessInput()
{
  SlurryFluidBase::PostProcessInput();

  GEOS_ERROR_IF( m_compressibility < 0.0, "An invalid value of fluid compressibility ("
                                          << m_compressibility << ") is specified" );

  GEOS_ERROR_IF( m_referenceProppantDensity <= 0.0, "An invalid value of proppant density is specified" );

  GEOS_ERROR_IF( m_referenceDensity <= 0.0, "An invalid value of fluid reference density (" << m_compressibility << ") is specified" );

  GEOS_ERROR_IF( m_referenceViscosity <= 0.0, "An invalid value of fluid reference viscosity is specified" );

  GEOS_ERROR_IF( m_maxProppantVolumeFraction <= 0.0 || m_maxProppantVolumeFraction > 1.0, "An invalid value of maximum proppant volume fraction is specified" );  

  real64 dRho_dP, dRho_dC;
  real64 dVisc_dP, dVisc_dC;
  real64 dFluidRho_dP;
  
  Compute( m_referencePressure, m_referenceProppantVolumeFraction, m_referenceDensity, m_referenceFluidDensity, dRho_dP, dRho_dC, dFluidRho_dP, m_referenceViscosity, dVisc_dP, dVisc_dC );
  this->getWrapper< array2d<real64> >(viewKeyStruct::dDens_dPresString)->setDefaultValue( dRho_dP );
  this->getWrapper< array2d<real64> >(viewKeyStruct::dDens_dConcString)->setDefaultValue( dRho_dC );
  
  this->getWrapper< array2d<real64> >(viewKeyStruct::dVisc_dPresString)->setDefaultValue( dVisc_dP );
  this->getWrapper< array2d<real64> >(viewKeyStruct::dVisc_dConcString)->setDefaultValue( dVisc_dC );
  
}

void ProppantSlurryFluid::PointUpdate( real64 const & pressure, real64 const & concentration, localIndex const k, localIndex const q )
{
  Compute( pressure, concentration, m_density[k][q], m_fluidDensity[k][q], m_dDens_dPres[k][q], m_dDens_dConc[k][q], m_dFluidDens_dPres[k][q], m_viscosity[k][q], m_dVisc_dPres[k][q], m_dVisc_dConc[k][q]);

}

void ProppantSlurryFluid::Compute( real64 const & pressure,
				   real64 const & concentration,
				   real64 & density,
				   real64 & fluidDensity,
				   real64 & dDens_dPres,
				   real64 & dDens_dConc,
				   real64 & dFluidDens_dPres,
				   real64 & viscosity,
				   real64 & dVisc_dPres,
				   real64 & dVisc_dConc ) const
{

  // density

  fluidDensity = m_referenceDensity * exp(m_compressibility * (pressure - m_referencePressure));

  dFluidDens_dPres = m_compressibility * fluidDensity;
  
  density = (1.0 - concentration) * fluidDensity + concentration *   m_referenceProppantDensity;

  dDens_dPres = (1.0 - concentration) * dFluidDens_dPres;

  dDens_dConc = -fluidDensity + m_referenceProppantDensity;

  // Stokes-Einstein model

  viscosity = m_referenceViscosity * (1.0 + 2.5 * concentration);

  dVisc_dPres = 0.0;
  dVisc_dConc = m_referenceViscosity * 2.5;

}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, ProppantSlurryFluid, std::string const &, Group * const )

} /* namespace constitutive */

} /* namespace geosx */
